//! CAD exchange: STL (ASCII and binary) read/write, a pure-Rust STEP
//! B-rep reader, NURBS evaluation, and triangle-mesh clean-up (consistent
//! outward orientation).
//!
//! * [`step_parse`]: the STEP physical file (entities, units, placements);
//! * [`surfaces`]: curves and surfaces (analytic, revolution, extrusion,
//!   B-spline) with evaluation and point inversion;
//! * [`trim`]: topology (shells, faces, loops, edges), edge sampling and
//!   the trimmed parameter-space domain of each face, with seam handling;
//! * [`tessellate`]: triangulation of the faces.

use foiley_core::Vec3;

pub mod step_parse;
pub mod surfaces;
pub mod tessellate;
pub mod trim;

pub use surfaces::{Curve, Frame, Surface};
pub use trim::{Domain, Model, ParamLoop, StepOptions};

pub type Tri = [Vec3; 3];

/// Clamped knot vector -> basis functions N_{i,p}(u) for i in 0..n.
pub fn bspline_basis(knots: &[f64], p: usize, u: f64, out: &mut Vec<f64>) {
    let n = knots.len().saturating_sub(p + 1);
    out.clear();
    out.resize(n, 0.0);
    if n == 0 {
        return;
    }
    let (lo, hi) = (knots[p], knots[n]);
    let u = u.clamp(lo, hi);
    // find span k with knots[k] <= u < knots[k+1] (last span for u == hi)
    let mut k = p;
    while k < n - 1 && u >= knots[k + 1] {
        k += 1;
    }
    // de Boor triangular scheme for the p+1 non-zero functions
    let mut nb = vec![0.0; p + 1];
    nb[0] = 1.0;
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    for j in 1..=p {
        left[j] = u - knots[k + 1 - j];
        right[j] = knots[k + j] - u;
        let mut saved = 0.0;
        for r in 0..j {
            let den = right[r + 1] + left[j - r];
            let tmp = if den.abs() > 1e-300 { nb[r] / den } else { 0.0 };
            nb[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        nb[j] = saved;
    }
    for j in 0..=p {
        let i = k - p + j;
        if i < n {
            out[i] = nb[j];
        }
    }
}

/// A (rational) B-spline surface: control net [u][v].
#[derive(Clone, Debug)]
pub struct NurbsSurface {
    pub pu: usize,
    pub pv: usize,
    pub ku: Vec<f64>,
    pub kv: Vec<f64>,
    pub ctrl: Vec<Vec<Vec3>>,
    pub w: Vec<Vec<f64>>,
}

impl NurbsSurface {
    pub fn domain(&self) -> ((f64, f64), (f64, f64)) {
        let nu = self.ctrl.len();
        let nv = self.ctrl.first().map(|r| r.len()).unwrap_or(0);
        ((self.ku[self.pu], self.ku[nu]), (self.kv[self.pv], self.kv[nv]))
    }
    pub fn eval(&self, u: f64, v: f64) -> Vec3 {
        let mut bu = Vec::new();
        let mut bv = Vec::new();
        bspline_basis(&self.ku, self.pu, u, &mut bu);
        bspline_basis(&self.kv, self.pv, v, &mut bv);
        let mut p = Vec3::ZERO;
        let mut ws = 0.0;
        for (i, a) in bu.iter().enumerate() {
            if *a == 0.0 {
                continue;
            }
            for (j, b) in bv.iter().enumerate() {
                if *b == 0.0 {
                    continue;
                }
                let w = a * b * self.w[i][j];
                p += self.ctrl[i][j] * w;
                ws += w;
            }
        }
        if ws.abs() > 1e-300 {
            p / ws
        } else {
            p
        }
    }
    /// Triangulate on a uniform parameter grid.
    pub fn tessellate(&self, nu: usize, nv: usize) -> Vec<Tri> {
        let ((u0, u1), (v0, v1)) = self.domain();
        let g: Vec<Vec<Vec3>> = (0..=nu).map(|i| (0..=nv).map(|j| self.eval(u0 + (u1 - u0) * i as f64 / nu as f64, v0 + (v1 - v0) * j as f64 / nv as f64)).collect()).collect();
        let mut t = Vec::new();
        for i in 0..nu {
            for j in 0..nv {
                t.push([g[i][j], g[i + 1][j], g[i + 1][j + 1]]);
                t.push([g[i][j], g[i + 1][j + 1], g[i][j + 1]]);
            }
        }
        t
    }
}

/// Clamped uniform knot vector for n control points of degree p.
pub fn clamped_uniform_knots(n: usize, p: usize) -> Vec<f64> {
    let p = p.min(n.saturating_sub(1));
    let m = n + p + 1;
    let inner = n - p;
    (0..m)
        .map(|j| {
            if j <= p {
                0.0
            } else if j >= n {
                1.0
            } else {
                (j - p) as f64 / inner as f64
            }
        })
        .collect()
}

fn tri_area(t: &Tri) -> f64 {
    0.5 * (t[1] - t[0]).cross(t[2] - t[0]).norm()
}

/// Remove degenerate triangles.
pub fn clean(tris: Vec<Tri>) -> Vec<Tri> {
    let size = bbox_diag(&tris).max(1e-12);
    tris.into_iter().filter(|t| t.iter().all(|p| p.is_finite()) && tri_area(t) > 1e-14 * size * size).collect()
}

pub fn bbox(tris: &[Tri]) -> (Vec3, Vec3) {
    let mut lo = Vec3::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Vec3::new(f64::MIN, f64::MIN, f64::MIN);
    for t in tris {
        for p in t {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
    }
    (lo, hi)
}

fn bbox_diag(tris: &[Tri]) -> f64 {
    if tris.is_empty() {
        return 0.0;
    }
    let (lo, hi) = bbox(tris);
    (hi - lo).norm()
}

/// Make triangle winding consistent across shared edges (BFS over the
/// edge graph) and outward (positive enclosed volume) per connected part.
pub fn orient(tris: &mut [Tri]) {
    use std::collections::HashMap;
    let tol = 1e-7 * bbox_diag(tris).max(1e-12);
    let key = |p: Vec3| ((p.x / tol).round() as i64, (p.y / tol).round() as i64, (p.z / tol).round() as i64);
    let n = tris.len();
    let mut edges: HashMap<((i64, i64, i64), (i64, i64, i64)), Vec<usize>> = HashMap::new();
    for (i, t) in tris.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (key(t[k]), key(t[(k + 1) % 3]));
            let e = if a < b { (a, b) } else { (b, a) };
            edges.entry(e).or_default().push(i);
        }
    }
    let mut seen = vec![false; n];
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut part = vec![start];
        seen[start] = true;
        let mut q = vec![start];
        while let Some(i) = q.pop() {
            for k in 0..3 {
                let (a, b) = (key(tris[i][k]), key(tris[i][(k + 1) % 3]));
                let e = if a < b { (a, b) } else { (b, a) };
                if let Some(nb) = edges.get(&e) {
                    for &j in nb {
                        if j == i || seen[j] {
                            continue;
                        }
                        // neighbour must traverse the shared edge in the opposite direction
                        let same_dir = (0..3).any(|m| key(tris[j][m]) == a && key(tris[j][(m + 1) % 3]) == b);
                        if same_dir {
                            tris[j].swap(1, 2);
                        }
                        seen[j] = true;
                        part.push(j);
                        q.push(j);
                    }
                }
            }
        }
        // outward: signed volume of the part must be positive
        let c = part.iter().fold(Vec3::ZERO, |a, &i| a + (tris[i][0] + tris[i][1] + tris[i][2]) / 3.0) / part.len() as f64;
        let vol: f64 = part.iter().map(|&i| (tris[i][0] - c).dot((tris[i][1] - c).cross(tris[i][2] - c))).sum();
        if vol < 0.0 {
            for &i in &part {
                tris[i].swap(1, 2);
            }
        }
    }
}

/// Read an STL file (binary or ASCII).
pub fn read_stl(bytes: &[u8]) -> Result<Vec<Tri>, String> {
    let is_ascii = bytes.len() >= 5 && bytes[..5].eq_ignore_ascii_case(b"solid") && {
        // binary files may also start with "solid": check the size rule
        let n = if bytes.len() >= 84 { u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize } else { usize::MAX };
        !(bytes.len() == 84 + 50 * n)
    };
    let tris = if is_ascii {
        let text = String::from_utf8_lossy(bytes);
        let mut v: Vec<Vec3> = Vec::new();
        for line in text.lines() {
            let l = line.trim_start();
            if let Some(r) = l.strip_prefix("vertex") {
                let p: Vec<f64> = r.split_whitespace().filter_map(|t| t.parse().ok()).collect();
                if p.len() != 3 {
                    return Err(format!("bad vertex line: {line}"));
                }
                v.push(Vec3::new(p[0], p[1], p[2]));
            }
        }
        if v.len() % 3 != 0 {
            return Err("vertex count is not a multiple of 3".into());
        }
        v.chunks(3).map(|c| [c[0], c[1], c[2]]).collect::<Vec<_>>()
    } else {
        if bytes.len() < 84 {
            return Err("file too short for a binary STL".into());
        }
        let n = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        if bytes.len() < 84 + 50 * n {
            return Err(format!("binary STL declares {n} triangles but the file is truncated"));
        }
        let f = |o: usize| f32::from_le_bytes([bytes[o], bytes[o + 1], bytes[o + 2], bytes[o + 3]]) as f64;
        (0..n)
            .map(|i| {
                let o = 84 + 50 * i + 12;
                [Vec3::new(f(o), f(o + 4), f(o + 8)), Vec3::new(f(o + 12), f(o + 16), f(o + 20)), Vec3::new(f(o + 24), f(o + 28), f(o + 32))]
            })
            .collect()
    };
    let t = clean(tris);
    if t.is_empty() {
        return Err("no triangles found".into());
    }
    Ok(t)
}

/// Binary STL.
pub fn write_stl(tris: &[Tri], name: &str) -> Vec<u8> {
    let mut b = Vec::with_capacity(84 + 50 * tris.len());
    let mut hdr = format!("Foiley {name}").into_bytes();
    hdr.resize(80, b' ');
    b.extend_from_slice(&hdr);
    b.extend_from_slice(&(tris.len() as u32).to_le_bytes());
    for t in tris {
        let n = (t[1] - t[0]).cross(t[2] - t[0]).normalized();
        for v in [n, t[0], t[1], t[2]] {
            for c in [v.x, v.y, v.z] {
                b.extend_from_slice(&(c as f32).to_le_bytes());
            }
        }
        b.extend_from_slice(&[0, 0]);
    }
    b
}

// ---------------------------------------------------------------- STEP

/// Read the surfaces of a STEP file as triangles (metres).
/// `density` = segments per control-polygon span.  Returns (triangles, notes).
pub fn read_step(text: &str, density: usize) -> Result<(Vec<Tri>, Vec<String>), String> {
    tessellate::legacy(text, density as f64)
}

/// STEP surfaces tessellated to roughly `target` triangles in total.
pub fn read_step_target(text: &str, target: usize) -> Result<(Vec<Tri>, Vec<String>), String> {
    let (probe, _) = tessellate::legacy(text, 1.0)?;
    let f = (target as f64 / probe.len().max(1) as f64).sqrt().clamp(0.15, 8.0);
    tessellate::legacy(text, f)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn basis_partition_of_unity() {
        let k = clamped_uniform_knots(6, 3);
        let mut b = Vec::new();
        for i in 0..=20 {
            bspline_basis(&k, 3, i as f64 / 20.0, &mut b);
            assert!((b.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        }
        bspline_basis(&k, 3, 0.0, &mut b);
        assert!((b[0] - 1.0).abs() < 1e-12);
        bspline_basis(&k, 3, 1.0, &mut b);
        assert!((b[5] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn stl_roundtrip_and_orientation() {
        // tetrahedron with one face flipped
        let (a, b, c, d) = (Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::Z);
        let mut t = vec![[a, c, b], [a, b, d], [b, c, d], [a, c, d]];
        orient(&mut t);
        let vol: f64 = t.iter().map(|t| t[0].dot(t[1].cross(t[2]))).sum::<f64>() / 6.0;
        assert!((vol - 1.0 / 6.0).abs() < 1e-12, "{vol}");
        let bytes = write_stl(&t, "tet");
        let r = read_stl(&bytes).unwrap();
        assert_eq!(r.len(), 4);
        let ascii = "solid x\nfacet normal 0 0 1\nouter loop\nvertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid\n";
        assert_eq!(read_stl(ascii.as_bytes()).unwrap().len(), 1);
        assert!(read_stl(b"garbage").is_err());
    }
    #[test]
    fn step_spline_surface() {
        let s = "ISO-10303-21;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,1.,0.));#3=CARTESIAN_POINT('',(1.,0.,0.));#4=CARTESIAN_POINT('',(1.,1.,1.));\n#5=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#2),(#3,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.);\nENDSEC;\nEND-ISO-10303-21;\n";
        let (t, _) = read_step(s, 4).unwrap();
        assert_eq!(t.len(), 32);
        let (lo, hi) = bbox(&t);
        assert!((hi.z - 1.0).abs() < 1e-12 && lo.z.abs() < 1e-12);
    }
}
