//! CAD exchange: STL (ASCII and binary) read/write, a STEP reader for
//! B-spline / NURBS surfaces and planar faces, NURBS evaluation, and
//! triangle-mesh clean-up (welding, consistent outward orientation).

use foiley_core::Vec3;

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

#[derive(Clone, Debug)]
enum Val {
    Ref(u64),
    Num(f64),
    Str(#[allow(dead_code)] String),
    Enum(String),
    List(Vec<Val>),
    Typed(#[allow(dead_code)] String, Box<Val>),
    Unset,
}

impl Val {
    fn list(&self) -> Option<&[Val]> {
        match self {
            Val::List(v) => Some(v),
            _ => None,
        }
    }
    fn num(&self) -> Option<f64> {
        match self {
            Val::Num(x) => Some(*x),
            Val::Typed(_, v) => v.num(),
            _ => None,
        }
    }
    fn rf(&self) -> Option<u64> {
        match self {
            Val::Ref(r) => Some(*r),
            _ => None,
        }
    }
}

/// One entity: simple (one part) or complex (several parts).
#[derive(Clone, Debug)]
struct Ent {
    parts: Vec<(String, Vec<Val>)>,
}

impl Ent {
    fn part(&self, name: &str) -> Option<&Vec<Val>> {
        self.parts.iter().find(|p| p.0 == name).map(|p| &p.1)
    }
}

struct Lexer<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> Lexer<'a> {
    fn ws(&mut self) {
        loop {
            while self.i < self.s.len() && self.s[self.i].is_ascii_whitespace() {
                self.i += 1;
            }
            if self.s[self.i..].starts_with(b"/*") {
                match self.s[self.i + 2..].windows(2).position(|w| w == b"*/") {
                    Some(k) => self.i += k + 4,
                    None => self.i = self.s.len(),
                }
            } else {
                break;
            }
        }
    }
    fn peek(&mut self) -> u8 {
        self.ws();
        *self.s.get(self.i).unwrap_or(&0)
    }
    fn ident(&mut self) -> String {
        self.ws();
        let st = self.i;
        while self.i < self.s.len() && (self.s[self.i].is_ascii_alphanumeric() || self.s[self.i] == b'_') {
            self.i += 1;
        }
        String::from_utf8_lossy(&self.s[st..self.i]).to_uppercase()
    }
    fn expect(&mut self, c: u8) -> Result<(), String> {
        if self.peek() == c {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("STEP: expected '{}' at byte {}", c as char, self.i))
        }
    }
    fn value(&mut self) -> Result<Val, String> {
        match self.peek() {
            b'#' => {
                self.i += 1;
                let st = self.i;
                while self.i < self.s.len() && self.s[self.i].is_ascii_digit() {
                    self.i += 1;
                }
                Ok(Val::Ref(std::str::from_utf8(&self.s[st..self.i]).unwrap_or("0").parse().unwrap_or(0)))
            }
            b'\'' => {
                self.i += 1;
                let mut out = Vec::new();
                while self.i < self.s.len() {
                    if self.s[self.i] == b'\'' {
                        if self.s.get(self.i + 1) == Some(&b'\'') {
                            out.push(b'\'');
                            self.i += 2;
                            continue;
                        }
                        self.i += 1;
                        break;
                    }
                    out.push(self.s[self.i]);
                    self.i += 1;
                }
                Ok(Val::Str(String::from_utf8_lossy(&out).to_string()))
            }
            b'.' => {
                self.i += 1;
                let st = self.i;
                while self.i < self.s.len() && self.s[self.i] != b'.' {
                    self.i += 1;
                }
                let e = String::from_utf8_lossy(&self.s[st..self.i]).to_string();
                self.i += 1;
                Ok(Val::Enum(e))
            }
            b'(' => {
                self.i += 1;
                let mut v = Vec::new();
                if self.peek() == b')' {
                    self.i += 1;
                    return Ok(Val::List(v));
                }
                loop {
                    v.push(self.value()?);
                    match self.peek() {
                        b',' => self.i += 1,
                        b')' => {
                            self.i += 1;
                            break;
                        }
                        _ => return Err(format!("STEP: bad list at byte {}", self.i)),
                    }
                }
                Ok(Val::List(v))
            }
            b'$' | b'*' => {
                self.i += 1;
                Ok(Val::Unset)
            }
            c if c == b'-' || c == b'+' || c.is_ascii_digit() => {
                let st = self.i;
                self.i += 1;
                while self.i < self.s.len() && (self.s[self.i].is_ascii_digit() || matches!(self.s[self.i], b'.' | b'E' | b'e' | b'-' | b'+')) {
                    self.i += 1;
                }
                let t = std::str::from_utf8(&self.s[st..self.i]).unwrap_or("0");
                let t = if t.ends_with('.') { format!("{t}0") } else { t.replace(".E", ".0E").replace(".e", ".0e") };
                t.parse::<f64>().map(Val::Num).map_err(|_| format!("STEP: bad number '{t}'"))
            }
            c if c.is_ascii_alphabetic() => {
                let name = self.ident();
                self.expect(b'(')?;
                self.i -= 1;
                let v = self.value()?;
                // typed parameter: NAME(value)
                let inner = match v {
                    Val::List(mut l) if l.len() == 1 => l.pop().unwrap(),
                    other => other,
                };
                Ok(Val::Typed(name, Box::new(inner)))
            }
            _ => Err(format!("STEP: unexpected character at byte {}", self.i)),
        }
    }
    fn args(&mut self) -> Result<Vec<Val>, String> {
        match self.value()? {
            Val::List(v) => Ok(v),
            _ => Err("STEP: expected an argument list".into()),
        }
    }
}

fn parse_step(text: &str) -> Result<std::collections::HashMap<u64, Ent>, String> {
    let start = text.find("DATA;").ok_or("not a STEP file (no DATA section)")? + 5;
    let mut lx = Lexer { s: text.as_bytes(), i: start };
    let mut map = std::collections::HashMap::new();
    loop {
        lx.ws();
        if lx.i >= lx.s.len() || lx.s[lx.i..].starts_with(b"ENDSEC") {
            break;
        }
        if lx.peek() != b'#' {
            // skip to next statement
            while lx.i < lx.s.len() && lx.s[lx.i] != b';' {
                lx.i += 1;
            }
            lx.i += 1;
            continue;
        }
        let id = match lx.value()? {
            Val::Ref(r) => r,
            _ => unreachable!(),
        };
        lx.expect(b'=')?;
        let mut parts = Vec::new();
        if lx.peek() == b'(' {
            lx.i += 1;
            while lx.peek() != b')' {
                let name = lx.ident();
                if name.is_empty() {
                    return Err(format!("STEP: bad complex entity #{id}"));
                }
                parts.push((name, lx.args()?));
            }
            lx.i += 1;
        } else {
            let name = lx.ident();
            parts.push((name, lx.args()?));
        }
        lx.expect(b';')?;
        map.insert(id, Ent { parts });
    }
    Ok(map)
}

fn length_unit_scale(map: &std::collections::HashMap<u64, Ent>) -> f64 {
    fn unit_scale(map: &std::collections::HashMap<u64, Ent>, id: u64, depth: usize) -> Option<f64> {
        if depth > 6 {
            return None;
        }
        let e = map.get(&id)?;
        e.part("LENGTH_UNIT")?;
        if let Some(si) = e.part("SI_UNIT") {
            let pre = match si.first() {
                Some(Val::Enum(p)) => match p.as_str() {
                    "MILLI" => 1e-3,
                    "CENTI" => 1e-2,
                    "DECI" => 1e-1,
                    "KILO" => 1e3,
                    "MICRO" => 1e-6,
                    _ => 1.0,
                },
                _ => 1.0,
            };
            return Some(pre);
        }
        if let Some(cb) = e.part("CONVERSION_BASED_UNIT") {
            let m = map.get(&cb.get(1)?.rf()?)?;
            let a = m.part("LENGTH_MEASURE_WITH_UNIT").or(m.part("MEASURE_WITH_UNIT"))?;
            let v = a.first()?.num()?;
            let s = unit_scale(map, a.get(1)?.rf()?, depth + 1).unwrap_or(1.0);
            return Some(v * s);
        }
        None
    }
    for e in map.values() {
        if let Some(g) = e.part("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
            if let Some(l) = g.first().and_then(|v| v.list()) {
                for u in l {
                    if let Some(s) = u.rf().and_then(|r| unit_scale(map, r, 0)) {
                        return s;
                    }
                }
            }
        }
    }
    1.0
}

/// Read the surfaces of a STEP file as triangles (metres).
/// `density` = segments per control-polygon span.  Returns (triangles, notes).
pub fn read_step(text: &str, density: usize) -> Result<(Vec<Tri>, Vec<String>), String> {
    read_step_scaled(text, density as f64)
}

/// STEP surfaces tessellated to roughly `target` triangles in total.
pub fn read_step_target(text: &str, target: usize) -> Result<(Vec<Tri>, Vec<String>), String> {
    let (probe, _) = read_step_scaled(text, 1.0)?;
    let f = (target as f64 / probe.len().max(1) as f64).sqrt().clamp(0.15, 8.0);
    read_step_scaled(text, f)
}

fn read_step_scaled(text: &str, density: f64) -> Result<(Vec<Tri>, Vec<String>), String> {
    let map = parse_step(text)?;
    let scale = length_unit_scale(&map);
    let point = |r: u64| -> Option<Vec3> {
        let e = map.get(&r)?;
        let a = e.part("CARTESIAN_POINT")?;
        let c = a.get(1)?.list()?;
        Some(Vec3::new(c.first()?.num()?, c.get(1).and_then(|v| v.num()).unwrap_or(0.0), c.get(2).and_then(|v| v.num()).unwrap_or(0.0)) * scale)
    };
    let mut tris = Vec::new();
    let mut notes = Vec::new();
    let mut n_surf = 0;
    let mut skipped: std::collections::BTreeMap<String, usize> = Default::default();
    let mut ids: Vec<&u64> = map.keys().collect();
    ids.sort();
    for id in ids {
        let e = &map[id];
        if let Some(a) = e.part("B_SPLINE_SURFACE_WITH_KNOTS") {
            // simple form: name, pu, pv, ctrl, form, uc, vc, si, umult, vmult, uknots, vknots, spec
            // complex form: B_SPLINE_SURFACE(pu, pv, ctrl, ...) + B_SPLINE_SURFACE_WITH_KNOTS(umult, vmult, uknots, vknots, spec)
            let (pu, pv, net, mults) = if let Some(b) = e.part("B_SPLINE_SURFACE") {
                (b.first().and_then(|v| v.num()), b.get(1).and_then(|v| v.num()), b.get(2), &a[..])
            } else {
                (a.get(1).and_then(|v| v.num()), a.get(2).and_then(|v| v.num()), a.get(3), &a[8.min(a.len())..])
            };
            let (Some(pu), Some(pv), Some(net)) = (pu, pv, net) else { continue };
            let rows = net.list().ok_or("bad control net")?;
            let ctrl: Option<Vec<Vec<Vec3>>> = rows.iter().map(|r| r.list().and_then(|c| c.iter().map(|p| p.rf().and_then(point)).collect())).collect();
            let Some(ctrl) = ctrl else { continue };
            let expand = |m: Option<&Val>, k: Option<&Val>| -> Option<Vec<f64>> {
                let m = m?.list()?;
                let k = k?.list()?;
                let mut out = Vec::new();
                for (mm, kk) in m.iter().zip(k) {
                    for _ in 0..mm.num()? as usize {
                        out.push(kk.num()?);
                    }
                }
                Some(out)
            };
            let (Some(ku), Some(kv)) = (expand(mults.first(), mults.get(2)), expand(mults.get(1), mults.get(3))) else { continue };
            let nu = ctrl.len();
            let nv = ctrl.first().map(|r| r.len()).unwrap_or(0);
            if ku.len() != nu + pu as usize + 1 || kv.len() != nv + pv as usize + 1 {
                notes.push(format!("surface #{id}: inconsistent knot vectors, skipped"));
                continue;
            }
            let w = match e.part("RATIONAL_B_SPLINE_SURFACE") {
                Some(r) => r.first().and_then(|l| l.list()).map(|rows| rows.iter().map(|row| row.list().map(|c| c.iter().map(|x| x.num().unwrap_or(1.0)).collect()).unwrap_or_default()).collect()).unwrap_or_else(|| vec![vec![1.0; nv]; nu]),
                None => vec![vec![1.0; nv]; nu],
            };
            let s = NurbsSurface { pu: pu as usize, pv: pv as usize, ku, kv, ctrl, w };
            let seg = |n: usize| (((n.max(2) - 1) as f64 * density).round() as usize).max(2);
            tris.extend(s.tessellate(seg(nu), seg(nv)));
            n_surf += 1;
        } else if let Some(a) = e.part("ADVANCED_FACE").or(e.part("FACE_SURFACE")) {
            // planar faces bounded by straight edges
            let surf = a.get(2).and_then(|v| v.rf()).and_then(|r| map.get(&r));
            let Some(surf) = surf else { continue };
            if surf.part("PLANE").is_none() {
                if surf.part("B_SPLINE_SURFACE_WITH_KNOTS").is_none() {
                    *skipped.entry(surf.parts[0].0.clone()).or_default() += 1;
                }
                continue;
            }
            let Some(bounds) = a.get(1).and_then(|v| v.list()) else { continue };
            for b in bounds {
                let Some(be) = b.rf().and_then(|r| map.get(&r)) else { continue };
                if be.part("FACE_OUTER_BOUND").is_none() && bounds.len() > 1 {
                    continue; // holes are ignored
                }
                let lp = be.parts[0].1.get(1).and_then(|v| v.rf()).and_then(|r| map.get(&r));
                let Some(el) = lp.and_then(|l| l.part("EDGE_LOOP")) else { continue };
                let mut poly = Vec::new();
                for oe in el.get(1).and_then(|v| v.list()).unwrap_or(&[]) {
                    let o = oe.rf().and_then(|r| map.get(&r)).and_then(|e| e.part("ORIENTED_EDGE"));
                    let Some(o) = o else { continue };
                    let fwd = !matches!(o.get(4), Some(Val::Enum(s)) if s == "F");
                    let ec = o.get(3).and_then(|v| v.rf()).and_then(|r| map.get(&r)).and_then(|e| e.part("EDGE_CURVE"));
                    let Some(ec) = ec else { continue };
                    let vtx = |k: usize| ec.get(k).and_then(|v| v.rf()).and_then(|r| map.get(&r)).and_then(|e| e.part("VERTEX_POINT")).and_then(|vp| vp.get(1)?.rf()).and_then(point);
                    let p = if fwd { vtx(1) } else { vtx(2) };
                    if let Some(p) = p {
                        poly.push(p);
                    }
                }
                for k in 1..poly.len().saturating_sub(1) {
                    tris.push([poly[0], poly[k], poly[k + 1]]);
                }
                n_surf += 1;
            }
        }
    }
    for (k, n) in skipped {
        notes.push(format!("{n} face(s) on {k} surfaces are not supported and were skipped"));
    }
    let mut t = clean(tris);
    if t.is_empty() {
        return Err("no supported surfaces found (B-spline/NURBS surfaces and planar faces are supported)".into());
    }
    orient(&mut t);
    notes.insert(0, format!("{n_surf} surface(s), unit scale {scale} m, {} triangles", t.len()));
    Ok((t, notes))
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
