//! Import of the legacy binary files of xflr5 v6 and flow5:
//!
//! * xflr5 v6 project files (`.xfl`, archive formats 200001 / 200002);
//! * flow5 project files (`.fl5`, archive formats 500001 ..= 501000);
//! * xflr5 / flow5 2D polar files (`.plr`: xflr5 "v1.99+" format 100000+ and
//!   the flow5 v7 format 500001).
//!
//! All of them are Qt `QDataStream`s.  The project files use the default
//! stream settings (big-endian, Qt >= 4.6 version, so `float`s are written as
//! 8-byte doubles); the `.plr` files are little-endian with the Qt 4.5
//! version (4-byte floats).  QString = u32 byte length (0xFFFFFFFF = null) +
//! UTF-16 in the stream's byte order; bool = 1 byte; int = i32; QColor = 11
//! bytes.
//!
//! The readers follow the C++ loading code of xflr5 v6.61 (`MainFrame::
//! serializeProjectXFL`, `Plane::serializePlaneXFL`, ...) and flow5
//! (`FileIO::serializeProjectFl5`, `serial::serialize*Fl5`) statement by
//! statement, so that every object preceding the ones we keep (plane polars,
//! plane operating points, foil operating points with boundary-layer data,
//! splines, STL planes, CAD fuselages...) is parsed and skipped exactly.
//! Reading stops once the foils, foil polars and planes have been read (the
//! trailing operating points, boats and sails are never needed).
//!
//! The reader is bounds-checked everywhere: malformed input gives an `Err`
//! naming the byte offset and the item being read, never a panic, and element
//! counts are checked against the remaining data before anything is allocated.

use std::fmt::Display;
use std::path::Path;

use foiley_core::{Spacing, Vec3};
use foiley_foil::{Foil, FoilPoint, FoilPolar};
use foiley_plane::geom::{FuseFrame, NurbsBody, NurbsFrame};
use foiley_plane::{Fuse, Plane, PointMass, Wing, WingKind, WingSection};

/// The objects read from a legacy file.
#[derive(Clone, Debug, Default)]
pub struct LegacyImport {
    pub foils: Vec<Foil>,
    pub polars: Vec<FoilPolar>,
    pub planes: Vec<Plane>,
    /// Non-fatal remarks: skipped objects, approximations, ignored settings.
    pub warnings: Vec<String>,
}

type Res<T> = Result<T, String>;

/// Largest string / array length accepted.
const MAX_LEN: usize = 50_000_000;

// ---------------------------------------------------------------------------
// bounds-checked QDataStream reader
// ---------------------------------------------------------------------------

struct Rd<'a> {
    b: &'a [u8],
    p: usize,
    le: bool,
    /// `float`s are serialized as 8-byte doubles (Qt >= 4.6 default precision)
    dbl_float: bool,
}

impl<'a> Rd<'a> {
    fn new(b: &'a [u8], le: bool, dbl_float: bool) -> Self {
        Rd { b, p: 0, le, dbl_float }
    }
    fn remaining(&self) -> usize {
        self.b.len() - self.p
    }
    fn fail<T>(&self, at: usize, what: &str, msg: impl Display) -> Res<T> {
        Err(format!("byte offset {at}: {msg} while reading {what}"))
    }
    fn take(&mut self, n: usize, what: &str) -> Res<&'a [u8]> {
        if n > self.remaining() {
            return self.fail(self.p, what, format_args!("unexpected end of data ({n} bytes needed, {} left)", self.remaining()));
        }
        let s = &self.b[self.p..self.p + n];
        self.p += n;
        Ok(s)
    }
    fn arr<const N: usize>(&mut self, what: &str) -> Res<[u8; N]> {
        let s = self.take(N, what)?;
        let mut a = [0u8; N];
        a.copy_from_slice(s);
        Ok(a)
    }
    fn i8(&mut self, what: &str) -> Res<i8> {
        Ok(self.take(1, what)?[0] as i8)
    }
    fn bool(&mut self, what: &str) -> Res<bool> {
        Ok(self.take(1, what)?[0] != 0)
    }
    fn i32(&mut self, what: &str) -> Res<i32> {
        let a = self.arr::<4>(what)?;
        Ok(if self.le { i32::from_le_bytes(a) } else { i32::from_be_bytes(a) })
    }
    fn u32(&mut self, what: &str) -> Res<u32> {
        let a = self.arr::<4>(what)?;
        Ok(if self.le { u32::from_le_bytes(a) } else { u32::from_be_bytes(a) })
    }
    fn f64(&mut self, what: &str) -> Res<f64> {
        let a = self.arr::<8>(what)?;
        Ok(if self.le { f64::from_le_bytes(a) } else { f64::from_be_bytes(a) })
    }
    /// A C++ `float`: 8 bytes in double-precision streams, else 4.
    fn flt(&mut self, what: &str) -> Res<f64> {
        if self.dbl_float {
            return self.f64(what);
        }
        let a = self.arr::<4>(what)?;
        Ok(if self.le { f32::from_le_bytes(a) } else { f32::from_be_bytes(a) } as f64)
    }
    fn flt_size(&self) -> usize {
        if self.dbl_float {
            8
        } else {
            4
        }
    }
    fn v3(&mut self, what: &str) -> Res<Vec3> {
        Ok(Vec3::new(self.f64(what)?, self.f64(what)?, self.f64(what)?))
    }
    fn v3f(&mut self, what: &str) -> Res<Vec3> {
        Ok(Vec3::new(self.flt(what)?, self.flt(what)?, self.flt(what)?))
    }
    /// QString: u32 byte count (0xFFFFFFFF = null) + UTF-16 code units.
    fn qstr(&mut self, what: &str) -> Res<String> {
        let at = self.p;
        let n = self.u32(what)?;
        if n == u32::MAX {
            return Ok(String::new());
        }
        let n = n as usize;
        if n > MAX_LEN || n % 2 != 0 {
            return self.fail(at, what, format_args!("invalid string length {n}"));
        }
        let s = self.take(n, what)?;
        let le = self.le;
        let units: Vec<u16> = s.chunks_exact(2).map(|c| if le { u16::from_le_bytes([c[0], c[1]]) } else { u16::from_be_bytes([c[0], c[1]]) }).collect();
        Ok(String::from_utf16_lossy(&units))
    }
    /// Legacy MFC string (`xfl::readString`): signed byte count + Latin-1.
    fn cstr(&mut self, what: &str) -> Res<String> {
        let n = self.i8(what)?;
        if n <= 0 {
            return Ok(String::new());
        }
        Ok(self.take(n as usize, what)?.iter().map(|&c| c as char).collect())
    }
    /// QColor / fl5Color / `xfl::readQColor`: 1 + 5 x 2 bytes.
    fn color(&mut self, what: &str) -> Res<()> {
        self.take(11, what).map(|_| ())
    }
    /// A non-negative element count whose elements take at least
    /// `min_size` bytes each (checked against the remaining data).
    fn count(&mut self, what: &str, min_size: usize) -> Res<usize> {
        let at = self.p;
        let n = self.i32(what)?;
        if n < 0 {
            return self.fail(at, what, format_args!("negative count {n}"));
        }
        let n = n as usize;
        if n > MAX_LEN || n.saturating_mul(min_size) > self.remaining() {
            return self.fail(at, what, format_args!("count {n} exceeds the data size"));
        }
        Ok(n)
    }
    /// A count read with C++ `for(i=0; i<n; i++)` semantics (negative = none),
    /// still checked against the remaining data.
    fn loose_count(&mut self, what: &str, min_size: usize) -> Res<usize> {
        let at = self.p;
        let n = self.i32(what)?.max(0) as usize;
        if n > MAX_LEN || n.saturating_mul(min_size) > self.remaining() {
            return self.fail(at, what, format_args!("count {n} exceeds the data size"));
        }
        Ok(n)
    }
    fn skip(&mut self, n: usize, what: &str) -> Res<()> {
        self.take(n, what).map(|_| ())
    }
    fn skip_i32s(&mut self, n: usize, what: &str) -> Res<()> {
        self.skip(n.saturating_mul(4), what)
    }
    fn skip_f64s(&mut self, n: usize, what: &str) -> Res<()> {
        self.skip(n.saturating_mul(8), what)
    }
    fn skip_flts(&mut self, n: usize, what: &str) -> Res<()> {
        self.skip(n.saturating_mul(self.flt_size()), what)
    }
    fn skip_bools(&mut self, n: usize, what: &str) -> Res<()> {
        self.skip(n, what)
    }
    /// Archive format number, checked against the accepted range.
    fn fmt(&mut self, what: &str, lo: i32, hi: i32) -> Res<i32> {
        let at = self.p;
        let f = self.i32(what)?;
        if f < lo || f > hi {
            return self.fail(at, what, format_args!("unsupported format number {f} (expected {lo}..={hi})"));
        }
        Ok(f)
    }
    /// flow5 "dynamic space allocation": int count + ints, int count + doubles.
    fn spares(&mut self, what: &str) -> Res<()> {
        let n = self.loose_count(what, 4)?;
        self.skip_i32s(n, what)?;
        let n = self.loose_count(what, 8)?;
        self.skip_f64s(n, what)
    }
}

// ---------------------------------------------------------------------------
// intermediate objects and conversions
// ---------------------------------------------------------------------------

struct RawFoil {
    name: String,
    desc: String,
    te_flap: bool,
    te_angle: f64,
    /// hinge x/c, y/t as fractions
    hinge: (f64, f64),
    x: Vec<f64>,
    y: Vec<f64>,
}

#[derive(Clone, Copy, Default)]
struct RawPt {
    alpha: f64,
    cd: f64,
    cdp: f64,
    cl: f64,
    cm: f64,
    xtr1: f64,
    xtr2: f64,
    cpmn: f64,
    re: f64,
    ctrl: f64,
}

struct RawPolar {
    foil: String,
    name: String,
    ptype: u8,
    re: f64,
    mach: f64,
    ncrit: f64,
    xtop: f64,
    xbot: f64,
    viscous: bool,
    te_flap: f64,
    pts: Vec<RawPt>,
}

fn push_foil(imp: &mut LegacyImport, f: RawFoil) {
    let name = if f.name.trim().is_empty() { format!("Foil {}", imp.foils.len() + 1) } else { f.name.trim().to_string() };
    let foil = Foil { name: name.clone(), x: f.x, y: f.y, hinge: if f.te_flap { Some(f.hinge) } else { None }, note: f.desc };
    if let Err(e) = foil.validate() {
        imp.warnings.push(format!("foil '{name}' skipped: {e}"));
        return;
    }
    if f.te_flap && f.te_angle.abs() > 1e-6 {
        imp.warnings.push(format!("foil '{name}': the trailing-edge flap deflection of {:.2} deg is not applied (base geometry imported, hinge kept)", f.te_angle));
    }
    // a later foil with the same name replaces the former one (as in xflr5)
    imp.foils.retain(|g| g.name != foil.name);
    imp.foils.push(foil);
}

fn push_polar(imp: &mut LegacyImport, p: RawPolar) {
    let mut polar = FoilPolar { name: p.name.trim().to_string(), foil: p.foil.trim().to_string(), polar_type: p.ptype, re: p.re, mach: p.mach, ncrit: p.ncrit, xtr_top: p.xtop, xtr_bot: p.xbot, viscous: p.viscous, ..Default::default() };
    let mut bad = 0;
    for q in &p.pts {
        let vals = [q.alpha, q.cl, q.cd, q.cdp, q.cm, q.xtr1, q.xtr2, q.cpmn, q.re, q.ctrl];
        if vals.iter().any(|v| !v.is_finite()) {
            bad += 1;
            continue;
        }
        let mach = if p.ptype == 2 && q.cl > 1e-9 { p.mach / q.cl.sqrt() } else { p.mach };
        let re = if q.re > 0.0 { q.re } else { p.re };
        polar.insert(FoilPoint {
            alpha: q.alpha,
            cl: q.cl,
            cd: q.cd,
            cdp: q.cdp,
            cdf: q.cd - q.cdp,
            cm: q.cm,
            xtr_top: q.xtr1,
            xtr_bot: q.xtr2,
            cpmin: q.cpmn,
            re,
            mach,
            flap: if p.ptype == 6 { q.ctrl } else { p.te_flap },
            ..Default::default()
        });
    }
    if bad > 0 {
        imp.warnings.push(format!("polar '{}' / '{}': {bad} point(s) with non-finite values skipped", polar.foil, polar.name));
    }
    imp.polars.retain(|q| !(q.foil == polar.foil && q.name == polar.name));
    imp.polars.push(polar);
}

fn xflr5_dist(k: i32) -> Spacing {
    match k {
        1 => Spacing::Cosine,
        2 => Spacing::MinusSine,
        -2 => Spacing::Sine,
        _ => Spacing::Uniform,
    }
}

fn flow5_dist(k: i32) -> Spacing {
    match k {
        1 => Spacing::Cosine,
        2 | 5 => Spacing::MinusSine,
        -2 | 6 => Spacing::Sine,
        3 | 4 => Spacing::Tanh,
        _ => Spacing::Uniform,
    }
}

fn panels(n: i32, lo: usize) -> usize {
    (n.max(0) as usize).max(lo)
}

/// Elliptic-frame approximation of control-point frames (x, points relative
/// to the fuse position): half width = max |y|, centre / half height from the
/// z extent.
fn ellipse_frames(frames: &[NurbsFrame]) -> Vec<FuseFrame> {
    let mut out: Vec<FuseFrame> = frames
        .iter()
        .filter(|f| !f.points.is_empty())
        .map(|f| {
            let hw = f.points.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
            let zmax = f.points.iter().map(|p| p.z).fold(f64::MIN, f64::max);
            let zmin = f.points.iter().map(|p| p.z).fold(f64::MAX, f64::min);
            FuseFrame { x: f.x, z: 0.5 * (zmax + zmin), half_width: hw, half_height: 0.5 * (zmax - zmin) }
        })
        .filter(|f| f.x.is_finite() && f.z.is_finite() && f.half_width.is_finite() && f.half_height.is_finite())
        .collect();
    out.sort_by(|a, b| a.x.total_cmp(&b.x));
    out.dedup_by(|a, b| (a.x - b.x).abs() < 1e-9);
    out
}

/// Orient frames as our NURBS bodies expect: right half (y >= 0), points
/// ordered from the top of the section to the bottom.
fn normalize_nurbs_frames(frames: &mut [NurbsFrame]) {
    let flip_y = frames.iter().flat_map(|f| f.points.iter()).map(|p| p.y).fold(0.0, f64::max) <= 0.0;
    for f in frames.iter_mut() {
        if flip_y {
            for p in &mut f.points {
                p.y = -p.y;
            }
        }
        if let (Some(a), Some(b)) = (f.points.first(), f.points.last()) {
            if a.z < b.z {
                f.points.reverse();
            }
        }
    }
}

struct FuseSpec {
    name: String,
    position: Vec3,
    frames: Vec<NurbsFrame>,
    flat: bool,
    degree_u: i32,
    degree_v: i32,
    nx: usize,
    nh: usize,
    mass: f64,
}

fn make_fuse(s: FuseSpec, warnings: &mut Vec<String>) -> Option<Fuse> {
    let mut frames: Vec<NurbsFrame> = s.frames.into_iter().filter(|f| f.points.len() >= 2 && f.x.is_finite() && f.points.iter().all(|p| p.is_finite())).collect();
    frames.sort_by(|a, b| a.x.total_cmp(&b.x));
    normalize_nurbs_frames(&mut frames);
    if frames.len() < 2 {
        warnings.push(format!("fuselage '{}' skipped: fewer than 2 usable frames", s.name));
        return None;
    }
    // all frames must have the same number of control points for a tensor-product surface
    let np = frames.iter().map(|f| f.points.len()).min().unwrap_or(0);
    if frames.iter().any(|f| f.points.len() != np) {
        warnings.push(format!("fuselage '{}': frames have different point counts; truncated to {np} points", s.name));
        for f in &mut frames {
            f.points.truncate(np);
        }
    }
    let (du, dv) = if s.flat { (1, 1) } else { (s.degree_u, s.degree_v) };
    let du = if (1..=5).contains(&du) { du as usize } else { 3 };
    let dv = if (1..=5).contains(&dv) { dv as usize } else { 3 };
    let ell = ellipse_frames(&frames);
    Some(Fuse {
        name: if s.name.trim().is_empty() { "Fuselage".into() } else { s.name.trim().to_string() },
        position: s.position,
        frames: ell,
        exponent: if s.flat { 4.0 } else { 2.0 },
        nx: s.nx.max(2),
        nh: s.nh.max(4),
        mass: s.mass,
        nurbs: Some(NurbsBody { frames, degree_u: du, degree_v: dv }),
        ..Fuse::default_pod()
    })
}

fn finish_plane(mut plane: Plane, warnings: &mut Vec<String>) -> Option<Plane> {
    if plane.name.trim().is_empty() {
        plane.name = "Imported plane".into();
    }
    if plane.wings.is_empty() {
        warnings.push(format!("plane '{}' skipped: no usable wing", plane.name));
        return None;
    }
    if let Err(e) = plane.validate() {
        warnings.push(format!("plane '{}': {e}", plane.name));
    }
    Some(plane)
}

fn empty_plane(name: String, description: String) -> Plane {
    Plane { name, description, wings: vec![], fuses: vec![], point_masses: vec![], ..Plane::default_plane("", "") }
}

fn pm(tag: String, mass: f64, position: Vec3) -> PointMass {
    PointMass { tag, mass, position }
}

// ---------------------------------------------------------------------------
// shared pieces
// ---------------------------------------------------------------------------

/// `LineStyle::serializeXfl` / `serial::serializeLineStyleXfl`.
fn style_xfl(r: &mut Rd) -> Res<()> {
    let w = "line style";
    r.skip_i32s(3, w)?;
    r.color(w)?;
    r.bool(w)?;
    Ok(())
}

/// `LineStyle::serializeFl5` (xflr5) / `serial::serializeLineStyleFl5` (flow5).
fn style_fl5(r: &mut Rd) -> Res<()> {
    let w = "line style";
    let k = r.i32(w)?;
    if k >= 500001 {
        r.i32(w)?; // stipple
    }
    r.i32(w)?; // width
    r.i32(w)?; // symbol
    r.color(w)?;
    r.bool(w)?;
    r.qstr(w)?; // tag
    Ok(())
}

// ---------------------------------------------------------------------------
// xflr5 v6 .xfl
// ---------------------------------------------------------------------------

/// Read an xflr5 v6 project file (`.xfl`).
pub fn read_xfl(bytes: &[u8]) -> Result<LegacyImport, String> {
    let mut r = Rd::new(bytes, false, true);
    let mut imp = LegacyImport::default();
    let fmt = r.fmt("xfl project header", 200001, 210000)?;
    r.skip_i32s(6, "unit indices")?;
    if fmt == 200001 {
        let w = "default plane polar";
        r.skip_i32s(2, w)?;
        r.skip_f64s(5, w)?; // mass, QInf, CoG
        r.skip_flts(2, w)?; // density, viscosity
        r.skip_f64s(2, w)?; // alpha, beta
        r.skip_bools(2, w)?;
    } else {
        xfl_wpolar(&mut r)?;
    }
    let n = r.count("number of planes", 1)?;
    for _ in 0..n {
        let (plane, w) = xfl_plane(&mut r)?;
        imp.warnings.extend(w);
        if let Some(p) = plane {
            imp.planes.push(p);
        }
    }
    let n = r.count("number of plane polars", 1)?;
    for _ in 0..n {
        xfl_wpolar(&mut r)?;
    }
    let n = r.count("number of plane operating points", 1)?;
    for _ in 0..n {
        xfl_popp(&mut r)?;
    }
    let n = r.count("number of foils", 1)?;
    for _ in 0..n {
        let f = xfl_foil(&mut r)?;
        push_foil(&mut imp, f);
    }
    let n = r.count("number of foil polars", 1)?;
    for _ in 0..n {
        let p = xfl_polar(&mut r)?;
        push_polar(&mut imp, p);
    }
    // foil operating points, the spline foil and the unit provisions follow:
    // nothing we need.
    Ok(imp)
}

/// `MainFrame::serializeFoilXFL`.
fn xfl_foil(r: &mut Rd) -> Res<RawFoil> {
    let w = "xfl foil";
    let fmt = r.fmt(w, 100000, 110000)?;
    let name = r.qstr(w)?;
    let desc = r.qstr(w)?;
    if fmt < 100007 {
        r.skip_i32s(2, w)?;
        r.color(w)?;
        r.bool(w)?;
        r.i8(w)?;
    } else {
        style_xfl(r)?;
    }
    r.skip_bools(2, w)?; // centre line, LE flap
    let te_flap = r.bool(w)?;
    r.skip_f64s(3, w)?; // LE flap angle, hinge
    let te_angle = r.f64(w)?;
    let xh = r.f64(w)?;
    let yh = r.f64(w)?;
    let n = r.count("xfl foil point count", 16)?;
    let (mut x, mut y) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        x.push(r.f64(w)?);
        y.push(r.f64(w)?);
    }
    Ok(RawFoil { name, desc, te_flap, te_angle, hinge: (xh / 100.0, yh / 100.0), x, y })
}

/// `MainFrame::serializePolarXFL`.
fn xfl_polar(r: &mut Rd) -> Res<RawPolar> {
    let w = "xfl foil polar";
    let fmt = r.fmt(w, 100000, 110000)?;
    let foil = r.qstr(w)?;
    let name = r.qstr(w)?;
    if fmt < 100005 {
        r.skip_i32s(2, w)?;
        r.color(w)?;
        r.skip_bools(2, w)?;
    } else {
        style_xfl(r)?;
    }
    let ptype = match r.i32(w)? {
        2 => 2,
        3 => 3,
        4 => 4,
        5 => 6, // flow5-written control polars
        _ => 1,
    };
    r.skip_i32s(2, w)?; // MaType, ReType
    let re = r.f64(w)?;
    let mach = r.f64(w)?;
    r.f64(w)?; // alpha spec
    let xtop = r.f64(w)?;
    let xbot = r.f64(w)?;
    let ncrit = r.f64(w)?;
    let n = r.count("xfl polar point count", 11 * r.flt_size())?;
    let mut pts = Vec::with_capacity(n);
    for _ in 0..n {
        let alpha = r.flt(w)?;
        let cd = r.flt(w)?;
        let cdp = r.flt(w)?;
        let cl = r.flt(w)?;
        let cm = r.flt(w)?;
        let xtr1 = r.flt(w)?;
        let xtr2 = r.flt(w)?;
        r.flt(w)?; // HMom
        let cpmn = r.flt(w)?;
        let pre = r.flt(w)?;
        r.flt(w)?; // XCp
        pts.push(RawPt { alpha, cd, cdp, cl, cm, xtr1, xtr2, cpmn, re: pre, ctrl: 0.0 });
    }
    if fmt < 100005 {
        r.i32(w)?;
    }
    r.skip_i32s(19, w)?;
    r.skip_f64s(50, w)?;
    Ok(RawPolar { foil, name, ptype, re, mach, ncrit, xtop, xbot, viscous: true, te_flap: 0.0, pts })
}

/// `WPolar::serializeWPlrXFL` (parsed and discarded).
fn xfl_wpolar(r: &mut Rd) -> Res<()> {
    let w = "xfl plane polar";
    let fmt = r.fmt(w, 200000, 205000)?;
    r.qstr(w)?;
    r.qstr(w)?;
    r.skip_f64s(3, w)?;
    if fmt < 200014 {
        r.skip_i32s(2, w)?;
        r.color(w)?;
        r.skip_bools(2, w)?;
    } else {
        style_xfl(r)?;
    }
    r.skip_i32s(2, w)?; // method, type
    r.skip_bools(7, w)?;
    r.skip_f64s(3, w)?; // height, density, viscosity
    r.i32(w)?; // reference dimensions
    r.bool(w)?; // auto inertia
    r.skip_f64s(8, w)?; // mass, CoG, inertia
    let at = r.p;
    let nc = r.loose_count("plane polar control count", 8)?;
    if nc > 100 {
        return r.fail(at, w, format_args!("{nc} controls (max 100)"));
    }
    r.skip_f64s(nc, w)?;
    r.i32(w)?; // wake panels
    r.skip_f64s(2, w)?;
    r.skip_f64s(3, w)?; // QInf, alpha, beta
    let at = r.p;
    let n = r.i32(w)?;
    if n.abs() > 10000 {
        return r.fail(at, w, format_args!("{n} data points"));
    }
    let n = n.max(0) as usize;
    r.skip_f64s(n * 36, w)?;
    r.skip_i32s(20, w)?;
    r.skip_f64s(35 + 4 + 4 + 7, w)?;
    Ok(())
}

/// `PlaneOpp::serializePOppXFL` (parsed and discarded).
fn xfl_popp(r: &mut Rd) -> Res<()> {
    let w = "xfl plane operating point";
    let fmt = r.fmt(w, 200000, 200100)?;
    r.qstr(w)?;
    r.qstr(w)?;
    if fmt < 200002 {
        r.skip_i32s(2, w)?;
        r.color(w)?;
        r.skip_bools(2, w)?;
    } else {
        style_xfl(r)?;
    }
    r.skip_bools(4, w)?;
    r.i32(w)?; // polar type
    let method = r.i32(w)?;
    let np = r.loose_count("plane operating point panel count", 0)?;
    r.i32(w)?; // NStation
    r.skip_f64s(5, w)?;
    if method != 1 {
        r.skip_flts(np.saturating_mul(3), w)?;
    }
    for _ in 0..4 {
        if r.i32(w)? != 0 {
            xfl_wingopp(r)?;
        }
    }
    r.skip_f64s(14 + 18, w)?;
    r.i32(w)?; // nControls
    r.skip_f64s(6 + 8 + 32 + 1 + 80, w)?;
    r.skip_i32s(20, w)?;
    r.skip_f64s(2 + 4 + 44, w)?;
    Ok(())
}

/// `WingOpp::serializeWingOppXFL` (parsed and discarded).
fn xfl_wingopp(r: &mut Rd) -> Res<()> {
    let w = "xfl wing operating point";
    r.i32(w)?;
    r.qstr(w)?;
    r.qstr(w)?;
    r.i32(w)?;
    r.bool(w)?;
    let ns = r.loose_count("wing operating point station count", 22 * 8)?;
    r.skip_i32s(3, w)?; // VLM panels, wake nodes, wake panels
    r.skip_f64s(2 + 6 + 3 + 2 + 3 + 2 + 3, w)?;
    r.skip_f64s(ns * 22, w)?;
    let nf = r.loose_count("wing operating point flap count", 8)?;
    r.skip_f64s(nf, w)?;
    r.skip_i32s(20, w)?;
    r.skip_f64s(50, w)?;
    Ok(())
}

struct XflWing {
    fmt: i32,
    name: String,
    symmetric: bool,
    sections: Vec<WingSection>,
    mass: f64,
    pms: Vec<PointMass>,
    ints: [i32; 20],
}

/// `Wing::serializeWingXFL`.
fn xfl_wing(r: &mut Rd) -> Res<XflWing> {
    let w = "xfl wing";
    let fmt = r.fmt(w, 100000, 100003)?;
    let name = r.qstr(w)?;
    r.qstr(w)?; // description
    r.color(w)?;
    let symmetric = r.bool(w)?;
    let n = r.count("xfl wing section count", 64)?;
    let mut sections = Vec::with_capacity(n);
    for _ in 0..n {
        let right = r.qstr(w)?;
        let left = r.qstr(w)?;
        let chord = r.f64(w)?;
        let y = r.f64(w)?;
        let offset = r.f64(w)?;
        let dihedral = r.f64(w)?;
        let twist = r.f64(w)?;
        let nx = r.i32(w)?;
        let ny = r.i32(w)?;
        let xd = r.i32(w)?;
        let yd = r.i32(w)?;
        let mut s = WingSection::new(y, chord, offset, &right);
        s.foil_left = left;
        s.dihedral = dihedral;
        s.twist = twist;
        s.nx = panels(nx, 1);
        s.ny = panels(ny, 1);
        s.xdist = xflr5_dist(xd);
        s.ydist = xflr5_dist(yd);
        sections.push(s);
    }
    let mass = r.f64(w)?;
    let n = r.count("xfl wing point mass count", 36)?;
    let mut pms = Vec::with_capacity(n);
    for _ in 0..n {
        let m = r.f64(w)?;
        let p = r.v3(w)?;
        let tag = r.qstr(w)?;
        pms.push(pm(tag, m, p));
    }
    let mut ints = [0i32; 20];
    for v in ints.iter_mut() {
        *v = r.i32(w)?;
    }
    r.skip_f64s(50, w)?;
    Ok(XflWing { fmt, name, symmetric, sections, mass, pms, ints })
}

struct XflBody {
    name: String,
    flat: bool,
    nx: i32,
    nh: i32,
    hpanels: Vec<i32>,
    xpanels: Vec<i32>,
    frames: Vec<NurbsFrame>,
    mass: f64,
    pms: Vec<PointMass>,
    udeg: i32,
    vdeg: i32,
}

/// `Body::serializeBodyXFL` with `Frame::serializeFrame`.
fn xfl_body(r: &mut Rd) -> Res<XflBody> {
    let w = "xfl body";
    r.fmt(w, 100000, 200000)?;
    let name = r.qstr(w)?;
    r.qstr(w)?;
    r.color(w)?;
    let flat = r.i32(w)? == 1;
    r.i32(w)?; // iRes
    let nx = r.i32(w)?;
    let nh = r.i32(w)?;
    r.f64(w)?; // bunch
    let n = r.count("xfl body side line count", 4)?;
    let mut hpanels = Vec::with_capacity(n);
    for _ in 0..n {
        hpanels.push(r.i32(w)?);
    }
    let n = r.count("xfl body frame count", 16)?;
    let mut frames = Vec::with_capacity(n);
    let mut xpanels = Vec::with_capacity(n);
    for _ in 0..n {
        xpanels.push(r.i32(w)?);
        let u = r.f64(w)?;
        r.fmt("xfl body frame", 1000, 1100)?;
        let np = r.count("xfl body frame point count", 3 * r.flt_size())?;
        let mut points = Vec::with_capacity(np);
        for _ in 0..np {
            let p = r.v3f(w)?;
            points.push(Vec3::new(u, p.y, p.z));
        }
        frames.push(NurbsFrame { x: u, points });
    }
    let mass = r.f64(w)?;
    let n = r.count("xfl body point mass count", 36)?;
    let mut pms = Vec::with_capacity(n);
    for _ in 0..n {
        let m = r.f64(w)?;
        let p = r.v3(w)?;
        let tag = r.qstr(w)?;
        pms.push(pm(tag, m, p));
    }
    r.skip_i32s(18, w)?;
    let udeg = r.i32(w)?;
    let vdeg = r.i32(w)?;
    r.skip_f64s(50, w)?;
    Ok(XflBody { name, flat, nx, nh, hpanels, xpanels, frames, mass, pms, udeg, vdeg })
}

/// `Plane::serializePlaneXFL`.  Returns the plane (None if unusable) and warnings.
fn xfl_plane(r: &mut Rd) -> Res<(Option<Plane>, Vec<String>)> {
    let w = "xfl plane";
    let fmt = r.fmt(w, 100001, 110000)?;
    let name = r.qstr(w)?.trim().to_string();
    let desc = r.qstr(w)?;
    if fmt >= 100002 {
        style_fl5(r)?;
    }
    let mut wings = Vec::with_capacity(4);
    for _ in 0..4 {
        wings.push(xfl_wing(r)?);
    }
    let biplane = r.bool(w)?;
    let stab = r.bool(w)?;
    let fin = r.bool(w)?;
    let double_fin = r.bool(w)?;
    let sym_fin = r.bool(w)?;
    r.bool(w)?;
    let fix = |v: f64| if v.is_nan() || v.abs() < 1e-6 || v.abs() > 1000.0 { 0.0 } else { v };
    let mut le = [(Vec3::ZERO, 0.0); 4];
    for l in le.iter_mut() {
        let p = r.v3(w)?;
        let tilt = r.f64(w)?;
        *l = (Vec3::new(fix(p.x), fix(p.y), fix(p.z)), fix(tilt));
    }
    let has_body = r.bool(w)?;
    let bx = r.f64(w)?;
    let bz = r.f64(w)?;
    let mut body = None;
    if has_body {
        let bname = r.qstr(w)?;
        let mut b = xfl_body(r)?;
        if b.name.trim().is_empty() {
            b.name = bname;
        }
        body = Some(b);
    }
    let n = r.count("xfl plane point mass count", 36)?;
    let mut pms = Vec::with_capacity(n);
    for _ in 0..n {
        let m = r.f64(w)?;
        let p = r.v3(w)?;
        let tag = r.qstr(w)?;
        pms.push(pm(tag, m, p));
    }
    r.skip_i32s(20, w)?;
    r.skip_f64s(50, w)?;

    // --- assemble ---
    let mut warnings = Vec::new();
    let mut plane = empty_plane(name, desc);
    plane.point_masses = pms;
    let present = [true, biplane, stab, fin];
    let slot_kinds = [WingKind::Main, WingKind::Other, WingKind::Elevator, WingKind::Fin];
    let slot_names = ["Main wing", "Second wing", "Elevator", "Fin"];
    for (iw, xw) in wings.into_iter().enumerate() {
        if !present[iw] {
            continue;
        }
        // plane format 100003+ (later xflr5 builds) keeps the wing's own type
        let own_type = if xw.fmt >= 100002 { xw.ints[18] } else { xw.ints[19] };
        let kind = if fmt >= 100003 {
            match own_type {
                0 => WingKind::Main,
                2 => WingKind::Elevator,
                3 => WingKind::Fin,
                1 | 4 => WingKind::Other,
                _ => slot_kinds[iw],
            }
        } else {
            slot_kinds[iw]
        };
        let wname = if xw.name.trim().is_empty() { slot_names[iw].to_string() } else { xw.name.trim().to_string() };
        if xw.sections.len() < 2 {
            warnings.push(format!("plane '{}': wing '{wname}' skipped (fewer than 2 sections)", plane.name));
            continue;
        }
        let mut wing = Wing::new(&wname, kind, xw.sections);
        wing.position = le[iw].0;
        wing.ry = le[iw].1;
        wing.symmetric = xw.symmetric;
        wing.mass = xw.mass;
        wing.point_masses = xw.pms; // relative to the wing position, as in xflr5
        let is_fin = kind == WingKind::Fin;
        if is_fin {
            wing.rx = 90.0;
            wing.two_sided = sym_fin;
            wing.close_inner = true;
        }
        if is_fin && double_fin {
            let mut a = wing.clone();
            let mut b = wing.clone();
            let y = wing.position.y.abs();
            a.position.y = -y;
            b.position.y = y;
            a.name = format!("{} (left)", wing.name);
            b.name = format!("{} (right)", wing.name);
            plane.wings.push(a);
            plane.wings.push(b);
        } else {
            plane.wings.push(wing);
        }
    }
    if let Some(b) = body {
        let bpos = Vec3::new(bx, 0.0, bz);
        for m in &b.pms {
            plane.point_masses.push(pm(m.tag.clone(), m.mass, m.position + bpos));
        }
        let (nx, nh) = if b.flat { (b.xpanels.iter().map(|&k| k.max(0) as usize).sum(), b.hpanels.iter().map(|&k| k.max(0) as usize).sum()) } else { (panels(b.nx, 0), panels(b.nh, 0)) };
        let spec = FuseSpec { name: b.name, position: bpos, frames: b.frames, flat: b.flat, degree_u: b.udeg, degree_v: b.vdeg, nx, nh, mass: b.mass };
        if let Some(f) = make_fuse(spec, &mut warnings) {
            plane.fuses.push(f);
        }
    }
    Ok((finish_plane(plane, &mut warnings), warnings))
}

// ---------------------------------------------------------------------------
// .plr polar files
// ---------------------------------------------------------------------------

/// Read a `.plr` polar file (xflr5 v1.99+ / v6 format, or flow5 v7 format).
pub fn read_plr(bytes: &[u8]) -> Result<LegacyImport, String> {
    match read_plr_stream(bytes) {
        Ok(r) => Ok(r),
        Err(e) => {
            // flow5 opens .plr files in text mode for writing: on Windows every
            // 0x0A byte then gets a 0x0D inserted before it.  Undo and retry.
            if bytes.windows(2).any(|w| w == b"\r\n") {
                let mut fixed = Vec::with_capacity(bytes.len());
                let mut i = 0;
                while i < bytes.len() {
                    if bytes[i] == b'\r' && bytes.get(i + 1) == Some(&b'\n') {
                        i += 1;
                        continue;
                    }
                    fixed.push(bytes[i]);
                    i += 1;
                }
                if let Ok(mut r) = read_plr_stream(&fixed) {
                    r.warnings.push("the file had Windows text-mode line-end translation (CR LF); it was undone".into());
                    return Ok(r);
                }
            }
            Err(e)
        }
    }
}

fn read_plr_stream(bytes: &[u8]) -> Res<LegacyImport> {
    let mut r = Rd::new(bytes, true, false);
    let mut imp = LegacyImport::default();
    let at = r.p;
    let n = r.i32("plr header")?;
    if n < 100000 {
        return r.fail(at, "plr header", format_args!("obsolete or unknown .plr format {n} (only xflr5 v1.99+ and flow5 files are supported)"));
    }
    if (500000..600000).contains(&n) {
        let nf = r.count("number of foils", 1)?;
        for _ in 0..nf {
            let f = fl5_foil(&mut r)?;
            push_foil(&mut imp, f);
        }
        let np = r.count("number of polars", 1)?;
        for _ in 0..np {
            let p = fl5_polar(&mut r)?;
            push_polar(&mut imp, p);
        }
    } else {
        let nf = r.count("number of foils", 1)?;
        for _ in 0..nf {
            let f = plr_foil(&mut r)?;
            push_foil(&mut imp, f);
        }
        let np = r.count("number of polars", 1)?;
        for _ in 0..np {
            let p = plr_polar(&mut r)?;
            push_polar(&mut imp, p);
        }
    }
    Ok(imp)
}

/// `xfl::serializeFoil` (xflr5) / `serial::serializeFoil` (flow5).
fn plr_foil(r: &mut Rd) -> Res<RawFoil> {
    let w = "plr foil";
    let fmt = r.fmt(w, 1000, 1010)?;
    let name = r.cstr(w)?;
    let desc = if fmt >= 1006 { r.cstr(w)? } else { String::new() };
    if fmt >= 1002 {
        r.skip_i32s(3, w)?; // stipple, width, color
    }
    if fmt >= 1003 {
        r.i32(w)?;
    }
    if fmt >= 1004 {
        r.skip_i32s(2, w)?;
    }
    if fmt >= 1005 {
        r.i32(w)?;
        r.skip_flts(3, w)?;
    }
    let te_flap = r.i32(w)? != 0;
    let te_angle = r.flt(w)?;
    let xh = r.flt(w)?;
    let yh = r.flt(w)?;
    r.skip_flts(3, w)?;
    let n = r.count("plr foil point count", 2 * r.flt_size())?;
    let (mut x, mut y) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        x.push(r.flt(w)?);
        y.push(r.flt(w)?);
    }
    if fmt >= 1001 {
        let n = r.loose_count("plr foil current point count", 2 * r.flt_size())?;
        r.skip_flts(2 * n, w)?;
    }
    // format 1007 (flow5) stores the hinge as fractions, older ones in %
    let k = if fmt >= 1007 { 1.0 } else { 0.01 };
    Ok(RawFoil { name, desc, te_flap, te_angle, hinge: (xh * k, yh * k), x, y })
}

/// `xfl::serializePolar` (xflr5) / `serial::serializePolarv6` (flow5).
fn plr_polar(r: &mut Rd) -> Res<RawPolar> {
    let w = "plr polar";
    let fmt = r.fmt(w, 1001, 1100)?;
    let at = r.p;
    let foil = r.cstr(w)?;
    let name = r.cstr(w)?;
    if foil.is_empty() || name.is_empty() {
        return r.fail(at, w, "empty foil or polar name");
    }
    let ptype = match r.i32(w)? {
        2 => 2,
        3 => 3,
        4 => 4,
        _ => 1,
    };
    r.skip_i32s(2, w)?;
    let re = r.i32(w)? as f64;
    let mach = r.flt(w)?;
    r.flt(w)?; // alpha spec
    let n = r.loose_count("plr polar point count", 10 * r.flt_size())?;
    let mut ncrit = r.flt(w)?;
    let mut xtop = r.flt(w)?;
    let mut xbot = r.flt(w)?;
    if fmt < 1005 {
        r.skip_i32s(3, w)?; // color, stipple, width
        if fmt >= 1002 {
            let at = r.p;
            let l = r.i32(w)?;
            if l != 0 && l != 1 {
                return r.fail(at, w, format_args!("invalid visibility flag {l}"));
            }
        }
        r.i32(w)?;
    } else {
        style_xfl(r)?;
    }
    let mut pts = Vec::with_capacity(n);
    for _ in 0..n {
        let alpha = r.flt(w)?;
        let cd = r.flt(w)?;
        let cdp = r.flt(w)?;
        let cl = r.flt(w)?;
        let cm = r.flt(w)?;
        let xtr1 = r.flt(w)?;
        let xtr2 = r.flt(w)?;
        r.flt(w)?; // HMom
        let cpmn = r.flt(w)?;
        let pre = r.flt(w)?;
        if fmt >= 1004 {
            r.flt(w)?; // XCp
        }
        // xflr5 ignores duplicate points (same alpha, or same Re for type 4)
        let dup = pts.iter().any(|q: &RawPt| if ptype == 4 { (q.re - pre).abs() < 0.1 } else { (q.alpha - alpha).abs() < 0.001 });
        if !dup {
            pts.push(RawPt { alpha, cd, cdp, cl, cm, xtr1, xtr2, cpmn, re: pre, ctrl: 0.0 });
        }
    }
    if fmt >= 1003 {
        ncrit = r.f64(w)?;
        xtop = r.f64(w)?;
        xbot = r.f64(w)?;
    }
    Ok(RawPolar { foil, name, ptype, re, mach, ncrit, xtop, xbot, viscous: true, te_flap: 0.0, pts })
}

// ---------------------------------------------------------------------------
// flow5 .fl5
// ---------------------------------------------------------------------------

/// Read a flow5 project file (`.fl5`).
pub fn read_fl5(bytes: &[u8]) -> Result<LegacyImport, String> {
    let mut r = Rd::new(bytes, false, true);
    let mut imp = LegacyImport::default();

    // --- project meta data ---
    let fmt = r.fmt("fl5 project header", 500001, 501000)?;
    fl5_polar(&mut r)?;
    if fmt <= 500006 {
        fl5_planepolar_v726(&mut r)?;
    } else {
        fl5_planepolar_v750(&mut r)?;
    }
    if fmt >= 500002 {
        r.qstr("project meta data")?;
    }
    if fmt >= 500003 {
        fl5_splinefoil(&mut r)?;
        fl5_bspline(&mut r)?;
        fl5_bspline(&mut r)?;
    }
    if fmt >= 500006 {
        fl5_bspline(&mut r)?;
        fl5_cubicspline(&mut r)?;
    }
    let n = r.loose_count("project meta data spares", 4)?;
    r.skip_i32s(n, "project meta data spares")?;
    if fmt >= 500002 {
        let n = r.loose_count("project meta data spares", 8)?;
        r.skip_f64s(n, "project meta data spares")?;
    }

    // --- 2D objects ---
    if fmt >= 500750 {
        r.i32("2d objects format")?;
    }
    let n = r.count("number of foils", 1)?;
    for _ in 0..n {
        let f = fl5_foil(&mut r)?;
        push_foil(&mut imp, f);
    }
    let n = r.count("number of foil polars", 1)?;
    for _ in 0..n {
        let p = fl5_polar(&mut r)?;
        push_polar(&mut imp, p);
    }
    let n = r.count("number of foil operating points", 1)?;
    for _ in 0..n {
        r.i32("foil operating point")?;
        fl5_opp(&mut r)?;
    }

    // --- planes ---
    let n = r.count("number of planes", 1)?;
    for _ in 0..n {
        let kind = if fmt >= 500004 { r.i32("plane type")? } else { 0 };
        match kind {
            0 => {
                let mut warnings = Vec::new();
                let plane = fl5_plane_xfl(&mut r, &mut warnings)?;
                imp.warnings.extend(warnings);
                if let Some(p) = plane {
                    imp.planes.push(p);
                }
            }
            1 => {
                let name = fl5_plane_stl(&mut r)?;
                imp.warnings.push(format!("plane '{name}' is an STL (triangulated) plane and was skipped"));
            }
            k => return r.fail(r.p - 4, "plane type", format_args!("unknown plane type {k}")),
        }
    }
    // plane polars, plane operating points and boats follow: not needed.
    Ok(imp)
}

/// `serial::serializeFoilFl5`.
fn fl5_foil(r: &mut Rd) -> Res<RawFoil> {
    let w = "fl5 foil";
    let fmt = r.fmt(w, 500000, 550000)?;
    let name = r.qstr(w)?;
    let desc = r.qstr(w)?;
    r.skip_i32s(3, w)?;
    r.color(w)?;
    r.bool(w)?; // visible
    r.skip_bools(2, w)?; // camber line, LE flap
    let te_flap = r.bool(w)?;
    r.skip_f64s(3, w)?;
    let te_angle = r.f64(w)?;
    let xh = r.f64(w)?;
    let yh = r.f64(w)?;
    if fmt >= 500753 {
        r.f64(w)?;
        r.i32(w)?;
    }
    let n = r.count("fl5 foil point count", 16)?;
    let (mut x, mut y) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        x.push(r.f64(w)?);
        y.push(r.f64(w)?);
    }
    // flow5 >= 7.50 stores flaps made permanent (deflection reset to 0)
    let te_angle = if fmt >= 500750 { 0.0 } else { te_angle };
    Ok(RawFoil { name, desc, te_flap, te_angle, hinge: (xh, yh), x, y })
}

/// `serial::serializePolarFl5`.
fn fl5_polar(r: &mut Rd) -> Res<RawPolar> {
    let w = "fl5 foil polar";
    let fmt = r.fmt(w, 500000, 501000)?;
    let foil = r.qstr(w)?;
    let name = r.qstr(w)?;
    style_fl5(r)?;
    let ptype = match r.i32(w)? {
        2 => 2,
        3 => 3,
        4 => 4,
        5 | 6 => 6,
        _ => 1,
    };
    let viscous = r.i32(w)? == 0;
    r.skip_i32s(2, w)?; // MaType, ReType
    let re = if fmt < 500750 {
        let _rho = r.f64(w)?;
        let nu = r.f64(w)?;
        let qinf = r.f64(w)?;
        let chord = r.f64(w)?;
        if nu.abs() > 0.0 {
            qinf * chord / nu
        } else {
            0.0
        }
    } else {
        r.f64(w)?
    };
    let mach = r.f64(w)?;
    r.f64(w)?; // aoa spec
    let xtop = r.f64(w)?;
    let xbot = r.f64(w)?;
    let ncrit = r.f64(w)?;
    let nc = r.loose_count("fl5 polar control count", 16)?;
    r.skip_f64s(2 * nc, w)?;
    r.i32(w)?; // nVariables
    let n = r.loose_count("fl5 polar point count", 16 * r.flt_size())?;
    let mut pts = Vec::with_capacity(n);
    for _ in 0..n {
        let mut v = [0.0; 16];
        for x in v.iter_mut() {
            *x = r.flt(w)?;
        }
        // Alpha Cd Cdp Cl Cm HMom Cpmn Re XCp Ctrl XTr1 XTr2 ...
        pts.push(RawPt { alpha: v[0], cd: v[1], cdp: v[2], cl: v[3], cm: v[4], cpmn: v[6], re: v[7], ctrl: v[9], xtr1: v[10], xtr2: v[11] });
    }
    let n = r.loose_count("fl5 polar spares", 4)?;
    r.skip_i32s(n, w)?;
    let nd = r.i32(w)?;
    let mut te_flap = if nd > 0 { r.f64(w)? } else { 0.0 };
    if ptype == 6 || !te_flap.is_finite() {
        te_flap = 0.0;
    }
    Ok(RawPolar { foil, name, ptype, re, mach, ncrit, xtop, xbot, viscous, te_flap, pts })
}

/// `serial::serializeOppFl5` (parsed and discarded).
fn fl5_opp(r: &mut Rd) -> Res<()> {
    let w = "fl5 foil operating point";
    let fmt = r.fmt(w, 500000, 550000)?;
    r.qstr(w)?;
    r.qstr(w)?;
    style_fl5(r)?;
    r.i32(w)?; // BL method
    r.skip_f64s(3, w)?;
    if fmt >= 500750 {
        r.f64(w)?;
    }
    r.skip_bools(2, w)?;
    r.skip_f64s(10, w)?;
    let fs = r.flt_size();
    if fmt >= 500003 {
        let n = r.loose_count("fl5 operating point node count", 5 + 4 * fs)?;
        for _ in 0..n {
            r.i32(w)?;
            r.bool(w)?;
            r.skip_flts(4, w)?;
        }
    }
    let n = r.loose_count("fl5 operating point surface point count", 2 * fs)?;
    if fmt >= 500002 {
        r.skip_flts(2 * n, w)?;
    }
    r.skip_flts(2 * n, w)?;
    if fmt >= 500004 {
        fl5_blxfoil(r)?;
    } else {
        fl5_bldata(r)?;
        fl5_bldata(r)?;
    }
    r.spares(w)
}

/// `serial::serializeBLXFoil`.
fn fl5_blxfoil(r: &mut Rd) -> Res<()> {
    let w = "fl5 XFoil boundary layer";
    r.i32(w)?;
    r.skip_i32s(2, w)?; // nside1, nside2
    let nd1 = r.i32(w)?;
    let nd2 = r.i32(w)?;
    let nd3 = r.i32(w)?;
    let fs = r.flt_size();
    for n in [nd1.saturating_add(1), nd2, nd3] {
        let n = n.max(0) as usize;
        if n.saturating_mul(2 * fs) > r.remaining() {
            return r.fail(r.p, w, format_args!("point count {n} exceeds the data size"));
        }
        r.skip_flts(2 * n, w)?;
    }
    r.spares(w)
}

/// `serial::serializeBLDataFl5`.
fn fl5_bldata(r: &mut Rd) -> Res<()> {
    let w = "fl5 boundary layer data";
    r.skip_i32s(5, w)?; // format, method, side, iLE, nTE
    r.bool(w)?;
    r.skip_f64s(8, w)?;
    let n = r.loose_count("fl5 boundary layer point count", 0)?;
    let nv = r.loose_count("fl5 boundary layer variable count", 0)?;
    let total = n.saturating_mul(nv);
    if total.saturating_mul(r.flt_size()) > r.remaining() {
        return r.fail(r.p, w, format_args!("{n} x {nv} values exceed the data size"));
    }
    r.skip_flts(total, w)?;
    r.spares(w)
}

/// `serial::serializeSplineFl5`.
fn fl5_spline(r: &mut Rd) -> Res<()> {
    let w = "fl5 spline";
    r.fmt(w, 500001, 500001)?;
    style_fl5(r)?;
    r.bool(w)?;
    r.i32(w)?; // output size
    let n = r.loose_count("fl5 spline control point count", 16)?;
    r.skip_f64s(2 * n, w)?;
    let n = r.loose_count("fl5 spline weight count", 8)?;
    r.skip_f64s(n, w)?;
    r.skip_bools(2, w)?;
    r.f64(w)?;
    r.skip_i32s(2, w)?;
    r.spares(w)
}

fn fl5_bspline(r: &mut Rd) -> Res<()> {
    fl5_spline(r)?;
    r.skip_i32s(3, "fl5 B-spline")
}

fn fl5_cubicspline(r: &mut Rd) -> Res<()> {
    fl5_spline(r)?;
    r.spares("fl5 cubic spline")
}

/// `serial::serializeSplineFoilFl5`.
fn fl5_splinefoil(r: &mut Rd) -> Res<()> {
    let w = "fl5 spline foil";
    r.fmt(w, 500000, 510000)?;
    style_fl5(r)?;
    fl5_bspline(r)?;
    fl5_bspline(r)?;
    r.skip_bools(5, w)?;
    r.spares(w)
}

/// `serial::serializeAngleControlFl5`.
fn fl5_angle_control(r: &mut Rd) -> Res<()> {
    let w = "fl5 angle control";
    r.i32(w)?;
    r.qstr(w)?;
    let n = r.loose_count("fl5 angle control value count", 8)?;
    r.skip_f64s(n, w)
}

/// A list of (name, min, max) control ranges.
fn fl5_ranges(r: &mut Rd, w: &str) -> Res<()> {
    let n = r.loose_count(w, 20)?;
    for _ in 0..n {
        r.qstr(w)?;
        r.skip_f64s(2, w)?;
    }
    Ok(())
}

/// `serial::serializeAeroForcesFl5`.
fn fl5_aeroforces(r: &mut Rd) -> Res<()> {
    let w = "fl5 aero forces";
    let fmt = r.fmt(w, 500000, 510000)?;
    r.skip_f64s(24, w)?;
    if fmt >= 500750 {
        r.skip_f64s(3, w)?;
    }
    Ok(())
}

/// `serial::serializePolar3dFl5v726`; returns the polar format.
fn fl5_polar3d_v726(r: &mut Rd) -> Res<i32> {
    let w = "fl5 3D polar (v7.26)";
    let fmt = r.fmt(w, 500000, 500100)?;
    r.qstr(w)?;
    style_fl5(r)?;
    r.i32(w)?; // method
    r.skip_bools(4, w)?;
    r.i32(w)?; // type
    let n = r.loose_count("fl5 extra drag count", 20)?;
    for _ in 0..n {
        r.qstr(w)?;
        r.skip_f64s(2, w)?;
    }
    r.bool(w)?;
    if fmt >= 500029 {
        r.bool(w)?;
    }
    r.skip_f64s(6, w)?; // height, ?, density, viscosity, beta, phi
    r.bool(w)?;
    r.skip_f64s(1 + 3 + 6, w)?; // mass, CoG, inertia
    r.skip_bools(2, w)?;
    r.skip_f64s(3, w)?;
    r.skip_i32s(2, w)?;
    r.skip_f64s(2, w)?;
    r.i32(w)?;
    let n_int = r.i32(w)?;
    if n_int >= 1 {
        r.i32(w)?;
    }
    if n_int >= 2 {
        r.i32(w)?;
    }
    // C++: for(int n=2; n<nIntExtra; n++) ar >> n;  (reads into the loop counter)
    let mut k = 2i32;
    while k < n_int {
        k = r.i32(w)?;
        k = k.saturating_add(1);
    }
    let n_dbl = r.i32(w)?;
    if n_dbl >= 2 {
        r.skip_f64s(2, w)?;
    }
    if n_dbl >= 4 {
        r.skip_f64s(2, w)?;
    }
    if fmt >= 500015 {
        r.f64(w)?;
    }
    if fmt >= 500016 {
        r.f64(w)?;
    }
    if fmt >= 500019 {
        r.i32(w)?;
    }
    Ok(fmt)
}

/// `serial::serializePolar3dFl5v750`.
fn fl5_polar3d_v750(r: &mut Rd) -> Res<()> {
    let w = "fl5 3D polar";
    r.fmt(w, 500750, 501000)?;
    r.qstr(w)?;
    style_fl5(r)?;
    r.i32(w)?;
    r.skip_bools(3, w)?;
    r.i32(w)?;
    let n = r.loose_count("fl5 extra drag count", 20)?;
    for _ in 0..n {
        r.qstr(w)?;
        r.skip_f64s(2, w)?;
    }
    r.skip_bools(2, w)?;
    r.skip_f64s(6, w)?;
    r.bool(w)?;
    r.skip_f64s(10, w)?;
    r.skip_bools(2, w)?;
    r.skip_f64s(3, w)?;
    r.bool(w)?;
    r.i32(w)?;
    r.skip_f64s(2, w)?;
    r.bool(w)?;
    r.skip_f64s(4, w)?;
    r.i32(w)?;
    r.skip_bools(10, w)?;
    r.skip_i32s(20, w)?;
    r.skip_f64s(20, w)?;
    Ok(())
}

/// `serial::serializePlanePolarFl5v726` (parsed and discarded).
fn fl5_planepolar_v726(r: &mut Rd) -> Res<()> {
    let w = "fl5 plane polar (v7.26)";
    let at = r.p;
    let fmt = fl5_polar3d_v726(r)?;
    if fmt < 500001 {
        return r.fail(at, w, format_args!("unsupported format number {fmt}"));
    }
    r.qstr(w)?;
    r.skip_f64s(2, w)?;
    r.bool(w)?;
    r.i32(w)?;
    r.skip_f64s(3, w)?;
    if fmt >= 500021 {
        let n = r.loose_count("fl5 AVL control count", 12)?;
        for _ in 0..n {
            fl5_angle_control(r)?;
        }
    }
    let nc = r.count("fl5 deprecated control count", 4)?;
    for _ in 0..nc {
        let k = r.loose_count("fl5 deprecated control gain count", 8)?;
        let per = if fmt >= 500007 { 2 } else { 1 };
        r.skip_f64s(k * per, w)?;
    }
    r.skip_f64s(7, w)?;
    r.i32(w)?;
    let n = r.loose_count("fl5 angle range count", 4)?;
    for _ in 0..n {
        fl5_ranges(r, w)?;
    }
    if fmt >= 500002 {
        fl5_ranges(r, w)?;
    }
    if fmt >= 500003 {
        let n = r.loose_count("fl5 range group count", 4)?;
        for _ in 0..n {
            fl5_ranges(r, w)?;
        }
    }
    fl5_ranges(r, w)?; // operating range
    fl5_ranges(r, w)?; // inertia range
    if fmt >= 500004 {
        r.i32(w)?;
    }
    r.skip_i32s(2, w)?;
    let n_spares = r.i32(w)?.max(0) as usize;
    let at = r.p;
    let n = r.i32(w)?;
    if n.abs() > 10000 {
        return r.fail(at, w, format_args!("{n} data points"));
    }
    for _ in 0..n.max(0) {
        if fmt < 500022 {
            r.skip_f64s(18, w)?;
        } else {
            fl5_aeroforces(r)?;
        }
        r.skip_f64s(6 + 16, w)?;
        r.skip_f64s(n_spares, w)?;
    }
    if fmt >= 500017 {
        r.bool(w)?;
        fl5_bspline(r)?;
    }
    let ni = r.i32(w)?;
    if ni >= 1 {
        r.i32(w)?;
    }
    if ni >= 2 {
        r.i32(w)?;
    }
    if r.i32(w)? >= 1 {
        r.f64(w)?;
    }
    Ok(())
}

/// `serial::serializePlanePolarFl5v750` (parsed and discarded).
fn fl5_planepolar_v750(r: &mut Rd) -> Res<()> {
    let w = "fl5 plane polar";
    fl5_polar3d_v750(r)?;
    r.qstr(w)?;
    r.skip_f64s(2, w)?;
    r.skip_bools(2, w)?;
    r.i32(w)?;
    r.skip_f64s(3, w)?;
    for _ in 0..2 {
        let n = r.loose_count("fl5 control count", 12)?;
        for _ in 0..n {
            fl5_angle_control(r)?;
        }
    }
    r.bool(w)?;
    let n = r.loose_count("fl5 angle range count", 4)?;
    for _ in 0..n {
        fl5_ranges(r, w)?;
    }
    fl5_ranges(r, w)?;
    fl5_ranges(r, w)?;
    r.skip_bools(2, w)?;
    r.i32(w)?;
    let n_spares = r.i32(w)?.max(0) as usize;
    let at = r.p;
    let n = r.i32(w)?;
    if n.abs() > 10000 {
        return r.fail(at, w, format_args!("{n} data points"));
    }
    for _ in 0..n.max(0) {
        fl5_aeroforces(r)?;
        r.skip_f64s(7 + 16, w)?;
        r.skip_f64s(n_spares, w)?;
    }
    r.bool(w)?;
    fl5_bspline(r)?;
    r.skip_bools(2, w)?;
    r.f64(w)?;
    r.skip_bools(10, w)?;
    r.skip_i32s(20, w)?;
    r.skip_f64s(20, w)?;
    Ok(())
}

struct RawInertia {
    mass: f64,
    pms: Vec<PointMass>,
}

/// `serial::serializeInertiaFl5`.
fn fl5_inertia(r: &mut Rd) -> Res<RawInertia> {
    let w = "fl5 inertia";
    r.fmt(w, 500001, 500001)?;
    let mass = r.f64(w)?;
    r.skip_f64s(9, w)?;
    let n = r.count("fl5 point mass count", 36)?;
    let mut pms = Vec::with_capacity(n);
    for _ in 0..n {
        let m = r.f64(w)?;
        let p = r.v3(w)?;
        let tag = r.qstr(w)?;
        pms.push(pm(tag, m, p));
    }
    r.spares(w)?;
    Ok(RawInertia { mass, pms })
}

struct RawPart {
    name: String,
    inertia: RawInertia,
}

/// `serial::serializePartFl5`.
fn fl5_part(r: &mut Rd) -> Res<RawPart> {
    let w = "fl5 part";
    let fmt = r.fmt(w, 500000, 501000)?;
    let name = r.qstr(w)?;
    r.qstr(w)?;
    style_fl5(r)?;
    let inertia = fl5_inertia(r)?;
    r.bool(w)?;
    if fmt >= 500002 {
        r.bool(w)?;
    }
    if fmt >= 500754 {
        for _ in 0..2 {
            r.skip_f64s(2, w)?;
            r.i32(w)?;
        }
    }
    let at = r.p;
    let n = r.i32(w)?;
    if !(0..=100).contains(&n) {
        return r.fail(at, w, format_args!("invalid spare count {n}"));
    }
    r.skip_i32s(n as usize, w)?;
    let at = r.p;
    let n = r.i32(w)?;
    if !(0..=100).contains(&n) {
        return r.fail(at, w, format_args!("invalid spare count {n}"));
    }
    r.skip_f64s(n as usize, w)?;
    Ok(RawPart { name, inertia })
}

/// `serial::serializeWingSectionFl5`.
fn fl5_wing_section(r: &mut Rd) -> Res<WingSection> {
    let w = "fl5 wing section";
    r.fmt(w, 500001, 500001)?;
    let right = r.qstr(w)?;
    let left = r.qstr(w)?;
    let chord = r.f64(w)?;
    let y = r.f64(w)?;
    let offset = r.f64(w)?;
    let dihedral = r.f64(w)?;
    let twist = r.f64(w)?;
    let nx = r.i32(w)?;
    let ny = r.i32(w)?;
    let xd = r.i32(w)?;
    let yd = r.i32(w)?;
    let n = r.loose_count("fl5 wing section spares", 12)?;
    r.skip_i32s(n, w)?;
    r.skip_f64s(n, w)?;
    let mut s = WingSection::new(y, chord, offset, &right);
    s.foil_left = left;
    s.dihedral = dihedral;
    s.twist = twist;
    s.nx = panels(nx, 1);
    s.ny = panels(ny, 1);
    s.xdist = flow5_dist(xd);
    s.ydist = flow5_dist(yd);
    Ok(s)
}

/// `serial::serializeWingXflFl5`; returns the wing (None if unusable).
fn fl5_wing(r: &mut Rd, plane: &str, warnings: &mut Vec<String>) -> Res<Option<Wing>> {
    let w = "fl5 wing";
    let part = fl5_part(r)?;
    r.fmt(w, 500001, 500002)?;
    let symmetric = r.bool(w)?;
    let n = r.count("fl5 wing section count", 60)?;
    let mut sections = Vec::with_capacity(n);
    for _ in 0..n {
        sections.push(fl5_wing_section(r)?);
    }
    let two_sided = r.bool(w)?;
    let close_inner = r.bool(w)?;
    r.bool(w)?; // reversed (deprecated)
    let kind = match r.i32(w)? {
        0 => WingKind::Main,
        2 => WingKind::Elevator,
        3 => WingKind::Fin,
        _ => WingKind::Other,
    };
    let pos = r.v3(w)?;
    let rot = r.v3(w)?;
    let n = r.loose_count("fl5 wing spares", 4)?;
    r.skip_i32s(n, w)?;
    let n = r.loose_count("fl5 wing spares", 8)?;
    r.skip_f64s(n, w)?;

    let name = if part.name.trim().is_empty() { kind.label().to_string() } else { part.name.trim().to_string() };
    if sections.len() < 2 {
        warnings.push(format!("plane '{plane}': wing '{name}' skipped (fewer than 2 sections)"));
        return Ok(None);
    }
    if rot.z.abs() > 1e-9 {
        warnings.push(format!("plane '{plane}': wing '{name}': rotation about z ({:.2} deg) ignored", rot.z));
    }
    let mut wing = Wing::new(&name, kind, sections);
    wing.position = pos;
    wing.symmetric = symmetric;
    wing.two_sided = two_sided;
    wing.close_inner = close_inner;
    // flow5 builds one-sided wings on the left (-y) half, we build them on the right
    wing.rx = if two_sided { rot.x } else { -rot.x };
    wing.ry = rot.y;
    wing.mass = part.inertia.mass;
    wing.point_masses = part.inertia.pms; // relative to the wing position, as in flow5
    Ok(Some(wing))
}

struct Fl5FuseBase {
    part: RawPart,
    pos: Vec3,
    rot: Vec3,
}

/// `serial::serializeFuseFl5` (incl. the part).
fn fl5_fuse_base(r: &mut Rd) -> Res<Fl5FuseBase> {
    let w = "fl5 fuselage";
    let part = fl5_part(r)?;
    let fmt = r.fmt(w, 500001, 500003)?;
    let pos = r.v3(w)?;
    let rot = r.v3(w)?;
    if fmt >= 500002 {
        // serializeOccParams
        r.i32(w)?;
        r.bool(w)?;
        r.skip_f64s(6, w)?;
        r.skip_i32s(4, w)?;
    }
    if fmt >= 500003 {
        r.f64(w)?;
    }
    r.skip_i32s(2, w)?;
    Ok(Fl5FuseBase { part, pos, rot })
}

enum Fl5Fuse {
    Xfl { base: Fl5FuseBase, flat: bool, frames: Vec<NurbsFrame>, udeg: i32, vdeg: i32, nx: i32, nh: i32, hpanels: Vec<i32>, xpanels: Vec<i32> },
    Tris { base: Fl5FuseBase, tris: Vec<[Vec3; 3]>, kind: &'static str },
}

/// `serial::serializeFuseXflFl5` with `serializeNURBSFl5` / `serializeFrameFl5`.
fn fl5_fuse_xfl(r: &mut Rd) -> Res<Fl5Fuse> {
    let w = "fl5 xfl fuselage";
    let base = fl5_fuse_base(r)?;
    r.fmt(w, 500001, 500001)?;
    let ftype = r.i32(w)?;
    // NURBS surface
    r.fmt("fl5 NURBS surface", 500001, 500001)?;
    r.color(w)?;
    let udeg = r.i32(w)?;
    let vdeg = r.i32(w)?;
    r.skip_f64s(4, w)?;
    let uaxis = r.i32(w)?;
    r.i32(w)?; // vaxis
    r.i32(w)?; // output
    let nf = r.count("fl5 NURBS frame count", 48)?;
    let mut frames = Vec::with_capacity(nf);
    let ui = if (0..=2).contains(&uaxis) { uaxis } else { 0 };
    let comp = |v: Vec3, i: i32| match i {
        1 => v.y,
        2 => v.z,
        _ => v.x,
    };
    for _ in 0..nf {
        r.fmt("fl5 frame", 500001, 500001)?;
        let fpos = r.v3(w)?;
        r.f64(w)?; // angle
        let np = r.count("fl5 frame point count", 24)?;
        let mut points = Vec::with_capacity(np);
        for _ in 0..np {
            points.push(r.v3(w)?);
        }
        r.spares(w)?;
        let x = comp(fpos, ui);
        for p in &mut points {
            // Frame::setuPosition: the frame's u coordinate is its position
            match ui {
                1 => p.y = x,
                2 => p.z = x,
                _ => p.x = x,
            }
        }
        frames.push(NurbsFrame { x, points });
    }
    r.spares(w)?;
    let nx = r.i32(w)?;
    let nh = r.i32(w)?;
    let side_lines = frames.first().map(|f| f.points.len()).unwrap_or(0);
    let mut hpanels = Vec::with_capacity(side_lines);
    for _ in 0..side_lines {
        hpanels.push(r.i32(w)?);
    }
    let mut xpanels = Vec::with_capacity(frames.len());
    for _ in 0..frames.len() {
        xpanels.push(r.i32(w)?);
    }
    r.spares(w)?;
    if ui != 0 {
        // a body built along another axis: keep x as the frame station
        for f in &mut frames {
            for p in &mut f.points {
                *p = Vec3::new(f.x, p.y, p.z);
            }
        }
    }
    Ok(Fl5Fuse::Xfl { base, flat: ftype == 1, frames, udeg, vdeg, nx, nh, hpanels, xpanels })
}

/// `serial::serializeFuseStlFl5`.
fn fl5_fuse_stl(r: &mut Rd) -> Res<Fl5Fuse> {
    let w = "fl5 STL fuselage";
    let base = fl5_fuse_base(r)?;
    let fmt = r.fmt(w, 500001, 500010)?;
    let dbl = fmt < 500002;
    let vs = if dbl { 8 } else { r.flt_size() };
    let n = r.count("fl5 STL triangle count", 9 * vs)?;
    let mut tris = Vec::with_capacity(n);
    for _ in 0..n {
        let mut t = [Vec3::ZERO; 3];
        for v in t.iter_mut() {
            *v = if dbl { r.v3(w)? } else { r.v3f(w)? };
        }
        tris.push(t);
    }
    r.spares(w)?;
    Ok(Fl5Fuse::Tris { base, tris, kind: "STL" })
}

/// `serial::serializeFuseOccFl5` (the BREP shapes are skipped; the panel
/// mesh stored after the fuse is used instead).
fn fl5_fuse_occ(r: &mut Rd) -> Res<Fl5Fuse> {
    let w = "fl5 CAD fuselage";
    let base = fl5_fuse_base(r)?;
    r.fmt(w, 500001, 500100)?;
    let n = r.count("fl5 CAD shape count", 4)?;
    for _ in 0..n {
        r.qstr(w)?;
    }
    r.spares(w)?;
    Ok(Fl5Fuse::Tris { base, tris: vec![], kind: "CAD" })
}

/// `serial::serializePanelsFl5`: triangles in fuse axes.
fn fl5_panels(r: &mut Rd) -> Res<Vec<[Vec3; 3]>> {
    let w = "fl5 fuselage panels";
    r.i32(w)?;
    let n = r.count("fl5 panel count", 9 * r.flt_size() + 1)?;
    let mut tris = Vec::with_capacity(n);
    for _ in 0..n {
        let t = [r.v3f(w)?, r.v3f(w)?, r.v3f(w)?];
        r.bool(w)?;
        tris.push(t);
    }
    r.spares(w)?;
    Ok(tris)
}

/// `serial::serializePlaneXflFl5`.
fn fl5_plane_xfl(r: &mut Rd, warnings: &mut Vec<String>) -> Res<Option<Plane>> {
    let w = "fl5 plane";
    let fmt = r.fmt(w, 500000, 500010)?;
    let nw = r.count("fl5 plane wing count", 100)?;
    let name = r.qstr(w)?.trim().to_string();
    let desc = r.qstr(w)?;
    r.skip_i32s(3, w)?;
    r.color(w)?;
    let pname = if name.is_empty() { "Imported plane".to_string() } else { name.clone() };
    let mut plane = empty_plane(name, desc);
    for _ in 0..nw {
        if let Some(wing) = fl5_wing(r, &pname, warnings)? {
            plane.wings.push(wing);
        }
    }
    let at = r.p;
    let nf = r.i32(w)?;
    if !(0..=10000).contains(&nf) {
        return r.fail(at, w, format_args!("invalid fuselage count {nf}"));
    }
    for _ in 0..nf {
        let at = r.p;
        let f = match r.i32("fl5 fuselage type")? {
            100001 | 100004 | 100005 | 100006 => fl5_fuse_xfl(r)?,
            100002 => fl5_fuse_occ(r)?,
            100003 => fl5_fuse_stl(r)?,
            k => return r.fail(at, "fl5 fuselage type", format_args!("unknown fuselage type {k}")),
        };
        let mesh = fl5_panels(r)?;
        match f {
            Fl5Fuse::Xfl { base, flat, frames, udeg, vdeg, nx, nh, hpanels, xpanels } => {
                let fname = base.part.name.clone();
                if base.rot.norm() > 1e-9 {
                    warnings.push(format!("plane '{pname}': fuselage '{fname}': rotations ignored"));
                }
                for m in &base.part.inertia.pms {
                    plane.point_masses.push(pm(m.tag.clone(), m.mass, m.position + base.pos));
                }
                let (nx, nh) = if flat { (xpanels.iter().map(|&k| k.max(0) as usize).sum(), hpanels.iter().map(|&k| k.max(0) as usize).sum()) } else { (panels(nx, 0), panels(nh, 0)) };
                let spec = FuseSpec { name: fname, position: base.pos, frames, flat, degree_u: udeg, degree_v: vdeg, nx, nh, mass: base.part.inertia.mass };
                if let Some(f) = make_fuse(spec, warnings) {
                    plane.fuses.push(f);
                }
            }
            Fl5Fuse::Tris { base, tris, kind } => {
                let fname = if base.part.name.trim().is_empty() { "Fuselage".to_string() } else { base.part.name.trim().to_string() };
                if base.rot.norm() > 1e-9 {
                    warnings.push(format!("plane '{pname}': fuselage '{fname}': rotations ignored"));
                }
                for m in &base.part.inertia.pms {
                    plane.point_masses.push(pm(m.tag.clone(), m.mass, m.position + base.pos));
                }
                let src = if tris.len() >= 4 { tris } else { mesh };
                let world: Vec<[Vec3; 3]> = src.into_iter().filter(|t| t.iter().all(|v| v.is_finite())).map(|t| [t[0] + base.pos, t[1] + base.pos, t[2] + base.pos]).collect();
                if world.len() < 4 {
                    warnings.push(format!("plane '{pname}': {kind} fuselage '{fname}' skipped (no usable triangulation)"));
                } else {
                    if kind == "CAD" {
                        warnings.push(format!("plane '{pname}': CAD fuselage '{fname}' imported as its stored triangle mesh"));
                    }
                    plane.fuses.push(Fuse::from_triangles(&fname, world, base.part.inertia.mass));
                }
            }
        }
    }
    if fmt >= 500002 {
        r.bool(w)?;
        let i = fl5_inertia(r)?;
        plane.point_masses.extend(i.pms);
    } else {
        let n = r.count("fl5 plane point mass count", 36)?;
        for _ in 0..n {
            let m = r.f64(w)?;
            let p = r.v3(w)?;
            let tag = r.qstr(w)?;
            plane.point_masses.push(pm(tag, m, p));
        }
    }
    let n = r.loose_count("fl5 plane spares", 4)?;
    r.skip_i32s(n, w)?;
    let n = r.loose_count("fl5 plane spares", 8)?;
    r.skip_f64s(n, w)?;
    Ok(finish_plane(plane, warnings))
}

/// `serial::serializeMeshFl5` (parsed and discarded).
fn fl5_mesh(r: &mut Rd) -> Res<()> {
    let w = "fl5 mesh";
    let fmt = r.i32(w)?;
    let fs = r.flt_size();
    let n = r.loose_count("fl5 mesh node count", 4 + 6 * fs + 9)?;
    for _ in 0..n {
        r.i32(w)?;
        r.skip_flts(6, w)?;
        r.bool(w)?;
        for _ in 0..2 {
            let k = r.loose_count("fl5 mesh node neighbour count", 4)?;
            r.skip_i32s(k, w)?;
        }
    }
    let per = 4 * 7 + 2 + if fmt >= 500002 { 4 } else { 0 } + 4;
    let n = r.loose_count("fl5 mesh panel count", per)?;
    r.skip(n * per, w)
}

/// `serial::serializePlaneStlFl5` (parsed and discarded); returns the name.
fn fl5_plane_stl(r: &mut Rd) -> Res<String> {
    let w = "fl5 STL plane";
    let fmt = r.fmt(w, 500000, 500100)?;
    let name = r.qstr(w)?.trim().to_string();
    r.qstr(w)?;
    if fmt >= 500004 {
        r.color(w)?;
    }
    if fmt < 500002 {
        r.skip_i32s(3, w)?;
        r.color(w)?;
    } else {
        style_fl5(r)?;
        r.bool(w)?;
        r.skip_f64s(7, w)?;
        fl5_inertia(r)?;
    }
    if fmt <= 500002 {
        let n = r.loose_count("fl5 STL plane triangle count", 9 * r.flt_size())?;
        r.skip_flts(9 * n, w)?;
        r.spares(w)?;
    } else {
        fl5_mesh(r)?;
        if fmt >= 500005 {
            r.fmt("fl5 triangulation", 500001, 500010)?;
            let n = r.loose_count("fl5 triangulation triangle count", 9 * r.flt_size())?;
            r.skip_flts(9 * n, w)?;
        }
        r.skip_i32s(2, w)?;
    }
    Ok(if name.is_empty() { "STL plane".into() } else { name })
}

// ---------------------------------------------------------------------------
// dispatch
// ---------------------------------------------------------------------------

/// Read any supported legacy file, choosing the reader from the extension
/// (`.xfl`, `.fl5`, `.plr`) or, failing that, from the file header.
pub fn read_any(path: &Path) -> Result<LegacyImport, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("").to_ascii_lowercase();
    match ext.as_str() {
        "xfl" => read_xfl(&bytes),
        "fl5" => read_fl5(&bytes),
        "plr" => read_plr(&bytes),
        "wpa" => Err("old xflr5 .wpa projects are not supported; open the file in xflr5 v6 and save it as .xfl".into()),
        _ => {
            if bytes.len() < 4 {
                return Err(format!("{}: file too short", path.display()));
            }
            let be = i32::from_be_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            let le = i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
            if (200001..=210000).contains(&be) {
                read_xfl(&bytes)
            } else if (500001..=501000).contains(&be) {
                read_fl5(&bytes)
            } else if (100000..600000).contains(&le) {
                read_plr(&bytes)
            } else {
                Err(format!("{}: not an xflr5 / flow5 .xfl, .fl5 or .plr file", path.display()))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// QDataStream-style writer mirroring the C++ store code.
    struct W {
        b: Vec<u8>,
        le: bool,
        dbl: bool,
    }

    impl W {
        fn new(le: bool, dbl: bool) -> W {
            W { b: Vec::new(), le, dbl }
        }
        fn i(&mut self, v: i32) -> &mut Self {
            let a = if self.le { v.to_le_bytes() } else { v.to_be_bytes() };
            self.b.extend_from_slice(&a);
            self
        }
        fn is(&mut self, n: usize, v: i32) -> &mut Self {
            for _ in 0..n {
                self.i(v);
            }
            self
        }
        fn bo(&mut self, v: bool) -> &mut Self {
            self.b.push(v as u8);
            self
        }
        fn bos(&mut self, n: usize, v: bool) -> &mut Self {
            for _ in 0..n {
                self.bo(v);
            }
            self
        }
        fn f(&mut self, v: f64) -> &mut Self {
            let a = if self.le { v.to_le_bytes() } else { v.to_be_bytes() };
            self.b.extend_from_slice(&a);
            self
        }
        fn fs(&mut self, n: usize, v: f64) -> &mut Self {
            for _ in 0..n {
                self.f(v);
            }
            self
        }
        /// C++ float
        fn fl(&mut self, v: f64) -> &mut Self {
            if self.dbl {
                return self.f(v);
            }
            let a = if self.le { (v as f32).to_le_bytes() } else { (v as f32).to_be_bytes() };
            self.b.extend_from_slice(&a);
            self
        }
        fn s(&mut self, v: &str) -> &mut Self {
            let u: Vec<u16> = v.encode_utf16().collect();
            let n = (u.len() * 2) as u32;
            let a = if self.le { n.to_le_bytes() } else { n.to_be_bytes() };
            self.b.extend_from_slice(&a);
            for c in u {
                let a = if self.le { c.to_le_bytes() } else { c.to_be_bytes() };
                self.b.extend_from_slice(&a);
            }
            self
        }
        fn cs(&mut self, v: &str) -> &mut Self {
            self.b.push(v.len() as u8);
            self.b.extend_from_slice(v.as_bytes());
            self
        }
        fn color(&mut self) -> &mut Self {
            self.b.push(1);
            for c in [255u16, 10, 20, 30, 0] {
                let a = if self.le { c.to_le_bytes() } else { c.to_be_bytes() };
                self.b.extend_from_slice(&a);
            }
            self
        }
    }

    /// 21-point symmetric test foil in Selig order.
    fn test_foil() -> (Vec<f64>, Vec<f64>) {
        let mut x = Vec::new();
        let mut y = Vec::new();
        for k in 0..21 {
            let t = std::f64::consts::PI * k as f64 / 20.0;
            let xc = 0.5 * (1.0 + t.cos());
            let yt = 0.6 * (0.2969 * xc.sqrt() - 0.126 * xc - 0.3516 * xc * xc + 0.2843 * xc.powi(3) - 0.1036 * xc.powi(4)) * 0.2;
            x.push(xc);
            y.push(if k <= 10 { yt } else { -yt });
        }
        (x, y)
    }

    const PTS: [[f64; 4]; 3] = [[-2.0, 0.1, 0.009, 0.004], [0.0, 0.3, 0.008, 0.003], [2.0, 0.5, 0.0085, 0.0035]]; // alpha, cl, cd, cdp

    // ---------------- xflr5 .xfl writers ----------------

    fn x_style_xfl(w: &mut W) {
        w.i(0).i(1).i(0).color().bo(true);
    }
    fn x_style_fl5(w: &mut W) {
        w.i(500001).i(0).i(1).i(0).color().bo(true).s("");
    }
    fn x_foil(w: &mut W, name: &str, flap: bool) {
        let (x, y) = test_foil();
        w.i(100007).s(name).s("test foil");
        x_style_xfl(w);
        w.bo(false).bo(false).bo(flap);
        w.f(0.0).f(20.0).f(50.0);
        w.f(5.0).f(75.0).f(50.0);
        w.i(x.len() as i32);
        for k in 0..x.len() {
            w.f(x[k]).f(y[k]);
        }
    }
    fn x_polar(w: &mut W, foil: &str, name: &str) {
        w.i(100005).s(foil).s(name);
        x_style_xfl(w);
        w.i(1).i(1).i(1);
        w.f(200000.0).f(0.1).f(0.0).f(1.0).f(1.0).f(9.0);
        w.i(PTS.len() as i32);
        for p in PTS {
            w.fl(p[0]).fl(p[2]).fl(p[3]).fl(p[1]).fl(-0.05);
            w.fl(0.6).fl(0.9).fl(0.0).fl(-0.8).fl(200000.0).fl(0.25);
        }
        w.is(19, 0).fs(50, 0.0);
    }
    /// sections: (y, chord, offset, dihedral, twist, foil)
    fn x_wing(w: &mut W, name: &str, secs: &[(f64, f64, f64, f64, f64, &str)], mass: f64, ty: i32) {
        w.i(100001).s(name).s("").color().bo(true);
        w.i(secs.len() as i32);
        for s in secs {
            w.s(s.5).s(s.5).f(s.1).f(s.0).f(s.2).f(s.3).f(s.4).i(7).i(9).i(1).i(-2);
        }
        w.f(mass);
        w.i(1).f(0.05).f(0.1).f(0.2).f(0.0).s("servo");
        w.i(1).is(18, 0).i(ty).fs(50, 0.0);
    }
    fn x_body(w: &mut W) {
        w.i(100006).s("Pod").s("").color().i(2).i(0).i(19).i(11).f(0.0);
        w.i(3).is(3, 4);
        let frames = [(0.0, 0.01), (0.3, 0.05), (0.9, 0.02)];
        w.i(frames.len() as i32);
        for (x, r) in frames {
            w.i(5).f(x).i(1000).i(3);
            w.fl(x).fl(0.0).fl(r);
            w.fl(x).fl(r).fl(0.0);
            w.fl(x).fl(0.0).fl(-r);
        }
        w.f(0.3);
        w.i(1).f(0.2).f(0.1).f(0.0).f(0.0).s("battery");
        w.i(1).is(17, 0).i(3).i(3).fs(50, 0.0);
    }
    fn x_plane(w: &mut W) {
        w.i(100002).s("Glider").s("an xfl plane");
        x_style_fl5(w);
        x_wing(w, "Main", &[(0.0, 0.25, 0.0, 2.0, 0.0, "F1"), (1.0, 0.15, 0.05, 0.0, -1.0, "F1")], 0.5, 0);
        x_wing(w, "", &[(0.0, 0.2, 0.0, 0.0, 0.0, "F1"), (0.5, 0.1, 0.0, 0.0, 0.0, "F1")], 0.0, 1);
        x_wing(w, "Stab", &[(0.0, 0.12, 0.0, 0.0, 0.0, "F1"), (0.3, 0.08, 0.02, 0.0, 0.0, "F1")], 0.05, 2);
        x_wing(w, "Fin", &[(0.0, 0.15, 0.0, 0.0, 0.0, "F1"), (0.2, 0.1, 0.05, 0.0, 0.0, "F1")], 0.03, 3);
        // biplane, stab, fin, double fin, sym fin, -
        w.bo(false).bo(true).bo(true).bo(true).bo(false).bo(false);
        for (p, t) in [([0.0, 0.0, 0.0], 1.5), ([0.0, 0.0, 0.0], 0.0), ([0.9, 0.0, 0.02], -1.0), ([0.95, 0.15, 0.0], 0.0)] {
            w.f(p[0]).f(p[1]).f(p[2]).f(t);
        }
        w.bo(true).f(-0.2).f(-0.01).s("Pod");
        x_body(w);
        w.i(1).f(0.4).f(-0.1).f(0.0).f(0.0).s("ballast");
        w.is(20, 0).fs(50, 0.0);
    }
    fn x_wpolar(w: &mut W) {
        w.i(200014).s("Glider").s("T1").f(0.5).f(0.2).f(2.0);
        x_style_xfl(w);
        w.i(2).i(1).bos(7, false).f(0.0).f(1.225).f(1.5e-5).i(1).bo(true).f(1.0).f(0.1).f(0.0).f(0.0).fs(4, 0.0);
        w.i(1).f(0.0).i(1).f(100.0).f(1.1).f(10.0).f(0.0).f(0.0);
        w.i(2).fs(2 * 36, 0.5);
        w.is(19, 0).i(0).fs(35, 0.0).fs(4, 0.0).fs(4, 0.0).fs(7, 0.0);
    }
    fn x_popp(w: &mut W) {
        w.i(200002).s("Glider").s("T1");
        x_style_xfl(w);
        w.bo(false).bo(false).bo(true).bo(false).i(1).i(2).i(2).i(2).fs(5, 1.0);
        w.fl(0.1).fl(0.2).fl(0.3).fl(0.1).fl(0.2).fl(0.3);
        w.i(1);
        // WingOpp
        w.i(200001).s("Main").s("T1").i(2).bo(false).i(2).i(2).i(0).i(1).f(0.0).f(1.1).fs(19, 0.0);
        w.fs(2 * 22, 0.3).i(1).f(0.0).is(20, 0).fs(50, 0.0);
        w.i(0).i(0).i(0);
        w.fs(14 + 18, 0.0).i(0).fs(6 + 8 + 32 + 1 + 80, 0.0).is(20, 0).fs(2 + 4 + 44, 0.0);
    }
    fn xfl_bytes() -> Vec<u8> {
        let mut w = W::new(false, true);
        w.i(200002).is(6, 0);
        x_wpolar(&mut w);
        w.i(1);
        x_plane(&mut w);
        w.i(1);
        x_wpolar(&mut w);
        w.i(1);
        x_popp(&mut w);
        w.i(2);
        x_foil(&mut w, "F1", false);
        x_foil(&mut w, "F2 flap", true);
        w.i(1);
        x_polar(&mut w, "F1", "T1_Re0.200_M0.10_N9.0");
        w.b
    }

    #[test]
    fn xfl_project() {
        let imp = read_xfl(&xfl_bytes()).unwrap();
        assert_eq!(imp.foils.len(), 2);
        assert_eq!(imp.foils[0].name, "F1");
        assert_eq!(imp.foils[0].x.len(), 21);
        assert_eq!(imp.foils[0].hinge, None);
        let h = imp.foils[1].hinge.unwrap();
        assert!((h.0 - 0.75).abs() < 1e-12 && (h.1 - 0.5).abs() < 1e-12);
        assert!(imp.warnings.iter().any(|w| w.contains("flap deflection")));

        assert_eq!(imp.polars.len(), 1);
        let p = &imp.polars[0];
        assert_eq!((p.foil.as_str(), p.polar_type), ("F1", 1));
        assert!((p.re - 2e5).abs() < 1e-6 && (p.ncrit - 9.0).abs() < 1e-12);
        assert_eq!(p.points.len(), 3);
        assert!((p.points[1].cl - 0.3).abs() < 1e-12);
        assert!((p.points[1].cdf - 0.005).abs() < 1e-12);

        assert_eq!(imp.planes.len(), 1);
        let pl = &imp.planes[0];
        assert_eq!(pl.name, "Glider");
        // main + elevator + double fin (2)
        assert_eq!(pl.wings.len(), 4, "{:?}", pl.wings.iter().map(|w| &w.name).collect::<Vec<_>>());
        let m = &pl.wings[0];
        assert_eq!(m.kind, WingKind::Main);
        assert!((m.ry - 1.5).abs() < 1e-12);
        assert!((m.sections[1].y - 1.0).abs() < 1e-12 && (m.sections[0].dihedral - 2.0).abs() < 1e-12);
        assert_eq!((m.sections[0].nx, m.sections[0].ny), (7, 9));
        assert_eq!((m.sections[0].xdist, m.sections[0].ydist), (Spacing::Cosine, Spacing::Sine));
        assert!((m.mass - 0.5).abs() < 1e-12 && m.point_masses.len() == 1);
        let e = &pl.wings[1];
        assert_eq!(e.kind, WingKind::Elevator);
        assert!((e.position.x - 0.9).abs() < 1e-12 && (e.ry + 1.0).abs() < 1e-12);
        let (fl, fr) = (&pl.wings[2], &pl.wings[3]);
        assert_eq!(fl.kind, WingKind::Fin);
        assert!((fl.rx - 90.0).abs() < 1e-12 && !fl.two_sided && fl.close_inner);
        assert!((fl.position.y + 0.15).abs() < 1e-12 && (fr.position.y - 0.15).abs() < 1e-12);
        assert_eq!(pl.fuses.len(), 1);
        let f = &pl.fuses[0];
        assert!((f.position.x + 0.2).abs() < 1e-12 && (f.position.z + 0.01).abs() < 1e-12);
        let nb = f.nurbs.as_ref().unwrap();
        assert_eq!(nb.frames.len(), 3);
        assert_eq!(nb.frames[1].points.len(), 3);
        assert_eq!(f.frames.len(), 3);
        assert!((f.frames[1].half_width - 0.05).abs() < 1e-12 && (f.frames[1].half_height - 0.05).abs() < 1e-12);
        assert_eq!((f.nx, f.nh), (19, 11));
        assert!((f.mass - 0.3).abs() < 1e-12);
        // plane mass + body mass (offset by the body position)
        assert_eq!(pl.point_masses.len(), 2);
        let bm = pl.point_masses[1].position;
        assert!((bm.x + 0.1).abs() < 1e-12 && bm.y.abs() < 1e-12 && (bm.z + 0.01).abs() < 1e-12, "{bm:?}");
        pl.validate().unwrap();
    }

    #[test]
    fn xfl_format_200001() {
        let mut w = W::new(false, true);
        w.i(200001).is(6, 0);
        w.i(1).i(2).f(1.0).f(10.0).f(0.0).f(0.0).f(0.0).fl(1.225).fl(1.5e-5).f(0.0).f(0.0).bo(false).bo(false);
        w.i(0).i(0).i(0).i(1);
        x_foil(&mut w, "F1", false);
        w.i(0);
        let imp = read_xfl(&w.b).unwrap();
        assert_eq!(imp.foils.len(), 1);
    }

    // ---------------- flow5 .fl5 writers ----------------

    fn f_style(w: &mut W) {
        w.i(500756).i(0).i(1).i(0).color().bo(true).s("tag");
    }
    fn f_foil(w: &mut W, name: &str, flap: bool) {
        let (x, y) = test_foil();
        w.i(500753).s(name).s("fl5 foil").i(0).i(10).i(0).color().bo(true);
        w.bo(false).bo(false).bo(flap);
        w.f(0.0).f(0.2).f(0.5).f(0.0).f(0.7).f(0.5);
        w.f(0.0).i(3);
        w.i(x.len() as i32);
        for k in 0..x.len() {
            w.f(x[k]).f(y[k]);
        }
    }
    fn f_polar(w: &mut W, foil: &str, name: &str, ty: i32) {
        w.i(500750).s(foil).s(name);
        f_style(w);
        w.i(ty).i(0).i(1).i(1);
        w.f(300000.0).f(0.0).f(0.0).f(1.0).f(1.0).f(9.0);
        w.i(0).i(12).i(PTS.len() as i32);
        for (k, p) in PTS.iter().enumerate() {
            w.fl(p[0]).fl(p[2]).fl(p[3]).fl(p[1]).fl(-0.05).fl(0.0).fl(-0.8).fl(300000.0).fl(0.25).fl(k as f64 * 5.0);
            w.fl(0.6).fl(0.9).fl(0.5).fl(0.6).fl(1.0).fl(1.0);
        }
        w.i(0).i(1).f(0.0);
    }
    fn f_spline(w: &mut W) {
        w.i(500001);
        f_style(w);
        w.bo(false).i(79).i(2).f(0.0).f(0.0).f(1.0).f(0.0).i(2).f(1.0).f(1.0).bo(false).bo(false).f(0.0).i(0).i(0).i(0).i(0);
    }
    fn f_bspline(w: &mut W) {
        f_spline(w);
        w.i(3).i(0).i(0);
    }
    fn f_planepolar(w: &mut W) {
        // Polar3d v750
        w.i(500750).s("T1");
        f_style(w);
        w.i(5).bo(true).bo(true).bo(false).i(1).i(1).s("wheels").f(0.01).f(0.5);
        w.bo(false).bo(false).f(0.0).f(0.0).f(1.225).f(1.5e-5).f(0.0).f(0.0).bo(true).f(1.0).fs(3, 0.0).fs(6, 0.0);
        w.bo(true).bo(false).f(9.0).f(1.0).f(1.0).bo(false).i(1).f(100.0).f(1.1).bo(false).fs(4, 0.0).i(10).bo(false);
        w.bos(9, false).is(20, 0).fs(20, 0.0);
        // PlanePolar v750
        w.s("Plane").f(10.0).f(0.0).bo(false).bo(true).i(1).f(0.5).f(0.2).f(2.0);
        w.i(1).i(500750).s("flaps").i(2).f(0.0).f(5.0);
        w.i(0).bo(false);
        w.i(1).i(1).s("a").f(0.0).f(1.0);
        w.i(4);
        for n in ["Vinf", "alpha", "beta", "phi"] {
            w.s(n).f(0.0).f(1.0);
        }
        w.i(0).bo(true).bo(false).i(1).i(0);
        w.i(1).i(500750).fs(27, 0.1).fs(7, 0.2).fs(16, 0.0);
        w.bo(false);
        f_bspline(w);
        w.bo(true).bo(true).f(0.0).bos(10, false).is(20, 0).fs(20, 0.0);
    }
    fn f_inertia(w: &mut W, mass: f64, pm: Option<(&str, f64, [f64; 3])>) {
        w.i(500001).f(mass).fs(3, 0.0).fs(6, 0.0);
        match pm {
            Some((t, m, p)) => {
                w.i(1).f(m).f(p[0]).f(p[1]).f(p[2]).s(t);
            }
            None => {
                w.i(0);
            }
        }
        w.i(5).is(5, 0).i(5).fs(5, 0.0);
    }
    fn f_part(w: &mut W, name: &str, mass: f64, pm: Option<(&str, f64, [f64; 3])>) {
        w.i(500757).s(name).s("");
        f_style(w);
        f_inertia(w, mass, pm);
        w.bo(true).bo(false).f(0.01).f(0.1).i(23).f(0.01).f(0.1).i(23).i(0).i(1).f(1.0);
    }
    fn f_wing(w: &mut W, name: &str, secs: &[(f64, f64, f64, &str)], two_sided: bool, ty: i32, pos: [f64; 3], rot: [f64; 3]) {
        f_part(w, name, 0.4, Some(("lead", 0.01, [0.1, 0.0, 0.0])));
        w.i(500001).bo(true).i(secs.len() as i32);
        for s in secs {
            w.i(500001).s(s.3).s(s.3).f(s.1).f(s.0).f(s.2).f(1.0).f(0.0).i(13).i(11).i(4).i(6).i(0);
        }
        w.bo(two_sided).bo(!two_sided).bo(false).i(ty);
        w.f(pos[0]).f(pos[1]).f(pos[2]).f(rot[0]).f(rot[1]).f(rot[2]);
        w.i(2).i(1).i(0).i(0);
    }
    fn f_fuse_nurbs(w: &mut W) {
        f_part(w, "Body", 0.25, Some(("rx", 0.1, [0.3, 0.0, 0.0])));
        w.i(500003).f(-0.25).f(0.0).f(0.0).fs(3, 0.0);
        w.i(500001).bo(false).fs(6, 0.1).is(4, 0);
        w.f(0.05).i(0).i(0);
        w.i(500001).i(2);
        // NURBS
        w.i(500001).color().i(3).i(2).f(0.0).f(0.0).f(1.0).f(1.0).i(0).i(1).i(0);
        let frames = [(0.0, 0.01), (0.4, 0.05), (1.0, 0.01)];
        w.i(frames.len() as i32);
        for (x, r) in frames {
            w.i(500001).f(x).f(0.0).f(0.0).f(0.0).i(4);
            w.f(x).f(0.0).f(r).f(x).f(r * 0.9).f(r * 0.5).f(x).f(r).f(-r * 0.3).f(x).f(0.0).f(-r);
            w.i(0).i(0);
        }
        w.i(0).i(0);
        w.i(25).i(9);
        w.is(4, 3).is(3, 8).i(0).i(0);
    }
    fn f_panels(w: &mut W, n: usize) {
        w.i(500001).i(n as i32);
        for k in 0..n {
            w.fl(k as f64).fl(0.0).fl(0.0).fl(0.0).fl(1.0).fl(0.0).fl(0.0).fl(0.0).fl(1.0).bo(true);
        }
        w.i(0).i(0);
    }
    fn f_plane(w: &mut W) {
        w.i(500003).i(2).s(" Flow plane ").s("desc").i(0).i(1).i(0).color();
        f_wing(w, "Wing", &[(0.0, 0.3, 0.0, "N1"), (1.5, 0.2, 0.05, "N1")], true, 0, [0.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        f_wing(w, "Fin", &[(0.0, 0.2, 0.0, "N1"), (0.3, 0.1, 0.05, "N1")], false, 3, [1.0, 0.0, 0.0], [-90.0, 0.0, 0.0]);
        w.i(1).i(100005);
        f_fuse_nurbs(w);
        f_panels(w, 2);
        w.bo(true);
        f_inertia(w, 0.0, Some(("nose", 0.2, [-0.2, 0.0, 0.0])));
        w.i(1).i(1).i(0);
    }
    fn f_plane_stl(w: &mut W) {
        w.i(500005).s("STL thing").s("").color();
        f_style(w);
        w.bo(false).fs(3, 1.0).f(2.0).fs(3, 1.0);
        f_inertia(w, 1.0, None);
        // mesh
        w.i(500002).i(3);
        for k in 0..3 {
            w.i(k).fl(0.0).fl(0.0).fl(0.0).fl(0.0).fl(0.0).fl(1.0).bo(false).i(1).i(0).i(1).i(0);
        }
        w.i(1).i(0).i(0).i(1).i(2).i(-1).i(-1).i(-1).bo(true).bo(false).i(-1).i(4);
        w.i(500002).i(1).fs(9, 0.0);
        w.i(0).i(0);
    }
    fn f_opp(w: &mut W) {
        w.i(0);
        w.i(500750).s("N1").s("T1");
        f_style(w);
        w.i(0).f(1e5).f(0.0).f(2.0).f(0.0).bo(true).bo(true).fs(10, 0.1);
        w.i(1).i(0).bo(false).fl(1.0).fl(0.0).fl(0.0).fl(1.0);
        w.i(2).fs(4, 0.2).fs(4, 0.3);
        w.i(500001).i(10).i(10).i(1).i(1).i(0).fs(4, 0.0).fs(2, 0.0).i(0).i(0);
        w.i(0).i(0);
    }
    fn fl5_bytes(with_stl: bool) -> Vec<u8> {
        let mut w = W::new(false, true);
        w.i(500754);
        f_polar(&mut w, "", "", 1);
        f_planepolar(&mut w);
        w.s("");
        // spline foil
        w.i(500001);
        f_style(&mut w);
        f_bspline(&mut w);
        f_bspline(&mut w);
        w.bos(5, false).i(0).i(0);
        f_bspline(&mut w);
        f_bspline(&mut w);
        f_bspline(&mut w);
        f_spline(&mut w);
        w.i(0).i(0);
        w.i(0).i(0);
        // 2d objects
        w.i(500750).i(2);
        f_foil(&mut w, "N1", false);
        f_foil(&mut w, "N1 flap", true);
        w.i(2);
        f_polar(&mut w, "N1", "T1", 1);
        f_polar(&mut w, "N1 flap", "T6", 6);
        w.i(1);
        f_opp(&mut w);
        // planes
        w.i(if with_stl { 2 } else { 1 });
        if with_stl {
            w.i(1);
            f_plane_stl(&mut w);
        }
        w.i(0);
        f_plane(&mut w);
        w.b
    }

    #[test]
    fn fl5_project() {
        let imp = read_fl5(&fl5_bytes(true)).unwrap();
        assert_eq!(imp.foils.len(), 2);
        assert_eq!(imp.foils[1].hinge, Some((0.7, 0.5)));
        assert_eq!(imp.polars.len(), 2);
        assert_eq!(imp.polars[0].points.len(), 3);
        assert!((imp.polars[0].re - 3e5).abs() < 1e-6 && imp.polars[0].viscous);
        let t6 = &imp.polars[1];
        assert_eq!(t6.polar_type, 6);
        assert_eq!(t6.points.iter().map(|p| p.flap).collect::<Vec<_>>(), vec![0.0, 5.0, 10.0]);
        assert!(imp.warnings.iter().any(|w| w.contains("STL")));
        assert_eq!(imp.planes.len(), 1);
        let pl = &imp.planes[0];
        assert_eq!(pl.name, "Flow plane");
        assert_eq!(pl.wings.len(), 2);
        let (m, f) = (&pl.wings[0], &pl.wings[1]);
        assert!(m.two_sided && (m.ry - 2.0).abs() < 1e-12 && (m.mass - 0.4).abs() < 1e-12);
        assert_eq!(m.point_masses.len(), 1);
        assert_eq!((m.sections[0].nx, m.sections[0].ny, m.sections[0].xdist, m.sections[0].ydist), (13, 11, Spacing::Tanh, Spacing::Sine));
        assert_eq!(f.kind, WingKind::Fin);
        assert!(!f.two_sided && f.close_inner && (f.rx - 90.0).abs() < 1e-12);
        assert_eq!(pl.fuses.len(), 1);
        let fu = &pl.fuses[0];
        assert!((fu.position.x + 0.25).abs() < 1e-12 && (fu.mass - 0.25).abs() < 1e-12);
        let nb = fu.nurbs.as_ref().unwrap();
        assert_eq!((nb.frames.len(), nb.frames[0].points.len(), nb.degree_u, nb.degree_v), (3, 4, 3, 2));
        assert!((fu.frames[1].half_width - 0.05).abs() < 1e-12);
        assert_eq!((fu.nx, fu.nh), (25, 9));
        // fuse point mass + plane point mass
        assert_eq!(pl.point_masses.len(), 2);
        assert!((pl.point_masses[0].position.x - 0.05).abs() < 1e-12);
        pl.validate().unwrap();
    }

    // ---------------- .plr writers ----------------

    fn plr_v6_bytes() -> Vec<u8> {
        let (x, y) = test_foil();
        let mut w = W::new(true, false);
        w.i(100002).i(1);
        w.i(1006).cs("P1").cs("legacy").i(0).i(1).i(0x123456).i(1).i(0).i(0);
        w.i(0).fl(0.0).fl(20.0).fl(50.0);
        w.i(1).fl(3.0).fl(80.0).fl(40.0).fl(1.0).fl(1.0).fl(9.0);
        w.i(x.len() as i32);
        for k in 0..x.len() {
            w.fl(x[k]).fl(y[k]);
        }
        w.i(0);
        w.i(1);
        w.i(1005).cs("P1").cs("T2 polar").i(2).i(2).i(2).i(150000).fl(0.0).fl(0.0).i(PTS.len() as i32 + 1).fl(7.0).fl(1.0).fl(0.5);
        x_style_xfl(&mut w);
        for p in PTS.iter().chain(std::iter::once(&PTS[0])) {
            w.fl(p[0]).fl(p[2]).fl(p[3]).fl(p[1]).fl(0.0).fl(0.6).fl(0.9).fl(0.0).fl(-0.8).fl(150000.0).fl(0.25);
        }
        w.f(7.0).f(1.0).f(0.5);
        w.b
    }

    fn plr_v7_bytes() -> Vec<u8> {
        let mut w = W::new(true, false);
        w.i(500001).i(1);
        f_foil(&mut w, "Q1", true);
        w.i(1);
        f_polar(&mut w, "Q1", "T1", 1);
        w.b
    }

    #[test]
    fn plr_files() {
        let imp = read_plr(&plr_v6_bytes()).unwrap();
        assert_eq!(imp.foils.len(), 1);
        assert_eq!(imp.foils[0].note, "legacy");
        let h = imp.foils[0].hinge.unwrap();
        assert!((h.0 - 0.8).abs() < 1e-6 && (h.1 - 0.4).abs() < 1e-6);
        let p = &imp.polars[0];
        assert_eq!((p.polar_type, p.points.len()), (2, 3)); // duplicate alpha dropped
        assert!((p.ncrit - 7.0).abs() < 1e-12 && (p.xtr_bot - 0.5).abs() < 1e-12 && (p.re - 150000.0).abs() < 1e-9);
        assert!((p.points[2].cl - 0.5).abs() < 1e-6);
        assert!((p.points[2].mach - 0.0).abs() < 1e-12);

        let v7 = plr_v7_bytes();
        let imp = read_plr(&v7).unwrap();
        assert_eq!(imp.foils.len(), 1);
        assert_eq!(imp.foils[0].hinge, Some((0.7, 0.5)));
        assert_eq!(imp.polars[0].points.len(), 3);
        assert!(imp.warnings.is_empty());

        // flow5 on Windows writes .plr files in text mode: \n -> \r\n
        assert!(v7.contains(&b'\n'));
        let mut crlf = Vec::new();
        for &c in &v7 {
            if c == b'\n' {
                crlf.push(b'\r');
            }
            crlf.push(c);
        }
        let imp = read_plr(&crlf).unwrap();
        assert_eq!(imp.polars[0].points.len(), 3);
        assert!(!imp.warnings.is_empty());
    }

    #[test]
    fn read_any_dispatch() {
        let dir = std::env::temp_dir().join(format!("foiley_legacy_io_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cases: Vec<(&str, Vec<u8>)> = vec![("a.xfl", xfl_bytes()), ("b.fl5", fl5_bytes(false)), ("c.plr", plr_v6_bytes()), ("d.bin", fl5_bytes(false)), ("e", plr_v7_bytes())];
        for (name, b) in &cases {
            let p = dir.join(name);
            std::fs::write(&p, b).unwrap();
            let imp = read_any(&p).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(!imp.foils.is_empty(), "{name}");
        }
        assert!(read_any(&dir.join("missing.xfl")).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn malformed_input_never_panics() {
        type Reader = fn(&[u8]) -> Result<LegacyImport, String>;
        let streams: Vec<(Reader, Vec<u8>)> = vec![(read_xfl, xfl_bytes()), (read_fl5, fl5_bytes(true)), (read_plr, plr_v6_bytes()), (read_plr, plr_v7_bytes())];
        for (f, b) in &streams {
            assert!(f(b).is_ok());
            // every strict prefix is truncated data: must be an error, not a panic
            for len in 0..b.len() {
                assert!(f(&b[..len]).is_err(), "prefix of {len} / {} bytes accepted", b.len());
            }
            // single-byte corruptions: any result is fine, a panic is not
            let mut seed = 0x9E3779B97F4A7C15u64;
            let mut rnd = || {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                seed
            };
            for _ in 0..3000 {
                let mut c = b.clone();
                let i = (rnd() as usize) % c.len();
                c[i] = (rnd() & 0xFF) as u8;
                if rnd() % 4 == 0 {
                    let j = (rnd() as usize) % c.len();
                    c[j] = 0xFF;
                }
                let _ = f(&c);
            }
        }
        // pure garbage
        let mut seed = 12345u64;
        for len in [0usize, 1, 3, 4, 7, 64, 1000, 10000] {
            let g: Vec<u8> = (0..len)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    (seed >> 33) as u8
                })
                .collect();
            assert!(read_xfl(&g).is_err());
            assert!(read_fl5(&g).is_err());
            assert!(read_plr(&g).is_err());
        }
        // huge counts must be rejected before allocating
        let mut w = W::new(false, true);
        w.i(200002).is(6, 0);
        x_wpolar(&mut w);
        w.i(i32::MAX);
        assert!(read_xfl(&w.b).unwrap_err().contains("count"));
        let mut w = W::new(true, false);
        w.i(500001).i(1).i(500753).i(0x7FFF_FFF0);
        assert!(read_plr(&w.b).is_err());
    }
}
