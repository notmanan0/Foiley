//! Raw-bit dumps of stored points (finding F015).
//!
//! Foiley and the instrumented reference XFoil (`bench/xfoil-rawdump.patch`,
//! built by `build_ref.sh`) write the same records for every point PLRADD
//! stores, each value as the 16 hex digits of its IEEE-754 bits:
//!
//! ```text
//! POINT <nbl(1)> <nbl(2)>
//! ALFA | CL | CD | CDP | CM | XTR1 | XTR2  <Z16>
//! BL <is> <ibl> <UEDG> <THET> <DSTR> <CTAU>      for ibl = 2..nbl(is)
//! ```
//!
//! The writer runs only in builds with the `rawdump` cargo feature, and
//! only when the `RAWDUMP` environment variable names a file (as in the
//! Fortran).  Parsing and the comparison are always available.  Only the
//! double-precision reference build is comparable bit for bit.

use crate::xfoil::XFoil;

/// The scalar records of a point, in file order.
pub const SCALARS: [&str; 7] = ["ALFA", "CL", "CD", "CDP", "CM", "XTR1", "XTR2"];
/// The BL arrays of a `BL` record, in file order.
pub const ARRAYS: [&str; 4] = ["UEDG", "THET", "DSTR", "CTAU"];

/// One dumped point: raw bits of the scalars and BL stations.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Point {
    pub nbl: [usize; 2],
    pub scalars: [u64; 7],
    /// (side, station, [UEDG, THET, DSTR, CTAU])
    pub bl: Vec<(usize, usize, [u64; 4])>,
}

/// The records of the current point of `xf` (BL arrays only when `bl`).
pub fn format_point(xf: &XFoil, bl: bool) -> String {
    let mut s = String::new();
    let (n1, n2) = if bl { (xf.nbl[1], xf.nbl[2]) } else { (0, 0) };
    s.push_str(&format!("POINT{n1:6}{n2:6}\n"));
    for (name, v) in SCALARS.iter().zip([xf.alfa, xf.cl, xf.cd, xf.cdp, xf.cm, xf.xoctr[1], xf.xoctr[2]]) {
        s.push_str(&format!("{name} {:016X}\n", v.to_bits()));
    }
    if bl {
        for is in 1..=2 {
            for ibl in 2..=xf.nbl[is] {
                let v = [xf.uedg[is][ibl], xf.thet[is][ibl], xf.dstr[is][ibl], xf.ctau[is][ibl]];
                s.push_str(&format!("BL{is:6}{ibl:6} {:016X} {:016X} {:016X} {:016X}\n", v[0].to_bits(), v[1].to_bits(), v[2].to_bits(), v[3].to_bits()));
            }
        }
    }
    s
}

#[cfg(feature = "rawdump")]
mod sink {
    use std::io::Write;
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{Mutex, OnceLock};

    static SINK: OnceLock<Option<Mutex<std::fs::File>>> = OnceLock::new();
    static BL_LIVE: AtomicBool = AtomicBool::new(true);

    pub fn point(xf: &crate::xfoil::XFoil) {
        let sink = SINK.get_or_init(|| {
            let name = std::env::var("RAWDUMP").ok().filter(|n| !n.trim().is_empty())?;
            std::fs::File::create(name.trim()).ok().map(Mutex::new)
        });
        if let Some(f) = sink {
            let rec = super::format_point(xf, BL_LIVE.load(Ordering::Relaxed));
            let mut f = f.lock().unwrap_or_else(|e| e.into_inner());
            let _ = f.write_all(rec.as_bytes());
            let _ = f.flush();
        }
    }

    pub fn set_bl_live(live: bool) {
        BL_LIVE.store(live, Ordering::Relaxed);
    }
}

/// Hook in PLRADD: append the current point to the `RAWDUMP` file.
#[cfg(feature = "rawdump")]
pub fn point(xf: &XFoil) {
    sink::point(xf)
}

/// Whether the session state holds the live BL arrays of the point being
/// stored (false while the parallel sweep replays coefficients only).
#[cfg(feature = "rawdump")]
pub fn set_bl_live(live: bool) {
    sink::set_bl_live(live)
}

/// Whether this build carries the hook (the `rawdump` feature).
pub const ENABLED: bool = cfg!(feature = "rawdump");

/// Parse a dump (Foiley's or the reference's).
pub fn parse(text: &str) -> Result<Vec<Point>, String> {
    let hex = |t: &str, line: usize| u64::from_str_radix(t, 16).map_err(|e| format!("line {line}: bad hex '{t}': {e}"));
    let mut out: Vec<Point> = Vec::new();
    for (k, l) in text.lines().enumerate() {
        let t: Vec<&str> = l.split_whitespace().collect();
        let line = k + 1;
        match t.first().copied() {
            None => {}
            Some("POINT") if t.len() == 3 => {
                let n = |s: &str| s.parse::<usize>().map_err(|e| format!("line {line}: {e}"));
                out.push(Point { nbl: [n(t[1])?, n(t[2])?], ..Default::default() });
            }
            Some("BL") if t.len() == 7 => {
                let p = out.last_mut().ok_or(format!("line {line}: BL before POINT"))?;
                let n = |s: &str| s.parse::<usize>().map_err(|e| format!("line {line}: {e}"));
                p.bl.push((n(t[1])?, n(t[2])?, [hex(t[3], line)?, hex(t[4], line)?, hex(t[5], line)?, hex(t[6], line)?]));
            }
            Some(name) if t.len() == 2 && SCALARS.contains(&name) => {
                let p = out.last_mut().ok_or(format!("line {line}: {name} before POINT"))?;
                p.scalars[SCALARS.iter().position(|s| *s == name).unwrap()] = hex(t[1], line)?;
            }
            Some(_) => return Err(format!("line {line}: unexpected record '{l}'")),
        }
    }
    Ok(out)
}

/// Distance in units in the last place between two f64 bit patterns
/// (0 = identical bits; +0 and -0 are 1 apart; NaN against anything else is
/// `u64::MAX`).
pub fn ulp_distance(a: u64, b: u64) -> u64 {
    if a == b {
        return 0;
    }
    if f64::from_bits(a).is_nan() || f64::from_bits(b).is_nan() {
        return u64::MAX;
    }
    // map the bits to a monotone integer line
    let ord = |x: u64| -> i128 {
        if x >> 63 == 1 {
            -((x & !(1u64 << 63)) as i128) - 1
        } else {
            x as i128
        }
    };
    (ord(a) - ord(b)).unsigned_abs().min(u64::MAX as u128) as u64
}

/// Result of comparing two dumps.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Comparison {
    /// points in the two dumps
    pub points: (usize, usize),
    /// per quantity (`SCALARS` then `ARRAYS`): max ulp distance over the
    /// common points and stations
    pub max_ulp: Vec<(&'static str, u64)>,
    /// per quantity: max |a - b| over the largest |a| of that quantity
    pub max_scaled: Vec<(&'static str, f64)>,
    /// structural differences (point counts, BL station counts)
    pub structure: Vec<String>,
}

impl Comparison {
    pub fn identical(&self) -> bool {
        self.structure.is_empty() && self.points.0 > 0 && self.max_ulp.iter().all(|(_, u)| *u == 0)
    }
    /// "bit-identical (N points)" or the max ulp of each quantity that differs.
    pub fn summary(&self) -> String {
        if self.identical() {
            return format!("bit-identical ({} points)", self.points.0);
        }
        let mut parts: Vec<String> = self.max_ulp.iter().filter(|(_, u)| *u > 0).map(|(n, u)| if *u == u64::MAX { format!("{n} NaN") } else { format!("{n} {u}") }).collect();
        if parts.is_empty() {
            parts.push("all common values identical".into());
        }
        let mut s = format!("max ulp: {}", parts.join(", "));
        let worst = self.max_scaled.iter().fold(("", 0.0f64), |a, b| if b.1 > a.1 { *b } else { a });
        if worst.1 > 0.0 {
            s.push_str(&format!(" (largest scaled difference {:.1e} in {})", worst.1, worst.0));
        }
        if !self.structure.is_empty() {
            s.push_str(&format!("; {}", self.structure.join("; ")));
        }
        s
    }
}

/// Compare two dumps point by point (in order) and station by station.
pub fn compare(a: &[Point], b: &[Point]) -> Comparison {
    let mut max = [0u64; 11];
    let (mut dif, mut mag) = ([0.0f64; 11], [0.0f64; 11]);
    let mut acc = |i: usize, x: u64, y: u64| {
        let (u, v) = (f64::from_bits(x), f64::from_bits(y));
        max[i] = max[i].max(ulp_distance(x, y));
        dif[i] = dif[i].max((u - v).abs());
        mag[i] = mag[i].max(u.abs());
    };
    let mut structure = Vec::new();
    if a.len() != b.len() {
        structure.push(format!("{} vs {} points", a.len(), b.len()));
    }
    if a.is_empty() || b.is_empty() {
        structure.push("empty dump".into());
    }
    for (k, (p, q)) in a.iter().zip(b).enumerate() {
        for i in 0..7 {
            acc(i, p.scalars[i], q.scalars[i]);
        }
        if p.nbl != q.nbl || p.bl.len() != q.bl.len() {
            structure.push(format!("point {}: BL stations {:?} vs {:?}", k + 1, p.nbl, q.nbl));
        }
        for ((is, ib, x), (js, jb, y)) in p.bl.iter().zip(&q.bl) {
            if (is, ib) != (js, jb) {
                structure.push(format!("point {}: station ({is}, {ib}) vs ({js}, {jb})", k + 1));
                break;
            }
            for i in 0..4 {
                acc(7 + i, x[i], y[i]);
            }
        }
    }
    let names: Vec<&'static str> = SCALARS.iter().chain(ARRAYS.iter()).copied().collect();
    let scaled = (0..11).map(|i| if dif[i] == 0.0 { 0.0 } else { dif[i] / mag[i].max(f64::MIN_POSITIVE) });
    Comparison { points: (a.len(), b.len()), max_ulp: names.iter().copied().zip(max).collect(), max_scaled: names.iter().copied().zip(scaled).collect(), structure }
}
