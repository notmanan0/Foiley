//! ISO 10303-21 (STEP physical file) reader: the DATA section as a map of
//! entities, plus typed accessors for the geometric primitives (points,
//! directions, vectors, placements) with the file's length unit applied.

use super::surfaces::Frame;
use foiley_core::Vec3;
use std::collections::HashMap;

/// One parameter value of an entity.
#[derive(Clone, Debug)]
pub enum Val {
    Ref(u64),
    Num(f64),
    Str(String),
    Enum(String),
    List(Vec<Val>),
    Typed(String, Box<Val>),
    Unset,
}

impl Val {
    pub fn list(&self) -> Option<&[Val]> {
        match self {
            Val::List(v) => Some(v),
            _ => None,
        }
    }
    pub fn num(&self) -> Option<f64> {
        match self {
            Val::Num(x) => Some(*x),
            Val::Typed(_, v) => v.num(),
            _ => None,
        }
    }
    pub fn rf(&self) -> Option<u64> {
        match self {
            Val::Ref(r) => Some(*r),
            _ => None,
        }
    }
    /// `.T.` / `.F.` (None for anything else).
    pub fn boolean(&self) -> Option<bool> {
        match self {
            Val::Enum(e) if e == "T" => Some(true),
            Val::Enum(e) if e == "F" => Some(false),
            _ => None,
        }
    }
    pub fn enumeration(&self) -> Option<&str> {
        match self {
            Val::Enum(e) => Some(e),
            _ => None,
        }
    }
    pub fn string(&self) -> Option<&str> {
        match self {
            Val::Str(s) => Some(s),
            _ => None,
        }
    }
    /// All references in a list.
    pub fn refs(&self) -> Vec<u64> {
        self.list().map(|l| l.iter().filter_map(|v| v.rf()).collect()).unwrap_or_default()
    }
}

/// One entity instance: simple (one part) or complex (several parts).
#[derive(Clone, Debug)]
pub struct Ent {
    pub parts: Vec<(String, Vec<Val>)>,
}

impl Ent {
    pub fn part(&self, name: &str) -> Option<&Vec<Val>> {
        self.parts.iter().find(|p| p.0 == name).map(|p| &p.1)
    }
    /// The first (for simple instances, the only) type name.
    pub fn kind(&self) -> &str {
        self.parts.first().map(|p| p.0.as_str()).unwrap_or("")
    }
    pub fn has(&self, name: &str) -> bool {
        self.parts.iter().any(|p| p.0 == name)
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
            if self.i < self.s.len() && self.s[self.i..].starts_with(b"/*") {
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
                std::str::from_utf8(&self.s[st..self.i]).ok().and_then(|t| t.parse().ok()).map(Val::Ref).ok_or_else(|| format!("STEP: bad reference at byte {st}"))
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
                let e = String::from_utf8_lossy(&self.s[st..self.i]).to_uppercase();
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
            b'"' => {
                // binary literal: skip
                self.i += 1;
                while self.i < self.s.len() && self.s[self.i] != b'"' {
                    self.i += 1;
                }
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

/// A parsed STEP file.
pub struct StepFile {
    ents: HashMap<u64, Ent>,
    /// metres per file length unit
    pub length_scale: f64,
    /// radians per file plane-angle unit
    pub angle_scale: f64,
}

impl StepFile {
    pub fn parse(text: &str) -> Result<StepFile, String> {
        let start = text.find("DATA;").ok_or("not a STEP file (no DATA section)")? + 5;
        let mut lx = Lexer { s: text.as_bytes(), i: start };
        let mut map = HashMap::new();
        loop {
            lx.ws();
            if lx.i >= lx.s.len() || lx.s[lx.i..].starts_with(b"ENDSEC") {
                break;
            }
            if lx.peek() != b'#' {
                // skip to the next statement
                while lx.i < lx.s.len() && lx.s[lx.i] != b';' {
                    lx.i += 1;
                }
                lx.i += 1;
                continue;
            }
            let id = match lx.value()? {
                Val::Ref(r) => r,
                _ => return Err(format!("STEP: bad entity id at byte {}", lx.i)),
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
        let mut f = StepFile { ents: map, length_scale: 1.0, angle_scale: 1.0 };
        f.length_scale = f.unit_scale("LENGTH_UNIT").unwrap_or(1.0);
        f.angle_scale = f.unit_scale("PLANE_ANGLE_UNIT").unwrap_or(1.0);
        Ok(f)
    }

    pub fn get(&self, id: u64) -> Option<&Ent> {
        self.ents.get(&id)
    }
    pub fn len(&self) -> usize {
        self.ents.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ents.is_empty()
    }
    /// Entity ids in increasing order.
    pub fn ids(&self) -> Vec<u64> {
        let mut v: Vec<u64> = self.ents.keys().copied().collect();
        v.sort_unstable();
        v
    }
    /// Ids of every entity that has a part `name`, in increasing order.
    pub fn ids_of(&self, name: &str) -> Vec<u64> {
        let mut v: Vec<u64> = self.ents.iter().filter(|(_, e)| e.has(name)).map(|(k, _)| *k).collect();
        v.sort_unstable();
        v
    }
    /// The arguments of part `name` of entity `id`.
    pub fn part(&self, id: u64, name: &str) -> Option<&Vec<Val>> {
        self.get(id)?.part(name)
    }
    pub fn kind(&self, id: u64) -> &str {
        self.get(id).map(|e| e.kind()).unwrap_or("")
    }

    /// CARTESIAN_POINT in metres (missing coordinates are 0).
    pub fn point(&self, id: u64) -> Option<Vec3> {
        let c = self.part(id, "CARTESIAN_POINT")?.get(1)?.list()?;
        let s = self.length_scale;
        Some(Vec3::new(c.first()?.num()?, c.get(1).and_then(|v| v.num()).unwrap_or(0.0), c.get(2).and_then(|v| v.num()).unwrap_or(0.0)) * s)
    }
    /// CARTESIAN_POINT as raw coordinates (parameter space: no unit applied).
    pub fn point_raw(&self, id: u64) -> Option<Vec3> {
        let c = self.part(id, "CARTESIAN_POINT")?.get(1)?.list()?;
        Some(Vec3::new(c.first()?.num()?, c.get(1).and_then(|v| v.num()).unwrap_or(0.0), c.get(2).and_then(|v| v.num()).unwrap_or(0.0)))
    }
    /// DIRECTION, normalised.
    pub fn direction(&self, id: u64) -> Option<Vec3> {
        let c = self.part(id, "DIRECTION")?.get(1)?.list()?;
        let d = Vec3::new(c.first()?.num()?, c.get(1).and_then(|v| v.num()).unwrap_or(0.0), c.get(2).and_then(|v| v.num()).unwrap_or(0.0));
        let n = d.norm();
        (n > 0.0 && n.is_finite()).then(|| d / n)
    }
    /// VECTOR (direction times magnitude); `length` applies the length unit.
    pub fn vector(&self, id: u64, length: bool) -> Option<Vec3> {
        let a = self.part(id, "VECTOR")?;
        let d = self.direction(a.get(1)?.rf()?)?;
        let m = a.get(2)?.num()?;
        Some(d * (m * if length { self.length_scale } else { 1.0 }))
    }
    /// AXIS2_PLACEMENT_3D (or _2D) as an orthonormal frame.
    pub fn axis2(&self, id: u64) -> Option<Frame> {
        let e = self.get(id)?;
        if let Some(a) = e.part("AXIS2_PLACEMENT_3D") {
            let o = self.point(a.get(1)?.rf()?)?;
            let z = a.get(2).and_then(|v| v.rf()).and_then(|r| self.direction(r)).unwrap_or(Vec3::Z);
            let x = a.get(3).and_then(|v| v.rf()).and_then(|r| self.direction(r));
            return Some(Frame::new(o, z, x));
        }
        if let Some(a) = e.part("AXIS2_PLACEMENT_2D") {
            let o = self.point_raw(a.get(1)?.rf()?)?;
            let x = a.get(2).and_then(|v| v.rf()).and_then(|r| self.direction(r));
            return Some(Frame::new(o, Vec3::Z, x));
        }
        None
    }
    /// AXIS1_PLACEMENT: (location, unit axis).
    pub fn axis1(&self, id: u64) -> Option<(Vec3, Vec3)> {
        let a = self.part(id, "AXIS1_PLACEMENT")?;
        let o = self.point(a.get(1)?.rf()?)?;
        let z = a.get(2).and_then(|v| v.rf()).and_then(|r| self.direction(r)).unwrap_or(Vec3::Z);
        Some((o, z))
    }

    /// Scale of the unit of kind `what` (LENGTH_UNIT / PLANE_ANGLE_UNIT)
    /// declared in a GLOBAL_UNIT_ASSIGNED_CONTEXT.
    fn unit_scale(&self, what: &str) -> Option<f64> {
        for id in self.ids_of("GLOBAL_UNIT_ASSIGNED_CONTEXT") {
            let g = self.part(id, "GLOBAL_UNIT_ASSIGNED_CONTEXT")?;
            for u in g.first().map(|v| v.refs()).unwrap_or_default() {
                if self.get(u).map(|e| e.has(what)).unwrap_or(false) {
                    if let Some(s) = self.unit_of(u, 0) {
                        return Some(s);
                    }
                }
            }
        }
        None
    }
    fn unit_of(&self, id: u64, depth: usize) -> Option<f64> {
        if depth > 6 {
            return None;
        }
        let e = self.get(id)?;
        if let Some(si) = e.part("SI_UNIT") {
            let pre = match si.first().and_then(|v| v.enumeration()) {
                Some("EXA") => 1e18,
                Some("PETA") => 1e15,
                Some("TERA") => 1e12,
                Some("GIGA") => 1e9,
                Some("MEGA") => 1e6,
                Some("KILO") => 1e3,
                Some("HECTO") => 1e2,
                Some("DECA") => 1e1,
                Some("DECI") => 1e-1,
                Some("CENTI") => 1e-2,
                Some("MILLI") => 1e-3,
                Some("MICRO") => 1e-6,
                Some("NANO") => 1e-9,
                _ => 1.0,
            };
            return Some(pre);
        }
        if let Some(cb) = e.part("CONVERSION_BASED_UNIT") {
            let m = self.get(cb.get(1)?.rf()?)?;
            let a = m.parts.iter().find(|p| p.0.ends_with("MEASURE_WITH_UNIT")).map(|p| &p.1)?;
            let v = a.first()?.num()?;
            let s = self.unit_of(a.get(1)?.rf()?, depth + 1).unwrap_or(1.0);
            return Some(v * s);
        }
        None
    }
}
