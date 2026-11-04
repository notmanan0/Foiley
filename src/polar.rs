//! Polar storage and XFoil polar files: POLWRIT / POLREAD / POLREF
//! (iopol.f), PINDEX.INC column definitions, HSORT, and Fortran
//! unformatted-sequential record helpers for the polar dump file.

use crate::fmt::{e_fmt, f_fmt, fixed};

// ---- PINDEX.INC ---------------------------------------------------------
pub const IAL: usize = 1;
pub const ICL: usize = 2;
pub const ICD: usize = 3;
pub const ICM: usize = 4;
pub const ICW: usize = 5;
pub const ICV: usize = 6;
pub const ICP: usize = 7;
pub const IMA: usize = 8;
pub const IRE: usize = 9;
pub const INC: usize = 10;
pub const ICH: usize = 11;
pub const IMC: usize = 12;
pub const ICDH: usize = 13;
pub const ICMDOT: usize = 14;
pub const IPTOT: usize = 14;

pub const JTP: usize = 1;
pub const JTN: usize = 2;
pub const JPTOT: usize = 2;

/// array limits from XFOIL.INC
pub const NAX: usize = 800;
pub const NPX: usize = 12;
pub const NFX: usize = 128;

pub const CPOLNAME: [&str; IPTOT + 1] =
    ["", "alpha", "CL", "CD", "CM", "CDw", "CDv", "CDp", "Mach", "Re", "Ncrit", "Chinge", "Cpmin", "CDh", "Cmdot"];

/// (kind, width, decimals) of CPOLFORM; kind 'F' or 'E'
pub const CPOLFORM: [(char, usize, usize); IPTOT + 1] = [
    ('F', 0, 0),
    ('F', 7, 3),
    ('F', 9, 4),
    ('F', 10, 5),
    ('F', 9, 4),
    ('F', 10, 5),
    ('F', 10, 5),
    ('F', 10, 5),
    ('F', 8, 4),
    ('E', 11, 3),
    ('F', 7, 3),
    ('F', 10, 5),
    ('F', 9, 4),
    ('F', 11, 5),
    ('F', 9, 5),
];

pub const CPOLSNAME: [&str; JPTOT + 1] = ["", "Xtrip", "Xtr"];
pub const CPOLSFORM: [(char, usize, usize); JPTOT + 1] = [('F', 0, 0), ('F', 9, 4), ('F', 9, 4)];

pub const VERSION: f64 = 6.97;

fn fmt_val(kind: char, w: usize, d: usize, v: f64) -> String {
    if kind == 'E' {
        e_fmt(v, w, d)
    } else {
        f_fmt(v, w, d)
    }
}

/// One converged operating point (used by the sequence solver, sweep.rs).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PolarPoint {
    pub alpha: f64,
    pub cl: f64,
    pub cd: f64,
    pub cdp: f64,
    pub cm: f64,
    pub cdv: f64,
    pub mach: f64,
    pub re: f64,
    pub ncrit: f64,
    pub xtr_top: f64,
    pub xtr_bot: f64,
    pub xtrip_top: f64,
    pub xtrip_bot: f64,
    pub chinge: f64,
    pub cpmin: f64,
    pub converged: bool,
}

/// One stored polar (CPOL, CPOLSD, CPOLXY and per-polar parameters).
#[derive(Clone, Debug)]
pub struct Polar {
    /// NAMEPOL (CHARACTER*48, stripped)
    pub name: String,
    /// CODEPOL (as read from a file, may carry trailing blanks)
    pub code: String,
    /// VERSPOL
    pub vers: f64,
    /// PFNAME, PFNAMX
    pub pfname: String,
    pub pfnamx: String,
    pub reyn: f64,
    pub mach: f64,
    pub acrit: f64,
    /// XSTRIPP(1..2) (index 0 unused)
    pub xstrip: [f64; 3],
    pub ptrat: f64,
    pub etap: f64,
    pub retyp: i32,
    pub matyp: i32,
    /// CPOL(ia, k): row ia-1, column k (1..IPTOT)
    pub cpol: Vec<[f64; IPTOT + 1]>,
    /// CPOLSD(ia, is, jp): row ia-1, [jp][is]
    pub cpolsd: Vec<[[f64; 3]; JPTOT + 1]>,
    /// NXYPOL and CPOLXY(.,1..2)  (1-based)
    pub nxy: usize,
    pub xyx: Vec<f64>,
    pub xyy: Vec<f64>,
}

impl Default for Polar {
    fn default() -> Self {
        Polar {
            name: String::new(),
            code: String::new(),
            vers: 0.0,
            pfname: String::new(),
            pfnamx: String::new(),
            reyn: 0.0,
            mach: 0.0,
            acrit: 0.0,
            xstrip: [0.0; 3],
            ptrat: 0.0,
            etap: 0.0,
            retyp: 0,
            matyp: 0,
            cpol: Vec::new(),
            cpolsd: Vec::new(),
            nxy: 0,
            xyx: vec![0.0; crate::xfoil::IQX + 1],
            xyy: vec![0.0; crate::xfoil::IQX + 1],
        }
    }
}

impl Polar {
    /// NAPOL
    pub fn na(&self) -> usize {
        self.cpol.len()
    }
}

/// Reference polar (POLREF data).
#[derive(Clone, Debug, Default)]
pub struct RefPolar {
    pub name: String,
    /// NDREF(k), k = 1..4  (index 0 unused)
    pub nd: [usize; 5],
    /// CPOLREF(i, 1..2, k): data[k][i-1] = (x, y)
    pub data: [Vec<(f64, f64)>; 5],
}

/// POLWRIT: returns the text of the polar file / listing.
#[allow(clippy::too_many_arguments)]
pub fn polwrit(
    pol: &Polar,
    lhead: bool,
    ia1: i64,
    ia2: i64,
    ipol: &[usize],
    nipol: usize,
    jpol: &[usize],
    njpol: usize,
    code: &str,
    version: f64,
) -> String {
    let nbl = 1usize;
    let mut out = String::new();
    if lhead {
        out.push_str("  \n");
        out.push_str(&format!("       {}         Version{}\n", code, f_fmt(version, 5, 2)));
        out.push_str("  \n");
        out.push_str(&format!(" Calculated polar for: {}\n", fixed(&pol.name, 48)));
        out.push_str("  \n");
        let line1 = match pol.retyp {
            1 => " Reynolds number fixed       ",
            2 => " Reynolds number ~ 1/sqrt(CL)",
            3 => " Reynolds number ~ 1/CL      ",
            _ => "                             ",
        };
        let line2 = match pol.matyp {
            1 => "   Mach number fixed         ",
            2 => "   Mach number ~ 1/sqrt(CL)  ",
            3 => "   Mach number ~ 1/CL        ",
            _ => "                             ",
        };
        out.push_str(&format!(" {}{}{}{}\n", i_fmt(pol.retyp as i64, 1), i_fmt(pol.matyp as i64, 2), line1, line2));
        out.push_str("  \n");
        out.push_str(&format!(" xtrf = {} (top)    {} (bottom)  \n", f_fmt(pol.xstrip[1], 7, 3), f_fmt(pol.xstrip[2], 9, 3)));
        out.push_str(&format!(
            " Mach = {}     Re = {} e 6     Ncrit = {}\n",
            f_fmt(pol.mach, 7, 3),
            f_fmt(pol.reyn / 1.0e6, 9, 3),
            f_fmt(pol.acrit, 7, 3)
        ));
        if pol.ptrat != 0.0 {
            out.push_str(&format!(" pi_p = {}     eta_p = {}\n", f_fmt(pol.ptrat, 7, 4), f_fmt(pol.etap, 9, 4)));
        }
        out.push_str("  \n");

        let mut linel = vec![b' '; 256];
        let mut lined = vec![b' '; 256];
        let put = |buf: &mut Vec<u8>, a: usize, s: &[u8]| {
            for (i, &ch) in s.iter().enumerate() {
                if a - 1 + i < buf.len() {
                    buf[a - 1 + i] = ch;
                }
            }
        };
        let dashes = |n: usize| vec![b'-'; n.min(32)];
        let mut kl = 1usize;
        let mut kd = 1usize;
        for &ip in ipol.iter().take(nipol) {
            if ip == 0 || ip > IPTOT {
                continue;
            }
            let nform = CPOLFORM[ip].1;
            let nm = CPOLNAME[ip].as_bytes();
            let nname = nm.len();
            let nblank = ((nform as i64 - nname as i64 + 2) / 2).max(0) as usize;
            put(&mut linel, kl + 1 + nblank, nm);
            kl += nform;
            put(&mut lined, kd + 2, &dashes(nform - 1));
            kd += nform;
        }
        for &jp in jpol.iter().take(njpol) {
            if jp == 0 || jp > JPTOT {
                continue;
            }
            let nform = CPOLSFORM[jp].1;
            let nm = CPOLSNAME[jp];
            let nname = nm.len();
            let nblank = ((nform as i64 - nname as i64 - 2) / 2).max(0) as usize;
            for _n in 1..=nbl {
                put(&mut linel, kl + 1 + nblank, format!("Top_{}", nm).as_bytes());
                kl += nform;
                put(&mut lined, kd + 2, &dashes(nform - 1));
                kd += nform;
                put(&mut linel, kl + 1 + nblank, format!("Bot_{}", nm).as_bytes());
                kl += nform;
                put(&mut lined, kd + 2, &dashes(nform - 1));
                kd += nform;
            }
        }
        // LINEL / LINED are CHARACTER*128
        let kl = kl.min(128);
        let kd = kd.min(128);
        out.push_str(std::str::from_utf8(&linel[..kl]).unwrap());
        out.push('\n');
        out.push_str(std::str::from_utf8(&lined[..kd]).unwrap());
        out.push('\n');
    }

    if ia2 >= ia1 && ia1 >= 1 {
        for ia in ia1..=ia2 {
            let ia = ia as usize;
            if ia > pol.cpol.len() {
                break;
            }
            let row = &pol.cpol[ia - 1];
            let mut s = String::from(" ");
            for &ip in ipol.iter().take(nipol) {
                if ip == 0 || ip > IPTOT {
                    continue;
                }
                let (k, w, d) = CPOLFORM[ip];
                s += &fmt_val(k, w, d, row[ip]);
            }
            for &jp in jpol.iter().take(njpol) {
                if jp == 0 || jp > JPTOT {
                    continue;
                }
                let (k, w, d) = CPOLSFORM[jp];
                for is in 1..=2 * nbl {
                    s += &fmt_val(k, w, d, pol.cpolsd[ia - 1][jp][is]);
                }
            }
            s.push('\n');
            out.push_str(&s);
        }
    }
    out
}

/// Fortran Iw
pub fn i_fmt(v: i64, w: usize) -> String {
    let s = v.to_string();
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// Fortran INDEX (1-based, 0 if not found)
pub fn index(s: &str, pat: &str) -> usize {
    match s.find(pat) {
        Some(k) => k + 1,
        None => 0,
    }
}

/// Fortran substring LINE(a:b) of a CHARACTER*len variable (1-based, blank padded)
pub fn substr(s: &str, a: usize, b: usize) -> String {
    let bytes = s.as_bytes();
    let mut out = String::new();
    if a == 0 || b < a {
        return out;
    }
    for k in a..=b {
        out.push(if k - 1 < bytes.len() { bytes[k - 1] as char } else { ' ' });
    }
    out
}

/// Parse one Fortran list-directed real (first item of string).
pub fn list_real1(s: &str) -> Option<f64> {
    let tok = s.split([' ', ',', '\t', '/']).find(|t| !t.is_empty())?;
    parse_fortran_real(tok)
}

pub fn parse_fortran_real(tok: &str) -> Option<f64> {
    let t = tok.replace(['d', 'D'], "e");
    let t = t.trim();
    if t.is_empty() {
        return None;
    }
    if let Ok(v) = t.parse::<f64>() {
        return Some(v);
    }
    // Fortran also accepts exponent without letter, e.g. "1.0+05"
    let b = t.as_bytes();
    for k in (1..b.len()).rev() {
        if (b[k] == b'+' || b[k] == b'-') && b[k - 1] != b'e' && b[k - 1] != b'E' {
            let s2 = format!("{}e{}", &t[..k], &t[k..]);
            return s2.parse::<f64>().ok();
        }
    }
    None
}

/// Parse one Fortran list-directed integer.
pub fn list_int1(s: &str) -> Option<i64> {
    let tok = s.split([' ', ',', '\t', '/']).find(|t| !t.is_empty())?;
    tok.trim().parse::<i64>().ok()
}

/// POLREAD: reads polar save file text into `pol` (in place, like the
/// Fortran which overwrites only the quantities found). Returns ERROR.
pub fn polread(txt: &str, pol: &mut Polar) -> bool {
    let mut lhead = true;
    pol.cpol.clear();
    pol.cpolsd.clear();
    let mut nbl = 1usize;
    let mut nipol = 0usize;
    let mut ipol = [0usize; IPTOT + 1];
    let mut ispol = [[0usize; JPTOT + 1]; 3];
    pol.retyp = 0;
    pol.matyp = 0;
    let mut kf = 0usize;
    let mut rinp = [0.0f64; IPTOT + 2 * JPTOT + 1];

    for raw in txt.lines() {
        // CHARACTER*128 LINE
        let line: String = raw.chars().take(128).collect();
        if line.trim_matches(' ').is_empty() {
            continue;
        }
        if lhead {
            let mut ldlab = true;
            let k = index(&line, "Version");
            if k != 0 {
                let k1 = line.bytes().position(|c| c != b' ').map(|p| p + 1).unwrap_or(129);
                if k > k1 {
                    pol.code = fixed(&substr(&line, k1, k - 1), 48);
                    if let Some(v) = list_real1(&substr(&line, k + 7, 128)) {
                        pol.vers = v;
                    }
                }
                ldlab = false;
            }
            let kfn = index(&line, "for:");
            if kfn != 0 {
                kf = kfn;
                pol.name = substr(&line, kf + 5, 128).chars().take(48).collect::<String>().trim().to_string();
                ldlab = false;
            }
            let ke = index(&line, "elements");
            if ke > 0 {
                if ke >= 5 {
                    if let Some(v) = list_int1(&substr(&line, ke - 4, ke - 1)) {
                        nbl = v.max(1) as usize;
                        pol.name = substr(&line, kf + 5, ke - 4).chars().take(48).collect::<String>().trim().to_string();
                        if 2 * nbl > 2 {
                            nbl = 1;
                            crate::outln!(" POLREAD: Number of elements set to array limit {:11}", nbl);
                        }
                    }
                }
                ldlab = false;
            }
            let kr = index(&line, "Reynolds number");
            let km = index(&line, "Mach number");
            if kr != 0 {
                let kend = if km > kr { km - 1 } else { 128 };
                let sub = substr(&line, kr, kend);
                if sub.contains("fixed") {
                    pol.retyp = 1;
                } else if sub.contains("1/sqrt(CL)") {
                    pol.retyp = 2;
                } else if sub.contains("1/CL") {
                    pol.retyp = 3;
                }
                ldlab = false;
            }
            if km != 0 {
                let kend = if kr > km { kr - 1 } else { 128 };
                let sub = substr(&line, km, kend);
                if sub.contains("fixed") {
                    pol.matyp = 1;
                } else if sub.contains("1/sqrt(CL)") {
                    pol.matyp = 2;
                } else if sub.contains("1/CL") {
                    pol.matyp = 3;
                }
                ldlab = false;
            }
            let k = index(&line, "xtrf");
            if k != 0 {
                let kt = index(&line, "(top)");
                let kb = index(&line, "(bottom)");
                let kel = index(&line, "element ");
                let ks = index(&line, "(suc");
                let kp = index(&line, "(pre");
                let nel = if kel != 0 { list_int1(&substr(&line, kel + 7, kel + 12)) } else { Some(1) };
                if let Some(nel) = nel {
                    if nel as usize <= nbl && nel >= 1 {
                        pol.xstrip[1] = 1.0;
                        pol.xstrip[2] = 1.0;
                        'tr: {
                            if kt > 0 {
                                match list_real1(&substr(&line, k + 6, kt - 1)) {
                                    Some(v) => pol.xstrip[1] = v,
                                    None => break 'tr,
                                }
                            }
                            if kb > kt {
                                match list_real1(&substr(&line, kt + 5, kb - 1)) {
                                    Some(v) => pol.xstrip[2] = v,
                                    None => break 'tr,
                                }
                            }
                            if ks > 0 {
                                match list_real1(&substr(&line, k + 6, ks - 1)) {
                                    Some(v) => pol.xstrip[1] = v,
                                    None => break 'tr,
                                }
                            }
                            if kp > ks {
                                if let Some(v) = list_real1(&substr(&line, ks + 5, kp - 1)) {
                                    pol.xstrip[2] = v;
                                }
                            }
                        }
                    }
                }
                ldlab = false;
            }
            let k = index(&line, "Mach =");
            if k != 0 {
                if let Some(v) = list_real1(&substr(&line, k + 6, 128)) {
                    pol.mach = v;
                }
                ldlab = false;
            }
            let k = index(&line, "Re =");
            if k != 0 {
                if let Some(v) = list_real1(&substr(&line, k + 4, 128)) {
                    pol.reyn = v * 1.0e6;
                }
                ldlab = false;
            }
            let k = index(&line, "Ncrit =");
            if k != 0 {
                if let Some(v) = list_real1(&substr(&line, k + 7, 128)) {
                    pol.acrit = v;
                }
                ldlab = false;
            }
            let k = index(&line, "pi_p =");
            if k != 0 {
                if let Some(v) = list_real1(&substr(&line, k + 6, 128)) {
                    pol.ptrat = v;
                }
                ldlab = false;
            }
            let k = index(&line, "eta_p =");
            if k != 0 {
                if let Some(v) = list_real1(&substr(&line, k + 7, 128)) {
                    pol.etap = v;
                }
                ldlab = false;
            }
            if ldlab && nipol == 0 {
                let nt = IPTOT + 2 * JPTOT;
                let mut itmp0 = vec![0usize; nt + 1];
                for ip in 1..=IPTOT {
                    itmp0[ip] = index(&line, CPOLNAME[ip]);
                }
                for jp in 1..=JPTOT {
                    let k1 = index(&line, &format!("Top {}", CPOLSNAME[jp]));
                    let k2 = index(&line, &format!("Top_{}", CPOLSNAME[jp]));
                    itmp0[IPTOT + jp] = k1.max(k2);
                    let k1 = index(&line, &format!("Bot {}", CPOLSNAME[jp]));
                    let k2 = index(&line, &format!("Bot_{}", CPOLSNAME[jp]));
                    itmp0[IPTOT + jp + JPTOT] = k1.max(k2);
                }
                let mut itmp = itmp0.clone();
                for _ in 1..=nt {
                    for ip in 1..nt {
                        if itmp[ip] > itmp[ip + 1] {
                            itmp.swap(ip, ip + 1);
                        }
                    }
                }
                for ipt in 1..=nt {
                    if itmp[ipt] > 0 {
                        nipol += 1;
                        for ip in 1..=IPTOT {
                            if itmp[ipt] == itmp0[ip] {
                                ipol[ip] = nipol;
                            }
                        }
                        for jp in 1..=JPTOT {
                            if itmp[ipt] == itmp0[IPTOT + jp] {
                                ispol[1][jp] = nipol;
                            }
                            if itmp[ipt] == itmp0[IPTOT + JPTOT + jp] {
                                ispol[2][jp] = nipol;
                            }
                        }
                    }
                }
            }
            if line.contains("-----") {
                lhead = false;
            }
        } else {
            if pol.cpol.len() >= NAX {
                continue;
            }
            let ninp = IPTOT + 2 * JPTOT;
            let v = match crate::geom::getflt(&line, ninp) {
                Some(v) => v,
                None => return true,
            };
            // (RINP keeps values of the previous line beyond those read)
            for (k, &x) in v.iter().enumerate() {
                rinp[k + 1] = x;
            }
            let mut row = [0.0f64; IPTOT + 1];
            for ip in 1..=IPTOT {
                row[ip] = rinp[ipol[ip]];
            }
            let mut sd = [[0.0f64; 3]; JPTOT + 1];
            for jp in 1..=JPTOT {
                sd[jp][1] = rinp[ispol[1][jp]];
                sd[jp][2] = rinp[ispol[2][jp]];
            }
            let acl = row[ICL].max(0.001);
            // (the Fortran tests IPOL(KP)/ISPOL(1,KP) against the parameter
            //  indices; for XFoil-written files Re/Mach/Ncrit/trip come from
            //  the header)
            let mut lire = false;
            let mut lima = false;
            let mut linc = false;
            let mut ljtp = false;
            for kp in 1..=nipol.min(IPTOT) {
                if ipol[kp] == IRE {
                    lire = true;
                }
                if ipol[kp] == IMA {
                    lima = true;
                }
                if ipol[kp] == INC {
                    linc = true;
                }
                if kp <= JPTOT && ispol[1][kp] == JTP {
                    ljtp = true;
                }
            }
            if !lire {
                match pol.retyp {
                    1 => row[IRE] = pol.reyn,
                    2 => row[IRE] = pol.reyn / acl.sqrt(),
                    3 => row[IRE] = pol.reyn / acl,
                    _ => {}
                }
            }
            if !lima {
                match pol.matyp {
                    1 => row[IMA] = pol.mach,
                    2 => row[IMA] = pol.mach / acl.sqrt(),
                    3 => row[IMA] = pol.mach / acl,
                    _ => {}
                }
            }
            if !linc {
                row[INC] = pol.acrit;
            }
            if !ljtp {
                sd[JTP][1] = pol.xstrip[1];
                sd[JTP][2] = pol.xstrip[2];
            }
            pol.cpol.push(row);
            pol.cpolsd.push(sd);
        }
    }
    false
}

/// POLREF: reads reference polar data file. Returns None on error.
pub fn polref(txt: &str) -> Option<RefPolar> {
    let mut lines: Vec<&str> = txt.lines().collect();
    let mut r = RefPolar::default();
    if lines.is_empty() {
        return None;
    }
    if lines[0].starts_with('#') {
        r.name = substr(lines[0], 2, 80).chars().take(48).collect::<String>();
        lines.remove(0);
    }
    let mut rd = ListReader::new(lines.iter().map(|s| s.to_string()).collect());
    for k in 1..=4 {
        let mut i = 1;
        while i <= NFX {
            match rd.read_reals(2) {
                ListRead::Ok(v) => {
                    r.data[k].push((v[0], v[1]));
                    if v[0] == 999.0 {
                        break;
                    }
                }
                ListRead::Eof => break,
                ListRead::Err => return None,
            }
            i += 1;
        }
        r.nd[k] = i - 1;
        r.data[k].truncate(i - 1);
    }
    Some(r)
}

pub enum ListRead<T> {
    Ok(T),
    Eof,
    Err,
}

/// Minimal Fortran list-directed sequential reader (one READ consumes whole records).
pub struct ListReader {
    lines: Vec<String>,
    pub pos: usize,
}

impl ListReader {
    pub fn new(lines: Vec<String>) -> Self {
        ListReader { lines, pos: 0 }
    }
    fn tokens(line: &str) -> (Vec<String>, bool) {
        // returns tokens and whether a '/' terminator was seen
        let mut toks = Vec::new();
        let mut cur = String::new();
        let mut slash = false;
        let mut last_comma = false;
        for c in line.chars() {
            if c == '/' {
                slash = true;
                break;
            }
            if c == ' ' || c == '\t' || c == ',' {
                if !cur.is_empty() {
                    toks.push(std::mem::take(&mut cur));
                    last_comma = false;
                }
                if c == ',' {
                    if last_comma {
                        toks.push(String::new()); // null value
                    }
                    last_comma = true;
                }
            } else {
                cur.push(c);
            }
        }
        if !cur.is_empty() {
            toks.push(cur);
        }
        (toks, slash)
    }
    /// Read n tokens (as strings) starting at the next record.
    pub fn read_tokens(&mut self, n: usize) -> ListRead<Vec<String>> {
        let mut out = Vec::new();
        if self.pos >= self.lines.len() {
            return ListRead::Eof;
        }
        while out.len() < n {
            if self.pos >= self.lines.len() {
                return ListRead::Eof;
            }
            let (t, slash) = Self::tokens(&self.lines[self.pos]);
            self.pos += 1;
            for tok in t {
                if out.len() < n {
                    out.push(tok);
                }
            }
            if slash {
                break;
            }
        }
        ListRead::Ok(out)
    }
    pub fn read_reals(&mut self, n: usize) -> ListRead<Vec<f64>> {
        match self.read_tokens(n) {
            ListRead::Ok(t) => {
                let mut v = Vec::new();
                for tok in &t {
                    if tok.is_empty() {
                        v.push(f64::NAN);
                        continue;
                    }
                    let tok = match tok.find('*') {
                        Some(k) => &tok[k + 1..],
                        None => tok.as_str(),
                    };
                    match parse_fortran_real(tok) {
                        Some(x) => v.push(x),
                        None => return ListRead::Err,
                    }
                }
                ListRead::Ok(v)
            }
            ListRead::Eof => ListRead::Eof,
            ListRead::Err => ListRead::Err,
        }
    }
}

/// HSORT (Numerical Recipes heapsort): returns 1-based index array.
pub fn hsort(n: usize, a: &[f64]) -> Vec<usize> {
    // a is 1-based
    let mut indx: Vec<usize> = (0..=n).collect();
    if n <= 1 {
        return indx;
    }
    let mut l = n / 2 + 1;
    let mut ir = n;
    loop {
        let indxt;
        let q;
        if l > 1 {
            l -= 1;
            indxt = indx[l];
            q = a[indxt];
        } else {
            indxt = indx[ir];
            q = a[indxt];
            indx[ir] = indx[1];
            ir -= 1;
            if ir == 1 {
                indx[1] = indxt;
                return indx;
            }
        }
        let mut i = l;
        let mut j = l + l;
        while j <= ir {
            if j < ir && a[indx[j]] < a[indx[j + 1]] {
                j += 1;
            }
            if q < a[indx[j]] {
                indx[i] = indx[j];
                i = j;
                j += j;
            } else {
                j = ir + 1;
            }
        }
        indx[i] = indxt;
    }
}

/// Append one Fortran unformatted sequential record (gfortran 4-byte markers).
pub fn ufrec(buf: &mut Vec<u8>, payload: &[u8]) {
    let n = payload.len() as u32;
    buf.extend_from_slice(&n.to_le_bytes());
    buf.extend_from_slice(payload);
    buf.extend_from_slice(&n.to_le_bytes());
}

/// Split a Fortran unformatted sequential file into records. None if malformed.
pub fn ufread(data: &[u8]) -> (Vec<Vec<u8>>, bool) {
    let mut recs = Vec::new();
    let mut p = 0usize;
    while p < data.len() {
        if p + 4 > data.len() {
            return (recs, false);
        }
        let n = u32::from_le_bytes([data[p], data[p + 1], data[p + 2], data[p + 3]]) as usize;
        if p + 8 + n > data.len() {
            return (recs, false);
        }
        let m = u32::from_le_bytes([data[p + 4 + n], data[p + 5 + n], data[p + 6 + n], data[p + 7 + n]]) as usize;
        if m != n {
            return (recs, false);
        }
        recs.push(data[p + 4..p + 4 + n].to_vec());
        p += 8 + n;
    }
    (recs, true)
}

/// Formatted-file record terminators as written by the (Windows) gfortran
/// runtime: "\r\n" on Windows, "\n" elsewhere.
pub fn eol(s: &str) -> String {
    if cfg!(windows) {
        s.replace('\n', "\r\n")
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hs() {
        let a = [0.0, 3.0, 1.0, 2.0, 1.0];
        let ix = hsort(4, &a);
        let v: Vec<f64> = (1..=4).map(|i| a[ix[i]]).collect();
        assert_eq!(v, vec![1.0, 1.0, 2.0, 3.0]);
    }
}
