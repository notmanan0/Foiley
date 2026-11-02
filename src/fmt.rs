//! Fortran-compatible numeric formatting (F, E, G edit descriptors as
//! produced by gfortran), so output files match XFoil byte-for-byte where
//! possible.

/// Fortran Fw.d
pub fn f_fmt(v: f64, w: usize, d: usize) -> String {
    if !v.is_finite() {
        return format!("{:>w$}", if v.is_nan() { "NaN" } else if v > 0.0 { "Infinity" } else { "-Infinity" }, w = w);
    }
    let mut s = format!("{:.*}", d, v);
    if d == 0 {
        s.push('.'); // gfortran prints "500000." for F7.0
    }
    // gfortran prints "-0.000" for small negatives; Rust does too.
    // gfortran drops the leading zero if the field is too narrow
    if s.len() > w {
        if s.starts_with("0.") {
            s = s[1..].to_string();
        } else if s.starts_with("-0.") {
            s = format!("-{}", &s[2..]);
        }
    }
    if s.len() > w {
        return "*".repeat(w);
    }
    format!("{:>w$}", s, w = w)
}

/// Fortran Ew.d (gfortran style: 0.ddddE+xx)
pub fn e_fmt(v: f64, w: usize, d: usize) -> String {
    if !v.is_finite() {
        return format!("{:>w$}", if v.is_nan() { "NaN" } else { "Infinity" }, w = w);
    }
    let (mant, exp) = if v == 0.0 {
        (0.0, 0i32)
    } else {
        let a = v.abs();
        let mut e = a.log10().floor() as i32 + 1;
        let mut m = a / 10f64.powi(e);
        // rounding may push mantissa to 1.0
        let scale = 10f64.powi(d as i32);
        if (m * scale).round() >= scale {
            e += 1;
            m = a / 10f64.powi(e);
        }
        (m, e)
    };
    let digits = format!("{:.*}", d, mant); // "0.dddd"
    let digits = &digits[1..]; // ".dddd"
    let es = if exp.abs() < 100 {
        format!("E{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    } else {
        format!("{}{:03}", if exp < 0 { '-' } else { '+' }, exp.abs())
    };
    let sign = if v < 0.0 { "-" } else { "" };
    let mut s = format!("{}0{}{}", sign, digits, es);
    if s.len() > w {
        s = format!("{}{}{}", sign, digits, es);
    }
    if s.len() > w {
        return "*".repeat(w);
    }
    format!("{:>w$}", s, w = w)
}

/// Fortran Gw.d (gfortran style)
pub fn g_fmt(v: f64, w: usize, d: usize) -> String {
    let a = v.abs();
    if v == 0.0 {
        // gfortran: 0.000000 followed by 4 blanks for exponent field
        let s = format!("{:.*}", d - 1, 0.0);
        return format!("{:>width$}    ", s, width = w - 4);
    }
    if a >= 0.1 - 0.5 * 10f64.powi(-(d as i32)) / 1.0 && a < 10f64.powi(d as i32) - 0.5 {
        // F format with d - k decimals where k = exponent
        let mut k = a.log10().floor() as i32 + 1;
        // account for rounding up
        let dec = (d as i32 - k).max(0) as usize;
        let r = format!("{:.*}", dec, a);
        if r.parse::<f64>().unwrap_or(a) >= 10f64.powi(k) {
            k += 1;
        }
        let dec = (d as i32 - k).max(0) as usize;
        let s = format!("{:.*}", dec, v);
        return format!("{:>width$}    ", s, width = w - 4);
    }
    e_fmt(v, w, d)
}

/// Fortran list-directed REAL output for a single value (approximate gfortran).
pub fn list_real(v: f64) -> String {
    // gfortran default for REAL*8: up to 17 significant digits, use G-like
    if v == 0.0 {
        return "   0.0000000000000000     ".to_string();
    }
    let a = v.abs();
    if (0.1..1e8).contains(&a) {
        let k = a.log10().floor() as i32 + 1;
        let dec = (17 - k).max(0) as usize;
        format!("{:>width$.dec$}", v, width = 24, dec = dec)
    } else {
        e_fmt(v, 25, 17)
    }
}

/// Pad/truncate string to fixed width (CHARACTER*n semantics).
pub fn fixed(s: &str, n: usize) -> String {
    let mut t: String = s.chars().take(n).collect();
    while t.chars().count() < n {
        t.push(' ');
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fmts() {
        assert_eq!(f_fmt(-4.0, 7, 3), " -4.000");
        assert_eq!(f_fmt(0.00770, 10, 5), "   0.00770");
        assert_eq!(e_fmt(1.234e-5, 10, 4), "0.1234E-04");
        assert_eq!(e_fmt(-1.234e-5, 10, 4), "-.1234E-04");
        assert_eq!(e_fmt(1.0e6, 11, 3), "  0.100E+07");
        assert_eq!(g_fmt(1.0, 15, 7), "   1.000000    ");
        assert_eq!(g_fmt(0.5, 15, 7), "  0.5000000    ");
    }
}
