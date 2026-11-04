//! GDES geometry design menu (port of xgdes.f, xtcam.f, the GDES-relevant
//! parts of pntops.f/modify.f and the geometry helpers from xgeom.f).
//!
//! The menu operates on the buffer airfoil (`xb, yb, nb, ...`); EXEC copies
//! it to the current airfoil with ABCOPY.  Plotting is not available:
//! plot-only commands are accepted and ignored, and commands that need
//! cursor input print a note and return.

use crate::cli::{Input, Session};
use crate::fmt::{f_fmt, fixed};
use crate::gdes_geom::*;
use crate::geom::{geopar, getflt, lefind, norm, sopps};
use crate::spline::{curv, deval, scalc, segspl, seval, sinvrt, splina};
use crate::xfoil::{XFoil, DTOR, IBX, IQX, PI};

/// Number of points in thickness/camber arrays (NTX).
const NTX: usize = 2 * IBX;

/// Persistent GDES state (kept across GDES sessions).
#[derive(Clone, Debug)]
pub struct GdesState {
    /// y-symmetry forcing flag (LGSYM)
    pub lgsym: bool,
    /// modified-contour slope matching flag (LGSLOP)
    pub lgslop: bool,
    /// modified-camber slope matching flag (LCSLOP)
    pub lcslop: bool,
    /// thickness/camber plotting flag (LPLCAM)
    pub lplcam: bool,
    /// DCp plot flag in CAMB (LDCPLOT)
    pub ldcplot: bool,
    /// plot toggles (no effect without graphics)
    pub lgtick: bool,
    pub lggrid: bool,
    pub lgparm: bool,
    /// tick-mark size (GTICK) and plot-object size (SIZE)
    pub gtick: f64,
    pub size: f64,
    /// default overlay filename (ONAME)
    pub oname: String,
    // camber / thickness arrays (XCM, YCM, ... in /CR14/)
    pub xcm: Vec<f64>,
    pub ycm: Vec<f64>,
    pub ycmp: Vec<f64>,
    pub xtk: Vec<f64>,
    pub ytk: Vec<f64>,
    pub ytkp: Vec<f64>,
    pub ncm: usize,
    pub ntk: usize,
    // CAMB added camber / loading arrays (/CR18/)
    pub xcadd: Vec<f64>,
    pub ycadd: Vec<f64>,
    pub ycaddp: Vec<f64>,
    pub xpadd: Vec<f64>,
    pub ypadd: Vec<f64>,
    pub ypaddp: Vec<f64>,
    pub xcam: Vec<f64>,
    pub ycam: Vec<f64>,
    pub ycamp: Vec<f64>,
    pub pcam: Vec<f64>,
    pub pcamp: Vec<f64>,
    pub ncam: usize,
    pub ncadd: usize,
    pub npadd: usize,
}

impl Default for GdesState {
    fn default() -> Self {
        let v = |n: usize| vec![0.0; n + 2];
        GdesState {
            lgsym: false,
            lgslop: true,
            lcslop: true,
            lplcam: false,
            ldcplot: false,
            lgtick: true,
            lggrid: true,
            lgparm: true,
            gtick: 0.0005,
            size: 10.0,
            oname: String::new(),
            xcm: v(2 * IBX),
            ycm: v(2 * IBX),
            ycmp: v(2 * IBX),
            xtk: v(2 * IBX),
            ytk: v(2 * IBX),
            ytkp: v(2 * IBX),
            ncm: 0,
            ntk: 0,
            xcadd: v(NTX),
            ycadd: v(NTX),
            ycaddp: v(NTX),
            xpadd: v(NTX),
            ypadd: v(NTX),
            ypaddp: v(NTX),
            xcam: v(NTX),
            ycam: v(NTX),
            ycamp: v(NTX),
            pcam: v(NTX),
            pcamp: v(NTX),
            ncam: 0,
            ncadd: 0,
            npadd: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

/// Parsed command arguments: RINPUT(1..20), IINPUT(1..20), NINPUT.
struct Args {
    r: Vec<f64>,
    i: Vec<i64>,
    n: usize,
}

fn parse_args(comarg: &str) -> Args {
    // ignore everything after a "!" (GETINT/GETFLT)
    let a = match comarg.find('!') {
        Some(k) => &comarg[..k],
        None => comarg,
    };
    let mut r = vec![0.0; 21];
    let mut i = vec![0i64; 21];
    let vals = getflt(a, 20).unwrap_or_default();
    for (k, v) in vals.iter().enumerate() {
        r[k + 1] = *v;
    }
    // GETINT: list-directed integer read; any non-integer token is an error
    let toks: Vec<&str> = a.split([' ', ',', '\t']).filter(|t| !t.is_empty()).take(20).collect();
    // (items before the first bad one are assigned, as in a gfortran list-directed READ)
    for (k, t) in toks.iter().enumerate() {
        match t.parse::<i64>() {
            Ok(v) => i[k + 1] = v,
            Err(_) => break,
        }
    }
    Args { r, i, n: vals.len() }
}

fn cursor_note(com: &str) {
    crate::outln!(" {}: cursor input is not available (no graphics).  No action taken.", com.trim());
}

/// READR: read up to n reals into var, leaving unchanged on <return>.  Returns error flag.
fn readr(input: &mut Input, var: &mut [f64]) -> bool {
    let line = input.ask("");
    if input.eof {
        return false;
    }
    match getflt(&line, 40) {
        Some(v) => {
            for (k, x) in v.iter().enumerate() {
                if k < var.len() {
                    var[k] = *x;
                }
            }
            false
        }
        None => true,
    }
}

/// READI: read up to n integers, leaving unchanged on <return>.  Returns error flag.
fn readi(input: &mut Input, var: &mut [i64]) -> bool {
    let line = input.ask("");
    if input.eof {
        return false;
    }
    let toks: Vec<&str> = line.split([' ', ',', '\t']).filter(|t| !t.is_empty()).collect();
    let mut tmp = Vec::new();
    for t in toks {
        match t.parse::<i64>() {
            Ok(v) => tmp.push(v),
            Err(_) => return true,
        }
    }
    for (k, x) in tmp.iter().enumerate() {
        if k < var.len() {
            var[k] = *x;
        }
    }
    false
}

/// List-directed read of two reals (may span lines).  None on EOF.
fn read2(input: &mut Input, prompt: &str) -> Option<(f64, f64)> {
    'outer: loop {
        let mut v: Vec<f64> = Vec::new();
        let mut first = true;
        while v.len() < 2 {
            let line = input.ask(if first { prompt } else { "" });
            first = false;
            if input.eof {
                return None;
            }
            match getflt(&line, 2) {
                Some(a) => v.extend(a),
                None => continue 'outer,
            }
        }
        return Some((v[0], v[1]));
    }
}

/// Spline the buffer airfoil (SCALC + SEGSPL x,y).
fn bspline(xf: &mut XFoil) {
    let nb = xf.nb;
    scalc(&xf.xb, &xf.yb, &mut xf.sb, nb);
    segspl(&xf.xb, &mut xf.xbp, &xf.sb, nb);
    segspl(&xf.yb, &mut xf.ybp, &xf.sb, nb);
}

/// GEOPAR on the buffer airfoil, setting the /CR14/ buffer parameters.
fn bgeopar(xf: &mut XFoil) {
    let g = geopar(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, xf.nb, true);
    xf.sble = g.sle;
    xf.chordb = g.chord;
    xf.areab = g.area;
    xf.radble = g.radle;
    xf.angbte = g.angte;
    xf.ei11ba = g.ei11a;
    xf.ei22ba = g.ei22a;
    xf.apx1ba = g.apx1a;
    xf.apx2ba = g.apx2a;
    xf.ei11bt = g.ei11t;
    xf.ei22bt = g.ei22t;
    xf.apx1bt = g.apx1t;
    xf.apx2bt = g.apx2t;
    xf.thickb = g.thick;
    xf.cambrb = g.cambr;
}

/// Respline buffer and recompute its geometric parameters.
fn bupdate(xf: &mut XFoil) {
    bspline(xf);
    bgeopar(xf);
}

/// GETCAM on the buffer airfoil into the persistent camber/thickness arrays.
fn getcam_b(st: &mut GdesState, xf: &XFoil) {
    let GdesState { xcm, ycm, ncm, xtk, ytk, ntk, .. } = st;
    getcam(xcm, ycm, ncm, xtk, ytk, ntk, &xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, xf.nb);
}

/// GETMAX for camber and thickness.  Returns (cxmax, cymax, txmax, tymax).
fn getmax_ct(st: &mut GdesState) -> (f64, f64, f64, f64) {
    let (cx, cy) = getmax(&st.xcm, &st.ycm, &mut st.ycmp, st.ncm);
    let (tx, ty) = getmax(&st.xtk, &st.ytk, &mut st.ytkp, st.ntk);
    (cx, cy, tx, ty)
}

fn print_tc(tymax: f64, txmax: f64, cymax: f64, cxmax: f64) {
    crate::outln!(" Max thickness = {}  at x = {}", f_fmt(2.0 * tymax, 8, 4), f_fmt(txmax, 7, 3));
    crate::outln!(" Max camber    = {}  at x = {}", f_fmt(cymax, 8, 4), f_fmt(cxmax, 7, 3));
}

/// GOFINI: only its non-graphical side effect (GETCAM) is kept.
fn gofini(s: &mut Session) {
    getcam_b(&mut s.gdes_state, &s.xf);
}

/// PLTCAM: only its non-graphical side effect (GETCAM) is kept.
fn pltcam(s: &mut Session) {
    if s.gdes_state.lplcam {
        getcam_b(&mut s.gdes_state, &s.xf);
    }
}

/// PLOTG: only its non-graphical side effects are kept.
fn plotg(s: &mut Session) {
    pltcam(s);
}

fn enable_tcplot(s: &mut Session) {
    if !s.gdes_state.lplcam {
        crate::outln!(" Enabling camber,thickness plotting");
        s.gdes_state.lplcam = true;
        gofini(s);
        plotg(s);
    }
}

// ---------------------------------------------------------------------------
// GDES menu
// ---------------------------------------------------------------------------

const GDES_HELP: &str = "
   <cr>     Return to Top Level
   !        Redo previous command

   GSET     Set buffer  airfoil <== current airfoil
   eXec     Set current airfoil <== buffer  airfoil
   SYMM     Toggle y-symmetry flag

   ADEG r   Rotate about origin (degrees)
   ARAD r   Rotate about origin (radians)
   Tran rr  Translate
   Scal r   Scale about origin
   LINS rr. Linearly-varying y scale
   DERO     Derotate (set chord line level)

   TGAP rr  Change trailing edge gap
   LERA rr  Change leading edge radius

   TCPL     Toggle thickness and camber plotting
   TFAC rr  Scale existing thickness and camber
   TSET rr  Set new thickness and camber
   HIGH rr  Move camber and thickness highpoints
  .CAMB     Modify camber shape directly or via loading

   BEND     Display structural properties of buffer airfoil

   Flap rrr Deflect trailing edge flap

   Modi     Modify contour via cursor
   SLOP     Toggle modified-contour slope matching flag

   CORN     Double point with cursor (set sharp corner)
   ADDP     Add    point with cursor or keyboard x,y
   MOVP     Move   point with cursor or keyboard x,y
   DELP     Delete point with cursor

   UNIT     Normalize buffer airfoil to unit chord
   Dist     Determine distance between 2 cursor points
   CLIS     List curvatures
   CPLO     Plot curvatures
   CANG     List panel corner angles
   CADD ri. Add points at corners exceeding angle threshold

   Plot     Replot buffer airfoil
   INPL     Replot buffer airfoil without scaling (in inches)
   Blow     Blowup plot region
   Rese     Reset plot scale and origin
   Wind     Plot window adjust via cursor and keys

   TSIZ r   Change tick-mark size
   TICK     Toggle node tick-mark plotting
   GRID     Toggle grid plotting
   GPAR     Toggle geometric parameter plotting
   Over f   Overlay disk file airfoil

   SIZE r   Change absolute plot-object size
  .ANNO     Annotate plot
   HARD     Hardcopy current plot

   NAME s   Specify new airfoil name
   NINC     Increment name version number";

/// GDES menu entry point.
pub fn gdes(s: &mut Session) {
    if s.xf.nb == 0 {
        crate::outln!();
        crate::outln!(" ***  No airfoil available  ***");
        return;
    }
    s.gdes_state.lplcam = false;
    crate::outln!();
    crate::outln!(" You are working with the buffer airfoil");
    gofini(s);
    plotg(s);

    let mut comand = String::from("****");
    let mut comarg = String::new();
    loop {
        // label 500
        let comold = comand.clone();
        let argold = comarg.clone();
        // label 501
        loop {
            let prompt = if s.gdes_state.lgsym { ".GDESs" } else { ".GDES" };
            let (c, a) = s.input.askc(prompt);
            if s.input.eof {
                return;
            }
            comand = c;
            comarg = a;
            if comand.starts_with('!') {
                if comold == "****" {
                    crate::outln!(" Previous .GDES command not valid");
                    continue;
                }
                comand = comold.clone();
                comarg = argold.clone();
            }
            break;
        }

        if comand.is_empty() {
            if !s.xf.lgsame {
                crate::outln!();
                crate::outln!(" Buffer airfoil is not identical to current airfoil");
            }
            return;
        }

        let a = parse_args(&comarg);
        let (rinput, ninput) = (&a.r, a.n);

        match comand.as_str() {
            "?" => crate::outln!("{}", GDES_HELP),
            "Z" | "U" => {}
            "GSET" => {
                let xf = &mut s.xf;
                xf.nb = xf.n;
                for i in 1..=xf.nb {
                    xf.xb[i] = xf.x[i];
                    xf.yb[i] = xf.y[i];
                }
                xf.lgsame = true;
                bupdate(xf);
                plotg(s);
                if s.gdes_state.lgsym {
                    zercam(s);
                }
            }
            "EXEC" | "X" => s.xf.abcopy(true),
            "SYMM" => {
                s.gdes_state.lgsym = !s.gdes_state.lgsym;
                if s.gdes_state.lgsym {
                    crate::outln!(" y-symmetry forcing enabled.");
                    zercam(s);
                } else {
                    crate::outln!(" y-symmetry forcing disabled.");
                }
            }
            "ADEG" | "ARAD" => {
                let arad = if comand == "ADEG" {
                    let adeg = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter angle change (deg)^", 0.0) };
                    s.xf.adeg = adeg;
                    adeg * PI / 180.0
                } else if ninput >= 1 {
                    rinput[1]
                } else {
                    s.input.askr("Enter angle change (rad)^", 0.0)
                };
                let xf = &mut s.xf;
                let nb = xf.nb;
                rotate(&mut xf.xb, &mut xf.yb, nb, arad);
                segspl(&xf.xb, &mut xf.xbp, &xf.sb, nb);
                segspl(&xf.yb, &mut xf.ybp, &xf.sb, nb);
                xf.apx1ba -= arad;
                xf.apx2ba -= arad;
                xf.apx1bt -= arad;
                xf.apx2bt -= arad;
            }
            "TRAN" | "T" => {
                let (delx, dely);
                if ninput >= 2 {
                    delx = rinput[1];
                    dely = rinput[2];
                } else if ninput >= 1 {
                    delx = rinput[1];
                    dely = s.input.askr("Enter delta(y)^", 0.0);
                } else {
                    delx = s.input.askr("Enter delta(x)^", 0.0);
                    dely = s.input.askr("Enter delta(y)^", 0.0);
                }
                let xf = &mut s.xf;
                for i in 1..=xf.nb {
                    xf.xb[i] += delx;
                    xf.yb[i] += dely;
                }
                xf.lgsame = false;
            }
            "SCAL" | "S" => {
                let fac = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter scale factor (0 for separate x,y scales)^", 1.0) };
                let mut xxfac = fac;
                let mut yyfac = fac;
                if fac == 0.0 {
                    if ninput >= 3 {
                        xxfac = rinput[2];
                        yyfac = rinput[3];
                    } else {
                        xxfac = s.input.askr("Enter x scale factor^", 1.0);
                        yyfac = s.input.askr("Enter y scale factor^", 1.0);
                    }
                }
                let xf = &mut s.xf;
                let nb = xf.nb;
                for i in 1..=nb {
                    xf.xb[i] *= xxfac;
                    xf.yb[i] *= yyfac;
                }
                if xxfac * yyfac < 0.0 {
                    for i in 1..=nb / 2 {
                        xf.xb.swap(i, nb - i + 1);
                        xf.yb.swap(i, nb - i + 1);
                    }
                }
                bupdate(xf);
                xf.lgsame = false;
            }
            "LINS" => {
                let mut ninp = ninput;
                let (mut xoc1, mut yfac1, mut xoc2, mut yfac2);
                loop {
                    if ninp >= 4 {
                        xoc1 = rinput[1];
                        yfac1 = rinput[2];
                        xoc2 = rinput[3];
                        yfac2 = rinput[4];
                    } else {
                        match read2(&mut s.input, "Location 1...  enter  x/c, y-scale :  ") {
                            Some((a, b)) => {
                                xoc1 = a;
                                yfac1 = b;
                            }
                            None => return,
                        }
                        match read2(&mut s.input, "Location 2...  enter  x/c, y-scale :  ") {
                            Some((a, b)) => {
                                xoc2 = a;
                                yfac2 = b;
                            }
                            None => return,
                        }
                    }
                    if (xoc1 - xoc2).abs() < 1.0e-5 {
                        crate::outln!(" x/c locations 1 and 2 must be different");
                        ninp = 0;
                        continue;
                    }
                    break;
                }
                let xf = &mut s.xf;
                let nb = xf.nb;
                xf.sble = lefind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
                let xle = seval(xf.sble, &xf.xb, &xf.xbp, &xf.sb, nb);
                let _yle = seval(xf.sble, &xf.yb, &xf.ybp, &xf.sb, nb);
                let xte = 0.5 * (xf.xb[1] + xf.xb[nb]);
                for i in 1..=nb {
                    let xoc = (xf.xb[i] - xle) / (xte - xle);
                    let fr1 = (xoc2 - xoc) / (xoc2 - xoc1);
                    let fr2 = (xoc - xoc1) / (xoc2 - xoc1);
                    let yyfac = fr1 * yfac1 + fr2 * yfac2;
                    xf.yb[i] *= yyfac;
                }
                bupdate(xf);
                xf.lgsame = false;
            }
            "DERO" => {
                let xf = &mut s.xf;
                let nb = xf.nb;
                xf.sble = lefind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
                let xle = seval(xf.sble, &xf.xb, &xf.xbp, &xf.sb, nb);
                let yle = seval(xf.sble, &xf.yb, &xf.ybp, &xf.sb, nb);
                let xte = 0.5 * (xf.xb[1] + xf.xb[nb]);
                let yte = 0.5 * (xf.yb[1] + xf.yb[nb]);
                let arad = (yte - yle).atan2(xte - xle);
                rotate(&mut xf.xb, &mut xf.yb, nb, arad);
                crate::outln!("\nRotating buffer airfoil by {} deg.", f_fmt(arad / DTOR, 8, 3));
                bupdate(xf);
                xf.lgsame = false;
            }
            "TGAP" => tgap(s, rinput, ninput),
            "LERA" => lerad(s, rinput, ninput),
            "TCPL" => {
                s.gdes_state.lplcam = !s.gdes_state.lplcam;
                gofini(s);
                plotg(s);
            }
            "TFAC" | "TSET" | "HIGH" => {
                enable_tcplot(s);
                match comand.as_str() {
                    "TFAC" => tcscal(s, rinput, ninput),
                    "TSET" => tcset(s, rinput, ninput),
                    _ => hipnt(s, rinput, ninput),
                }
                pltcam(s);
                s.xf.lgsame = false;
            }
            "CAMB" => {
                if s.gdes_state.lgsym {
                    crate::outln!(" Disabling symmetry enforcement.");
                    s.gdes_state.lgsym = false;
                }
                camb(s);
                if s.input.eof {
                    return;
                }
            }
            "BEND" => bendump(s.xf.nb, &s.xf.xb, &s.xf.yb),
            "CANG" => {
                cang_p(&s.xf.xb, &s.xf.yb, s.xf.nb, 2);
            }
            "CADD" => cadd(s, &a),
            "CLIS" => {
                let xf = &s.xf;
                clis(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, xf.nb);
            }
            "CPLO" => pltcrv(s),
            "FLAP" | "F" => flap(s, rinput, ninput),
            "MODI" | "M" => cursor_note(&comand),
            "SLOP" => {
                s.gdes_state.lgslop = !s.gdes_state.lgslop;
                if s.gdes_state.lgslop {
                    crate::outln!(" Modified segment will be made tangent at endpoints");
                } else {
                    crate::outln!(" Modified segment will not be made tangent at endpoints");
                }
            }
            "CORN" => {
                if s.xf.nb == 2 * IQX {
                    crate::outln!(" Buffer airfoil arrays will overflow.  No action taken.");
                } else {
                    cursor_note(&comand);
                }
            }
            "ADDP" => {
                if s.xf.nb >= IBX {
                    crate::outln!();
                    crate::outln!(" Buffer airfoil arrays will overflow.  No action taken.");
                } else {
                    cursor_note(&comand);
                }
            }
            "DELP" | "MOVP" => cursor_note(&comand),
            "UNIT" => {
                let xf = &mut s.xf;
                let nb = xf.nb;
                norm(&mut xf.xb, &mut xf.xbp, &mut xf.yb, &mut xf.ybp, &mut xf.sb, nb);
                xf.lgsame = false;
                bupdate(xf);
            }
            "DIST" | "D" => cursor_note(&comand),
            "HARD" => {}
            "PLOT" | "P" | "INPL" => plotg(s),
            "BLOW" | "B" | "WIND" | "W" => cursor_note(&comand),
            "RESE" | "R" => {
                gofini(s);
                plotg(s);
            }
            "TSIZ" => {
                if ninput >= 1 {
                    s.gdes_state.gtick = rinput[1];
                } else {
                    crate::outln!(" Current tick-mark size (as fraction of perimeter) ={}", lreal(s.gdes_state.gtick));
                    s.gdes_state.gtick = s.input.askr("Enter new tick-mark size^", s.gdes_state.gtick);
                }
                plotg(s);
            }
            "TICK" => {
                s.gdes_state.lgtick = !s.gdes_state.lgtick;
                plotg(s);
            }
            "GRID" => {
                s.gdes_state.lggrid = !s.gdes_state.lggrid;
                plotg(s);
            }
            "GPAR" => {
                s.gdes_state.lgparm = !s.gdes_state.lgparm;
                plotg(s);
            }
            "SIZE" => {
                if ninput >= 1 {
                    s.gdes_state.size = rinput[1];
                } else {
                    crate::outln!(" Current plot-object size ={}", lreal(s.gdes_state.size));
                    s.gdes_state.size = s.input.askr("Enter new plot-object size^", s.gdes_state.size);
                }
                plotg(s);
            }
            "OVER" | "O" => over(s, &comarg),
            "ANNO" => cursor_note(&comand),
            "NAME" => {
                if comarg.trim().is_empty() {
                    nammod(s, 0, -1);
                } else {
                    s.xf.name = comarg.clone();
                }
                s.xf.name = s.xf.name.trim().to_string();
            }
            "NINC" => {
                nammod(s, 1, 1);
                s.xf.name = s.xf.name.trim().to_string();
            }
            "NDEC" => {
                nammod(s, -1, 1);
                s.xf.name = s.xf.name.trim().to_string();
            }
            "SINT" => {
                let xf = &mut s.xf;
                let nb = xf.nb;
                splnxy(&xf.xb, &mut xf.xbp, &xf.yb, &mut xf.ybp, &mut xf.sb, nb);
            }
            _ => {
                crate::outln!(" Command {} not recognized.  Type a \" ? \" for list.", fixed(&comand, 4));
                comand = "****".to_string();
            }
        }
        if s.input.eof {
            return;
        }
    }
}

// ---------------------------------------------------------------------------
// PLTCRV (xgeom.f): curvature plot submenu, non-graphical part
// ---------------------------------------------------------------------------

fn pltcrv(s: &mut Session) {
    let mut lcvexp = false;
    let mut cvexp: f64 = 1.0 / 3.0;
    'recalc: loop {
        let (cvle, cvmax, svmax, xvmax, cvmin, svmin, xvmin, sble, xle);
        {
            let xf = &s.xf;
            let nb = xf.nb;
            let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
            sble = xf.sble;
            let sbtot = 0.5 * (sb[nb] - sb[1]);
            xle = seval(sble, xb, xbp, sb, nb);
            let mut c = curv(sble, xb, xbp, yb, ybp, sb, nb) * sbtot;
            if lcvexp {
                c = c.powf(cvexp);
            }
            cvle = c;
            let (mut cmx, mut smx, mut xmx) = (cvle, sble, xle);
            let (mut cmn, mut smn, mut xmn) = (cvle, sble, xle);
            for i in 1..=nb {
                let mut cv = curv(sb[i], xb, xbp, yb, ybp, sb, nb) * sbtot;
                if lcvexp {
                    if cv > 0.0 {
                        cv = cv.powf(cvexp);
                    } else if cv == 0.0 {
                        cv = 0.0;
                    } else if cv < 0.0 {
                        let cvsgn = 1.0f64.copysign(cv);
                        cv = cvsgn * cv.abs().powf(cvexp);
                    }
                }
                if cv > cmx {
                    cmx = cv;
                    smx = sb[i];
                    xmx = xb[i];
                }
                if cv < cmn {
                    cmn = cv;
                    smn = sb[i];
                    xmn = xb[i];
                }
            }
            cvmax = cmx;
            svmax = smx;
            xvmax = xmx;
            cvmin = cmn;
            svmin = smn;
            xvmin = xmn;
        }
        loop {
            crate::outln!("  ");
            crate::outln!(" Airfoil curvature (yellow-upper, cyan-lower) ");
            if lcvexp {
                crate::outln!("  Range compressed using CV=(curvature)^n with n ={}", lreal(cvexp));
            }
            crate::outln!("  CVLE  = {}  at S = {}  at X = {}", lreal(cvle), lreal(sble), lreal(xle));
            crate::outln!("  CVmax = {}  at S = {}  at X = {}", lreal(cvmax), lreal(svmax), lreal(xvmax));
            crate::outln!("  CVmin = {}  at S = {}  at X = {}", lreal(cvmin), lreal(svmin), lreal(xvmin));
            crate::outln!("  ");
            crate::outln!(" Enter C for curvature    plot");
            crate::outln!(" Enter N for curvature**N plot");
            crate::outln!(" Hit <cr> to exit");
            let (c, a) = s.input.askc("..CPLO");
            if s.input.eof {
                return;
            }
            let ans = match c.chars().next() {
                Some(ch) => ch.to_ascii_uppercase(),
                None => return,
            };
            let r = getflt(&a, 1).unwrap_or_default();
            if ans == 'N' {
                cvexp = if !r.is_empty() { r[0] } else { s.input.askr("Enter curvature exponent (default 0.3)^", 0.3) };
                lcvexp = true;
                continue 'recalc;
            }
            if ans == 'C' {
                lcvexp = false;
                continue 'recalc;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// NAMMOD (xoper.f)
// ---------------------------------------------------------------------------

fn nammod(s: &mut Session, kdel: i64, kmod0: i64) {
    let mut name = s.xf.name.trim().to_string();
    let kb1 = name.find('[').map(|k| k + 1).unwrap_or(0);
    let kb2 = name.find(']').map(|k| k + 1).unwrap_or(0);
    let mut namdef = name.clone();
    let mut kmodp: i64 = 0;
    let mut read_err = false;
    if kb1 != 0 && kb2 != 0 && kb2 as i64 - kb1 as i64 > 1 {
        let inner = &name[kb1..kb2 - 1];
        match getflt(inner, 1).and_then(|v| v.first().copied()).filter(|v| v.fract() == 0.0) {
            Some(v) => {
                let kmod = (v as i64).abs();
                kmodp = (kmod + kdel) % 100;
                name = if kb1 >= 2 { name[..kb1 - 1].trim().to_string() } else { String::new() };
            }
            None => read_err = true,
        }
    } else if kmod0 > 0 {
        kmodp = kmod0 % 100;
    }
    if !read_err {
        if kmodp >= 10 {
            namdef = format!("{} [{:2}]", name, kmodp);
        } else if kmodp >= 1 {
            namdef = format!("{} [{:1}]", name, kmodp);
        }
    }
    crate::outln!();
    let ans = s.input.ask(&format!("Enter airfoil name or <return> for default:  {}", namdef));
    s.xf.name = if ans.trim().is_empty() { namdef } else { ans };
}

// ---------------------------------------------------------------------------
// OVER: overlay airfoil (only the non-graphical part)
// ---------------------------------------------------------------------------

fn over(s: &mut Session, fname1: &str) {
    let fname = if !fname1.trim().is_empty() {
        fname1.trim().to_string()
    } else if !s.gdes_state.oname.is_empty() {
        let f = s.input.ask(&format!("Enter filename:  {}", s.gdes_state.oname));
        if f.trim().is_empty() {
            s.gdes_state.oname.clone()
        } else {
            f.trim().to_string()
        }
    } else {
        s.input.ask("Enter filename^").trim().to_string()
    };
    let r = crate::geom::aread(&fname, 2 * IQX);
    if r.itype == 0 {
        return;
    }
    s.gdes_state.oname = fname;
    let nn = r.n;
    let mut w1 = r.x.clone();
    let mut w2 = r.y.clone();
    w1.resize(nn + 2, 0.0);
    w2.resize(nn + 2, 0.0);
    let mut w3 = vec![0.0; nn + 2];
    let mut w4 = vec![0.0; nn + 2];
    let mut w5 = vec![0.0; nn + 2];
    if s.xf.lnorm {
        norm(&mut w1, &mut w3, &mut w2, &mut w4, &mut w5, nn);
    } else {
        scalc(&w1, &w2, &mut w5, nn);
        segspl(&w1, &mut w3, &w5, nn);
        segspl(&w2, &mut w4, &w5, nn);
    }
    // GEOPAR of overlay airfoil (buffer parameters are restored afterwards)
    let _ = geopar(&w1, &w3, &w2, &w4, &w5, nn, true);
}

// ---------------------------------------------------------------------------
// TGAP, LERAD (xgdes.f)
// ---------------------------------------------------------------------------

// MIN(MAX(..)) as in the Fortran; f64::clamp differs for NaN
#[allow(clippy::manual_clamp)]
fn tgap(s: &mut Session, rinput: &[f64], ninput: usize) {
    let (sble, xble, yble, xbte, ybte, chbsq, gap, dxu, dyu);
    {
        let xf = &mut s.xf;
        let nb = xf.nb;
        xf.sble = lefind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
        sble = xf.sble;
        xble = seval(sble, &xf.xb, &xf.xbp, &xf.sb, nb);
        yble = seval(sble, &xf.yb, &xf.ybp, &xf.sb, nb);
        xbte = 0.5 * (xf.xb[1] + xf.xb[nb]);
        ybte = 0.5 * (xf.yb[1] + xf.yb[nb]);
        chbsq = (xbte - xble).powi(2) + (ybte - yble).powi(2);
        let dxn = xf.xb[1] - xf.xb[nb];
        let dyn_ = xf.yb[1] - xf.yb[nb];
        gap = (dxn * dxn + dyn_ * dyn_).sqrt();
        if gap > 0.0 {
            dxu = dxn / gap;
            dyu = dyn_ / gap;
        } else {
            dxu = -0.5 * (xf.ybp[nb] - xf.ybp[1]);
            dyu = 0.5 * (xf.xbp[nb] - xf.xbp[1]);
        }
    }
    let (gapnew, mut doc);
    if ninput >= 2 {
        gapnew = rinput[1];
        doc = rinput[2];
    } else if ninput >= 1 {
        gapnew = rinput[1];
        doc = s.input.askr("Enter blending distance/c (0..1)^", 1.0);
    } else {
        crate::outln!("\n Current gap ={}", f_fmt(gap, 9, 5));
        gapnew = s.input.askr("Enter new gap^", 0.0);
        doc = s.input.askr("Enter blending distance/c (0..1)^", 1.0);
    }
    doc = doc.max(0.0).min(1.0);
    let dgap = gapnew - gap;
    let xf = &mut s.xf;
    let nb = xf.nb;
    for i in 1..=nb {
        let xoc = ((xf.xb[i] - xble) * (xbte - xble) + (xf.yb[i] - yble) * (ybte - yble)) / chbsq;
        let tfac = if doc == 0.0 {
            if i == 1 || i == nb { 1.0 } else { 0.0 }
        } else {
            let arg = ((1.0 - xoc) * (1.0 / doc - 1.0)).min(15.0);
            (-arg).exp()
        };
        if xf.sb[i] <= sble {
            xf.xb[i] += 0.5 * dgap * xoc * tfac * dxu;
            xf.yb[i] += 0.5 * dgap * xoc * tfac * dyu;
        } else {
            xf.xb[i] -= 0.5 * dgap * xoc * tfac * dxu;
            xf.yb[i] -= 0.5 * dgap * xoc * tfac * dyu;
        }
    }
    xf.lgsame = false;
    bupdate(xf);
}

fn lerad(s: &mut Session, rinput: &[f64], ninput: usize) {
    let (rfac, mut doc);
    if ninput >= 2 {
        rfac = rinput[1];
        doc = rinput[2];
    } else if ninput >= 1 {
        rfac = rinput[1];
        doc = s.input.askr("Enter blending distance/c from LE^", 1.0);
    } else {
        rfac = s.input.askr("Enter approx. new/old LE radius scaling ratio^", 1.0);
        doc = s.input.askr("Enter blending distance/c from LE^", 1.0);
    }
    doc = doc.max(0.001);
    let xf = &mut s.xf;
    let nb = xf.nb;
    let mut w1 = vec![0.0; nb + 2];
    let mut w2 = vec![0.0; nb + 2];
    lerscl(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb, doc, rfac, &mut w1, &mut w2);
    for i in 1..=nb {
        xf.xb[i] = w1[i];
        xf.yb[i] = w2[i];
    }
    xf.lgsame = false;
    bupdate(xf);
    let mut cvmax: f64 = 0.0;
    for i in nb / 4..=(3 * nb) / 4 {
        let cv = curv(xf.sb[i], &xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
        cvmax = cv.abs().max(cvmax);
    }
    let radius = 1.0 / cvmax;
    crate::outln!("\n New LE radius = {}", f_fmt(radius, 7, 5));
}

// ---------------------------------------------------------------------------
// CADD
// ---------------------------------------------------------------------------

fn cadd(s: &mut Session, a: &Args) {
    let (_imax, amax) = cang_p(&s.xf.xb, &s.xf.yb, s.xf.nb, 2);
    crate::outln!();
    let nb = s.xf.nb;
    let mut xbmin = s.xf.xb[1];
    let mut xbmax = s.xf.xb[1];
    for i in 1..=nb {
        xbmin = xbmin.min(s.xf.xb[i]);
        xbmax = xbmax.max(s.xf.xb[i]);
    }
    let mut atol = 0.5 * amax;
    let mut ispl: i64 = 1;
    let mut xrf = [xbmin - 0.1 * (xbmax - xbmin), xbmax + 0.1 * (xbmax - xbmin)];
    let ninput = a.n;
    // stage: 70 = ask all, 71 = ask ispl.., 72 = ask limits, 74 = go
    let mut stage = if ninput == 0 {
        70
    } else if ninput <= 1 {
        atol = a.r[1];
        71
    } else if ninput <= 2 {
        atol = a.r[1];
        ispl = a.i[2];
        72
    } else if ninput <= 4 {
        atol = a.r[1];
        ispl = a.i[2];
        xrf = [a.r[3], a.r[4]];
        74
    } else {
        70
    };
    loop {
        match stage {
            70 => {
                crate::outln!(" Enter corner angle criterion for refinement (deg):{}", f_fmt(atol, 8, 3));
                let mut v = [atol];
                if readr(&mut s.input, &mut v) {
                    continue;
                }
                atol = v[0];
                stage = 71;
            }
            71 => {
                crate::outln!(" Enter type of spline parameter (1=uniform, 2=arclength):{}", i_fmt(ispl, 4));
                let mut v = [ispl];
                if readi(&mut s.input, &mut v) {
                    continue;
                }
                ispl = v[0];
                if ispl <= 0 {
                    return;
                }
                if ispl > 2 {
                    continue;
                }
                stage = 72;
            }
            72 => {
                crate::outln!(" Enter refinement x limits:{}{}", f_fmt(xrf[0], 10, 5), f_fmt(xrf[1], 10, 5));
                let mut v = xrf;
                if readr(&mut s.input, &mut v) {
                    continue;
                }
                xrf = v;
                stage = 74;
            }
            _ => break,
        }
        if s.input.eof {
            return;
        }
    }

    let xf = &mut s.xf;
    if ispl == 1 {
        xf.sb[1] = 0.0;
        for i in 2..=nb {
            if xf.xb[i] == xf.xb[i - 1] && xf.yb[i] == xf.yb[i - 1] {
                xf.sb[i] = xf.sb[i - 1];
            } else {
                xf.sb[i] = xf.sb[i - 1] + 1.0;
            }
        }
        segspl(&xf.xb, &mut xf.xbp, &xf.sb, nb);
        segspl(&xf.yb, &mut xf.ybp, &xf.sb, nb);
    }
    let mut w1 = vec![0.0; 6 * IQX + 2];
    let mut w2 = vec![0.0; 6 * IQX + 2];
    let nnew = arefine(&xf.xb, &xf.yb, &xf.sb, &xf.xbp, &xf.ybp, nb, atol, IBX, &mut w1, &mut w2, xrf[0], xrf[1]);
    let nbadd = nnew as i64 - nb as i64;
    crate::outln!(" Number of points added: {}", lint(nbadd));
    xf.nb = nnew;
    for i in 1..=nnew {
        xf.xb[i] = w1[i];
        xf.yb[i] = w2[i];
    }
    xf.lgsame = false;
    if nnew < 2 {
        return;
    }
    bupdate(xf);
    cang_p(&xf.xb, &xf.yb, xf.nb, 1);
}

// ---------------------------------------------------------------------------
// FLAP (xgdes.f)
// ---------------------------------------------------------------------------

/// GETXYF: get flap hinge location.  Returns (tops, bots).
fn getxyf(s: &mut Session, xfh: &mut f64, yfh: &mut f64) -> (f64, f64) {
    if *xfh == -999.0 {
        *xfh = s.input.askr("Enter flap hinge x location^", *xfh);
    }
    let xf = &s.xf;
    let n = xf.nb;
    let (x, xp, y, yp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
    let mut tops = sb[1] + (x[1] - *xfh);
    let mut bots = sb[n] - (x[n] - *xfh);
    sinvrt(&mut tops, *xfh, x, xp, sb, n);
    sinvrt(&mut bots, *xfh, x, xp, sb, n);
    let topy = seval(tops, y, yp, sb, n);
    let boty = seval(bots, y, yp, sb, n);
    crate::outln!("\n  Top    surface:  y ={}     y/t = 1.0", f_fmt(topy, 8, 4));
    crate::outln!("  Bottom surface:  y ={}     y/t = 0.0", f_fmt(boty, 8, 4));
    if *yfh == -999.0 {
        *yfh = s.input.askr("Enter flap hinge y location (or 999 to specify y/t)^", *yfh);
    }
    if *yfh == 999.0 {
        let yrel = s.input.askr("Enter flap hinge relative y/t location^", 0.0);
        *yfh = topy * yrel + boty * (1.0 - yrel);
    }
    (tops, bots)
}

fn flap(s: &mut Session, rinput: &[f64], ninput: usize) {
    let (mut xbf, mut ybf) = if ninput >= 2 { (rinput[1], rinput[2]) } else { (-999.0, -999.0) };
    let (tops, bots) = getxyf(s, &mut xbf, &mut ybf);
    s.xf.xbf = xbf;
    s.xf.ybf = ybf;
    let insid = inside(&s.xf.xb, &s.xf.yb, s.xf.nb, xbf, ybf);
    crate::outln!("\n Flap hinge: x,y ={}{}", f_fmt(xbf, 9, 5), f_fmt(ybf, 9, 5));

    let ddef = if ninput >= 3 { rinput[3] } else { s.input.askr("Enter flap deflection in degrees (+ down)^", 0.0) };
    let rdef = ddef * PI / 180.0;
    if rdef == 0.0 {
        return;
    }

    let xf = &mut s.xf;
    let mut nb = xf.nb;
    let (atop, abot);
    if insid {
        atop = (-rdef).max(0.0);
        abot = rdef.max(0.0);
    } else {
        let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
        let chx = deval(bots, xb, xbp, sb, nb) - deval(tops, xb, xbp, sb, nb);
        let chy = deval(bots, yb, ybp, sb, nb) - deval(tops, yb, ybp, sb, nb);
        let fvx = seval(bots, xb, xbp, sb, nb) + seval(tops, xb, xbp, sb, nb);
        let fvy = seval(bots, yb, ybp, sb, nb) + seval(tops, yb, ybp, sb, nb);
        let crsp = chx * (ybf - 0.5 * fvy) - chy * (xbf - 0.5 * fvx);
        if crsp > 0.0 {
            atop = rdef.max(0.0);
            abot = rdef.max(0.0);
        } else {
            atop = (-rdef).max(0.0);
            abot = (-rdef).max(0.0);
        }
    }

    let (st1, st2, sb1, sb2, xt1, yt1, xb1, yb1);
    {
        let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
        let r = sss(tops, atop, xbf, ybf, xb, xbp, yb, ybp, sb, nb, 1);
        st1 = r.0;
        st2 = r.1;
        let r = sss(bots, abot, xbf, ybf, xb, xbp, yb, ybp, sb, nb, 2);
        sb1 = r.0;
        sb2 = r.1;
        xt1 = seval(st1, xb, xbp, sb, nb);
        yt1 = seval(st1, yb, ybp, sb, nb);
        let xt2 = seval(st2, xb, xbp, sb, nb);
        let yt2 = seval(st2, yb, ybp, sb, nb);
        xb1 = seval(sb1, xb, xbp, sb, nb);
        yb1 = seval(sb1, yb, ybp, sb, nb);
        let xb2 = seval(sb2, xb, xbp, sb, nb);
        let yb2 = seval(sb2, yb, ybp, sb, nb);
        crate::outln!(
            "\n Top breaks: x,y =  {}{}    {}{}",
            f_fmt(xt1, 9, 5),
            f_fmt(yt1, 9, 5),
            f_fmt(xt2, 9, 5),
            f_fmt(yt2, 9, 5)
        );
        crate::outln!(
            " Bot breaks: x,y =  {}{}    {}{}",
            f_fmt(xb1, 9, 5),
            f_fmt(yb1, 9, 5),
            f_fmt(xb2, 9, 5),
            f_fmt(yb2, 9, 5)
        );
    }

    // find points adjacent to breaks
    let (mut it1, mut it2, mut ib1, mut ib2) = (0usize, 0usize, 0usize, 0usize);
    for i in 1..nb {
        let (si, sip) = (xf.sb[i], xf.sb[i + 1]);
        if si <= st1 && sip > st1 {
            it1 = i + 1;
        }
        if si < st2 && sip >= st2 {
            it2 = i;
        }
        if si <= sb1 && sip > sb1 {
            ib1 = i;
        }
        if si < sb2 && sip >= sb2 {
            ib2 = i + 1;
        }
    }

    let dsavg = (xf.sb[nb] - xf.sb[1]) / (nb - 1) as f64;
    let sfrac = 0.33333;

    let (mut xt1new, mut yt1new, mut xt2new, mut yt2new) = (0.0, 0.0, 0.0, 0.0);
    let (mut xb1new, mut yb1new, mut xb2new, mut yb2new) = (0.0, 0.0, 0.0, 0.0);
    let (mut lt1new, mut lt2new, mut lb1new, mut lb2new) = (false, false, false, false);
    {
        let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
        if atop != 0.0 {
            let st1p = st1 + sfrac * (sb[it1] - st1);
            let st1q = st1 + sfrac * (sb[it1 + 1] - st1);
            if sb[it1] < st1q {
                xt1new = seval(st1q, xb, xbp, sb, nb);
                yt1new = seval(st1q, yb, ybp, sb, nb);
                lt1new = false;
            } else {
                xt1new = seval(st1p, xb, xbp, sb, nb);
                yt1new = seval(st1p, yb, ybp, sb, nb);
                lt1new = true;
            }
            let st2p = st2 + sfrac * (sb[it2] - st2);
            let it2q = it2.saturating_sub(1).max(1);
            let st2q = st2 + sfrac * (sb[it2q] - st2);
            if sb[it2] > st2q {
                xt2new = seval(st2q, xb, xbp, sb, nb);
                yt2new = seval(st2q, yb, ybp, sb, nb);
                lt2new = false;
            } else {
                xt2new = seval(st2p, xb, xbp, sb, nb);
                yt2new = seval(st2p, yb, ybp, sb, nb);
                lt2new = true;
            }
        }
        if abot != 0.0 {
            let sb1p = sb1 + sfrac * (sb[ib1] - sb1);
            let sb1q = sb1 + sfrac * (sb[ib1 - 1] - sb1);
            if sb[ib1] > sb1q {
                xb1new = seval(sb1q, xb, xbp, sb, nb);
                yb1new = seval(sb1q, yb, ybp, sb, nb);
                lb1new = false;
            } else {
                xb1new = seval(sb1p, xb, xbp, sb, nb);
                yb1new = seval(sb1p, yb, ybp, sb, nb);
                lb1new = true;
            }
            let sb2p = sb2 + sfrac * (sb[ib2] - sb2);
            let ib2q = (ib2 + 1).min(nb);
            let sb2q = sb2 + sfrac * (sb[ib2q] - sb2);
            if sb[ib2] < sb2q {
                xb2new = seval(sb2q, xb, xbp, sb, nb);
                yb2new = seval(sb2q, yb, ybp, sb, nb);
                lb2new = false;
            } else {
                xb2new = seval(sb2p, xb, xbp, sb, nb);
                yb2new = seval(sb2p, yb, ybp, sb, nb);
                lb2new = true;
            }
        }
    }

    let sind = rdef.sin();
    let cosd = rdef.cos();

    // rotate flap points about the hinge point
    for i in 1..=nb {
        if i >= it1 && i <= ib1 {
            continue;
        }
        let xbar = xf.xb[i] - xbf;
        let ybar = xf.yb[i] - ybf;
        xf.xb[i] = xbf + xbar * cosd + ybar * sind;
        xf.yb[i] = ybf - xbar * sind + ybar * cosd;
    }

    let idif = it1 as i64 - it2 as i64 - 1;
    if idif > 0 {
        let idif = idif as usize;
        nb -= idif;
        it1 -= idif;
        ib1 -= idif;
        ib2 -= idif;
        for i in it2 + 1..=nb {
            xf.sb[i] = xf.sb[i + idif];
            xf.xb[i] = xf.xb[i + idif];
            xf.yb[i] = xf.yb[i + idif];
        }
    }

    let idif = ib2 as i64 - ib1 as i64 - 1;
    if idif > 0 {
        let idif = idif as usize;
        nb -= idif;
        ib2 -= idif;
        for i in ib1 + 1..=nb {
            xf.sb[i] = xf.sb[i + idif];
            xf.xb[i] = xf.xb[i + idif];
            xf.yb[i] = xf.yb[i + idif];
        }
    }

    if atop == 0.0 {
        let dsnew = rdef.abs() * ((xt1 - xbf).powi(2) + (yt1 - ybf).powi(2)).sqrt();
        let npadd = (1.5 * dsnew / dsavg + 1.0) as i64 as usize;
        if npadd != 0 {
            nb += npadd;
            it1 += npadd;
            ib1 += npadd;
            ib2 += npadd;
            let mut i = nb;
            while i >= it1 {
                xf.xb[i] = xf.xb[i - npadd];
                xf.yb[i] = xf.yb[i - npadd];
                i -= 1;
            }
            let dang = rdef / npadd as f64;
            let xbar = xt1 - xbf;
            let ybar = yt1 - ybf;
            for ip in 1..=npadd {
                let ang = dang * (ip as f64 - 0.5);
                let ca = ang.cos();
                let sa = ang.sin();
                xf.xb[it1 - ip] = xbf + xbar * ca + ybar * sa;
                xf.yb[it1 - ip] = ybf - xbar * sa + ybar * ca;
            }
        }
    } else {
        let mut npadd = 1;
        if lt2new {
            npadd += 1;
        }
        if lt1new {
            npadd += 1;
        }
        nb += npadd;
        it1 += npadd;
        ib1 += npadd;
        ib2 += npadd;
        let mut i = nb;
        while i >= it1 {
            xf.xb[i] = xf.xb[i - npadd];
            xf.yb[i] = xf.yb[i - npadd];
            i -= 1;
        }
        if lt1new {
            xf.xb[it1 - 1] = xt1new;
            xf.yb[it1 - 1] = yt1new;
            xf.xb[it1 - 2] = xt1;
            xf.yb[it1 - 2] = yt1;
        } else {
            xf.xb[it1] = xt1new;
            xf.yb[it1] = yt1new;
            xf.xb[it1 - 1] = xt1;
            xf.yb[it1 - 1] = yt1;
        }
        let xbar = xt2new - xbf;
        let ybar = yt2new - ybf;
        if lt2new {
            xf.xb[it2 + 1] = xbf + xbar * cosd + ybar * sind;
            xf.yb[it2 + 1] = ybf - xbar * sind + ybar * cosd;
        } else {
            xf.xb[it2] = xbf + xbar * cosd + ybar * sind;
            xf.yb[it2] = ybf - xbar * sind + ybar * cosd;
        }
    }

    if abot == 0.0 {
        let dsnew = rdef.abs() * ((xb1 - xbf).powi(2) + (yb1 - ybf).powi(2)).sqrt();
        let npadd = (1.5 * dsnew / dsavg + 1.0) as i64 as usize;
        if npadd != 0 {
            nb += npadd;
            ib2 += npadd;
            let mut i = nb;
            while i >= ib2 {
                xf.xb[i] = xf.xb[i - npadd];
                xf.yb[i] = xf.yb[i - npadd];
                i -= 1;
            }
            let dang = rdef / npadd as f64;
            let xbar = xb1 - xbf;
            let ybar = yb1 - ybf;
            for ip in 1..=npadd {
                let ang = dang * (ip as f64 - 0.5);
                let ca = ang.cos();
                let sa = ang.sin();
                xf.xb[ib1 + ip] = xbf + xbar * ca + ybar * sa;
                xf.yb[ib1 + ip] = ybf - xbar * sa + ybar * ca;
            }
        }
    } else {
        let mut npadd = 1;
        if lb2new {
            npadd += 1;
        }
        if lb1new {
            npadd += 1;
        }
        nb += npadd;
        ib2 += npadd;
        let mut i = nb;
        while i >= ib2 {
            xf.xb[i] = xf.xb[i - npadd];
            xf.yb[i] = xf.yb[i - npadd];
            i -= 1;
        }
        if lb1new {
            xf.xb[ib1 + 1] = xb1new;
            xf.yb[ib1 + 1] = yb1new;
            xf.xb[ib1 + 2] = xb1;
            xf.yb[ib1 + 2] = yb1;
        } else {
            xf.xb[ib1] = xb1new;
            xf.yb[ib1] = yb1new;
            xf.xb[ib1 + 1] = xb1;
            xf.yb[ib1 + 1] = yb1;
        }
        let xbar = xb2new - xbf;
        let ybar = yb2new - ybf;
        if lb2new {
            xf.xb[ib2 - 1] = xbf + xbar * cosd + ybar * sind;
            xf.yb[ib2 - 1] = ybf - xbar * sind + ybar * cosd;
        } else {
            xf.xb[ib2] = xbf + xbar * cosd + ybar * sind;
            xf.yb[ib2] = ybf - xbar * sind + ybar * cosd;
        }
    }

    xf.lgsame = false;

    // check new geometry for splinter segments
    let stol = 0.2;
    scheck(&mut xf.xb, &mut xf.yb, &mut nb, stol);
    xf.nb = nb;

    bupdate(xf);
    xf.lbflap = true;

    if s.gdes_state.lgsym {
        crate::outln!();
        crate::outln!(" Disabling symmetry enforcement");
        s.gdes_state.lgsym = false;
    }
}

// ---------------------------------------------------------------------------
// thickness / camber routines (xtcam.f)
// ---------------------------------------------------------------------------

/// ZERCAM: zero out camber of buffer airfoil.
fn zercam(s: &mut Session) {
    crate::outln!(" Setting current camber to zero.");
    thkcam(&mut s.xf, 1.0, 0.0);
    let xf = &mut s.xf;
    let nb = xf.nb;
    xf.nb = ysym(&mut xf.xb, &mut xf.xbp, &mut xf.yb, &mut xf.ybp, &mut xf.sb, 2 * IQX, nb, 1);
    bgeopar(xf);
}

/// TCBUF: report buffer airfoil thickness and camber.
fn tcbuf(s: &mut Session) {
    getcam_b(&mut s.gdes_state, &s.xf);
    let (cx, cy, tx, ty) = getmax_ct(&mut s.gdes_state);
    print_tc(ty, tx, cy, cx);
}

/// TCSCAL: scale thickness and camber by factors.
fn tcscal(s: &mut Session, rinput: &[f64], ninput: usize) {
    getcam_b(&mut s.gdes_state, &s.xf);
    let (cx, cy, tx, ty) = getmax_ct(&mut s.gdes_state);
    crate::outln!();
    print_tc(ty, tx, cy, cx);
    crate::outln!();
    let lgsym = s.gdes_state.lgsym;
    let tfac = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter new/old thickness scale factor^", 1.0) };
    let cfac = if ninput >= 2 {
        rinput[2]
    } else if lgsym {
        // CFAC is left undefined by the Fortran here; the reference build
        // behaves as CFAC = 0 (which also matches the intent)
        crate::outln!(" Symmetry enforced:  Maintaining zero camber.");
        0.0
    } else {
        s.input.askr("Enter new/old camber    scale factor^", 1.0)
    };
    thkcam(&mut s.xf, tfac, cfac);
    getcam_b(&mut s.gdes_state, &s.xf);
}

/// TCSET: set new max thickness and camber.
fn tcset(s: &mut Session, rinput: &[f64], ninput: usize) {
    getcam_b(&mut s.gdes_state, &s.xf);
    let (cx, cy, tx, ty) = getmax_ct(&mut s.gdes_state);
    crate::outln!();
    print_tc(ty, tx, cy, cx);
    crate::outln!();
    let lgsym = s.gdes_state.lgsym;
    let tnew = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter new max thickness <ret> to skip^", 999.0) };
    let cnew = if ninput >= 2 {
        rinput[2]
    } else if lgsym {
        // CNEW is left undefined by the Fortran here; the reference build
        // behaves as CNEW = 0 (zero camber)
        crate::outln!(" Symmetry enforced:  Maintaining zero camber.");
        0.0
    } else {
        s.input.askr("Enter new max  camber   <ret> to skip^", 999.0)
    };
    let mut cfac = 1.0;
    let mut tfac = 1.0;
    if cy != 0.0 && cnew != 999.0 {
        cfac = cnew / cy;
    }
    if ty != 0.0 && tnew != 999.0 {
        tfac = tnew / (2.0 * ty);
    }
    if tfac.abs() > 100.0 || cfac.abs() > 100.0 {
        crate::outln!("\n Questionable input...");
        crate::outln!(" Implied scaling factors are:{} x thickness", f_fmt(tfac, 13, 2));
        crate::outln!("                             {} x camber   ", f_fmt(cfac, 13, 2));
        let ok = s.input.askl("Apply scaling factors?^");
        if !ok {
            crate::outln!(" No action taken");
            return;
        }
    }
    thkcam(&mut s.xf, tfac, cfac);
    getcam_b(&mut s.gdes_state, &s.xf);
}

/// THKCAM: change buffer airfoil thickness and camber.
fn thkcam(xf: &mut XFoil, tfac: f64, cfac: f64) {
    let nb = xf.nb;
    xf.sble = lefind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
    let sble = xf.sble;
    let xle = seval(sble, &xf.xb, &xf.xbp, &xf.sb, nb);
    let yle = seval(sble, &xf.yb, &xf.ybp, &xf.sb, nb);
    let xte = 0.5 * (xf.xb[1] + xf.xb[nb]);
    let yte = 0.5 * (xf.yb[1] + xf.yb[nb]);
    let chord = ((xte - xle).powi(2) + (yte - yle).powi(2)).sqrt();
    let dxc = (xte - xle) / chord;
    let dyc = (yte - yle) / chord;
    let mut w1 = vec![0.0; nb + 2];
    let mut w2 = vec![0.0; nb + 2];
    {
        let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
        for i in 1..=nb {
            let sbopp = sopps(sb[i], xb, xbp, yb, ybp, sb, nb, sble);
            let xbopp = seval(sbopp, xb, xbp, sb, nb);
            let ybopp = seval(sbopp, yb, ybp, sb, nb);
            let xcavg = 0.5 * (xb[i] + xbopp) * dxc + 0.5 * (yb[i] + ybopp) * dyc;
            let ycavg = cfac * (0.5 * (yb[i] + ybopp) * dxc - 0.5 * (xb[i] + xbopp) * dyc);
            let xcdel = 0.5 * (xb[i] - xbopp) * dxc + 0.5 * (yb[i] - ybopp) * dyc;
            let ycdel = tfac * (0.5 * (yb[i] - ybopp) * dxc - 0.5 * (xb[i] - xbopp) * dyc);
            w1[i] = (xcavg + xcdel) * dxc - (ycavg + ycdel) * dyc;
            w2[i] = (ycavg + ycdel) * dxc + (xcavg + xcdel) * dyc;
        }
    }
    for i in 1..=nb {
        xf.xb[i] = w1[i];
        xf.yb[i] = w2[i];
    }
    xf.lgsame = false;
    bupdate(xf);
}

/// HIPNT: move thickness and camber highpoints.
fn hipnt(s: &mut Session, rinput: &[f64], ninput: usize) {
    {
        let xf = &s.xf;
        let nb = xf.nb;
        let xle = seval(xf.sble, &xf.xb, &xf.xbp, &xf.sb, nb);
        let yle = seval(xf.sble, &xf.yb, &xf.ybp, &xf.sb, nb);
        let xte = 0.5 * (xf.xb[1] + xf.xb[nb]);
        let yte = 0.5 * (xf.yb[1] + xf.yb[nb]);
        let arot = (yle - yte).atan2(xte - xle) / DTOR;
        if arot.abs() > 1.0 {
            crate::outln!("  ");
            crate::outln!(" Warning: HIGH does not work well on rotated foils");
            crate::outln!(" Current chordline angle: {}", lreal(arot));
            crate::outln!(" Proceeding anyway...");
        }
    }
    let (sbl, ybl);
    {
        let xf = &s.xf;
        sbl = xlfind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, xf.nb);
        ybl = seval(sbl, &xf.yb, &xf.ybp, &xf.sb, xf.nb);
    }
    getcam_b(&mut s.gdes_state, &s.xf);
    let (cxmax, cymax, txmax, tymax) = getmax_ct(&mut s.gdes_state);
    crate::outln!();
    print_tc(tymax, txmax, cymax, cxmax);
    crate::outln!();

    let lgsym = s.gdes_state.lgsym;
    let mut thpnt = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter new thickness highpoint x: ^", 0.0) };
    let mut chpnt = 0.0;
    if ninput >= 2 {
        chpnt = rinput[2];
    } else if lgsym {
        crate::outln!(" Symmetry enforced:  Maintaining zero camber.");
    } else {
        chpnt = s.input.askr("Enter new camber    highpoint x: ^", 0.0);
    }
    if thpnt <= 0.0 {
        thpnt = txmax;
    }
    if chpnt <= 0.0 {
        chpnt = cxmax;
    }

    let st = &mut s.gdes_state;
    let mut xfn = [0.0f64; 6];
    let mut yfn = [0.0f64; 6];
    let mut yfnp = [0.0f64; 6];
    if thpnt > 0.0 {
        let ntk = st.ntk;
        xfn[1] = st.xtk[1];
        xfn[2] = txmax;
        xfn[3] = st.xtk[ntk];
        yfn[1] = st.xtk[1];
        yfn[2] = thpnt;
        yfn[3] = st.xtk[ntk];
        splina(&yfn, &mut yfnp, &xfn, 3);
        for i in 1..=ntk {
            st.xtk[i] = seval(st.xtk[i], &yfn, &yfnp, &xfn, 3);
        }
    }
    if chpnt > 0.0 {
        let ncm = st.ncm;
        xfn[1] = st.xcm[1];
        xfn[2] = cxmax;
        xfn[3] = st.xcm[ncm];
        yfn[1] = st.xcm[1];
        yfn[2] = chpnt;
        yfn[3] = st.xcm[ncm];
        splina(&yfn, &mut yfnp, &xfn, 3);
        for i in 1..=ncm {
            st.xcm[i] = seval(st.xcm[i], &yfn, &yfnp, &xfn, 3);
        }
    }
    segspl(&st.ytk, &mut st.ytkp, &st.xtk, st.ntk);
    segspl(&st.ycm, &mut st.ycmp, &st.xcm, st.ncm);

    let xf = &mut s.xf;
    for i in 1..=xf.nb {
        let ycc = seval(xf.xb[i], &st.ycm, &st.ycmp, &st.xcm, st.ncm);
        let ytt = seval(xf.xb[i], &st.ytk, &st.ytkp, &st.xtk, st.ntk);
        xf.yb[i] = if xf.sb[i] <= sbl { ycc + ytt } else { ycc - ytt };
        xf.yb[i] += ybl;
    }
    xf.lgsame = false;
    bupdate(xf);
}

// ---------------------------------------------------------------------------
// CAMB submenu (xtcam.f)
// ---------------------------------------------------------------------------

const CAMB_HELP: &str = "
   <cr>    Return to GDES
   TFAC rr Scale existing thickness and camber
   TSET rr Set new thickness and camber
   HIGH rr Move camber and thickness highpoints
   WRTC    Write airfoil camber x/c,y/c to file

   RDAC    Read   added camber  x/c,y/c from file
   SETC    Set    added camber  x/c,y/c from camberline
   INPC    Input  added camber  x/c,y/c from keyboard
   MODC    Modify added camber  x/c,y/c with cursor
   INPP    Input  added loading x/c,DCp from keyboard
   MODP    Modify added loading x/c,DCp with cursor
   SLOP    Toggle modified-camber,dCp slope matching flag
   SCAL r  Scale the added camber
   CLR     Clear the added camber
   ADD     Add added camber to the existing camberline

   DCPL    Toggle DCp plot
   CPLI rr Change DCp axis plot limits

   Blow    Blowup plot region
   Rese    Reset plot scale and origin

   SIZE r   Change absolute plot-object size
  .ANNO     Annotate plot
   HARD     Hardcopy current plot";

/// Interpolate added camber y(x) to the dense XCAM array.
fn camb_interp_c(st: &mut GdesState) {
    for k in 1..=st.ncam {
        st.ycam[k] = seval(st.xcam[k], &st.ycadd, &st.ycaddp, &st.xcadd, st.ncadd);
        st.ycamp[k] = deval(st.xcam[k], &st.ycadd, &st.ycaddp, &st.xcadd, st.ncadd);
    }
}

/// Added lift and moment from added loading.  Returns (clx, cmx).
fn camb_clcm(st: &GdesState, xcmref: f64) -> (f64, f64) {
    let mut clx = 0.0;
    let mut cmx = 0.0;
    for k in 1..st.ncam {
        let dx = st.xcam[k + 1] - st.xcam[k];
        let xa = 0.5 * (st.xcam[k + 1] + st.xcam[k]);
        let pa = 0.5 * (st.pcam[k + 1] + st.pcam[k]);
        clx += pa * dx;
        cmx += pa * dx * (xcmref - xa);
    }
    (clx, cmx)
}

fn print_clcm(clx: f64, cmx: f64) {
    crate::outln!("\n Delta Cp loading gives delta CL = {}", f_fmt(clx, 7, 3));
    crate::outln!("                        delta CM = {}", f_fmt(cmx, 7, 3));
}

/// Read x/c,value pairs from the keyboard until a blank line.
/// Returns the clamped pairs.
// explicit IF chain as in the Fortran; f64::clamp keeps -0.0 instead of 0.0
#[allow(clippy::manual_clamp)]
fn read_pairs(input: &mut Input) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for _ in 0..2 * IQX {
        loop {
            let line = input.ask("");
            if input.eof || line.trim().is_empty() {
                return out;
            }
            match getflt(&line, 2) {
                Some(v) if v.len() >= 2 => {
                    let xx = if v[0] <= 0.0 {
                        0.0
                    } else if v[0] >= 1.0 {
                        1.0
                    } else {
                        v[0]
                    };
                    out.push((xx, v[1]));
                    break;
                }
                _ => crate::outln!(" try again"),
            }
        }
    }
    out
}

fn camb(s: &mut Session) {
    let lgparsave = s.gdes_state.lgparm;
    let mut lclear = true;
    let mut comold = String::from("****");
    let mut argold = String::new();

    if !s.gdes_state.lplcam {
        crate::outln!(" Enabling camber,thickness plotting");
        s.gdes_state.lplcam = true;
        gofini(s);
    }

    {
        let xf = &s.xf;
        let nb = xf.nb;
        let xle = seval(xf.sble, &xf.xb, &xf.xbp, &xf.sb, nb);
        let yle = seval(xf.sble, &xf.yb, &xf.ybp, &xf.sb, nb);
        let xte = 0.5 * (xf.xb[1] + xf.xb[nb]);
        let yte = 0.5 * (xf.yb[1] + xf.yb[nb]);
        let arot = (yle - yte).atan2(xte - xle) / DTOR;
        if arot.abs() > 1.0 {
            crate::outln!("  ");
            crate::outln!(" Warning: CAMB does not work well on rotated foils");
            crate::outln!(" Current chordline angle: {}", lreal(arot));
            crate::outln!(" Proceeding anyway...");
        }
    }

    s.gdes_state.ldcplot = false;
    s.gdes_state.lgparm = !s.gdes_state.ldcplot;
    crate::outln!("{}", CAMB_HELP);

    let (mut sbl, mut xbch, mut ybch, mut sbch) = (0.0, 1.0, 0.0, 1.0);
    let mut stage = 100;
    loop {
        if stage == 100 {
            {
                let xf = &s.xf;
                let nb = xf.nb;
                sbl = xlfind(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb);
                let xbl = seval(sbl, &xf.xb, &xf.xbp, &xf.sb, nb);
                let xbr = 0.5 * (xf.xb[1] + xf.xb[nb]);
                xbch = xbr - xbl;
                ybch = 0.0;
                sbch = (xbch * xbch + ybch * ybch).sqrt();
            }
            getcam_b(&mut s.gdes_state, &s.xf);
            let st = &mut s.gdes_state;
            crate::outln!(" xc0 xc1{}{}", lreal(st.xcm[1]), lreal(st.xcm[st.ncm]));
            st.ncam = 201.min(NTX);
            let (x0, x1) = (st.xcm[1], st.xcm[st.ncm]);
            for k in 1..=st.ncam {
                st.xcam[k] = x0 + (x1 - x0) * (k - 1) as f64 / (st.ncam - 1) as f64;
            }
            if lclear {
                st.ncadd = 2;
                st.xcadd[1] = st.xcam[1];
                st.xcadd[2] = st.xcam[st.ncam];
                st.ycadd[1] = 0.0;
                st.ycadd[2] = 0.0;
                st.npadd = 2;
                st.xpadd[1] = st.xcam[1];
                st.xpadd[2] = st.xcam[st.ncam];
                st.ypadd[1] = 0.0;
                st.ypadd[2] = 0.0;
                segspl(&st.ycadd, &mut st.ycaddp, &st.xcadd, st.ncadd);
                segspl(&st.ypadd, &mut st.ypaddp, &st.xpadd, st.npadd);
                for k in 1..=st.ncam {
                    st.ycam[k] = seval(st.xcam[k], &st.ycadd, &st.ycaddp, &st.xcadd, st.ncadd);
                    st.ycamp[k] = deval(st.xcam[k], &st.ycadd, &st.ycaddp, &st.xcadd, st.ncadd);
                    st.pcam[k] = seval(st.xcam[k], &st.ypadd, &st.ypaddp, &st.xpadd, st.npadd);
                    st.pcamp[k] = deval(st.xcam[k], &st.ypadd, &st.ypaddp, &st.xpadd, st.npadd);
                }
                lclear = false;
            }
            stage = 200;
        }
        if stage == 200 {
            crate::outln!();
            crate::outln!(" Buffer airfoil thickness and camber:");
            tcbuf(s);
            let st = &s.gdes_state;
            let mut xmx = 0.0;
            let mut ymx: f64 = 0.0;
            for k in 1..=st.ncam {
                if st.ycam[k].abs() > ymx.abs() {
                    xmx = st.xcam[k];
                    ymx = st.ycam[k];
                }
            }
            let chrd = st.xcam[st.ncam] - st.xcam[1];
            let ale = deval(st.xcam[1], &st.ycam, &st.ycamp, &st.xcam, st.ncam).atan() / DTOR;
            let ate = deval(st.xcam[st.ncam], &st.ycam, &st.ycamp, &st.xcam, st.ncam).atan() / DTOR;
            crate::outln!("\n Added camber line incidence at LE = {}  deg.", f_fmt(ale, 6, 2));
            crate::outln!(" Added camber line incidence at TE = {}  deg.", f_fmt(ate, 6, 2));
            crate::outln!(" Max added camber y/c = {}  at x/c = {}", f_fmt(ymx / chrd, 8, 4), f_fmt(xmx / chrd, 7, 3));
            stage = 300;
        }
        if stage == 300 {
            s.gdes_state.lgparm = !s.gdes_state.ldcplot;
            plotg(s);
        }
        stage = 500;

        // ---- menu
        let (mut comand, mut comarg) = s.input.askc("..CAMB");
        if s.input.eof {
            s.gdes_state.lgparm = lgparsave;
            return;
        }
        if comand.starts_with('!') {
            if comold == "****" {
                crate::outln!(" Previous ..CAMB command not valid");
                continue;
            }
            comand = comold.clone();
            comarg = argold.clone();
        }
        if comand.is_empty() {
            s.gdes_state.lgparm = lgparsave;
            return;
        }
        let a = parse_args(&comarg);
        let (rinput, ninput) = (&a.r, a.n);

        match comand.as_str() {
            "?" => crate::outln!("{}", CAMB_HELP),
            "Z" | "U" => {}
            "TFAC" => {
                tcscal(s, rinput, ninput);
                stage = 100;
            }
            "TSET" => {
                tcset(s, rinput, ninput);
                stage = 100;
            }
            "HIGH" => {
                hipnt(s, rinput, ninput);
                stage = 100;
            }
            "WRTC" => {
                let fname = s.input.ask("Enter output camber filename^").trim().to_string();
                if s.input.eof {
                    return;
                }
                if std::path::Path::new(&fname).exists() {
                    crate::outln!();
                    crate::outln!(" Output file exists.  Overwrite?  Y");
                    let ans = s.input.ask("");
                    if ans.starts_with(['N', 'n']) {
                        crate::outln!(" Current camber not saved.");
                        continue;
                    }
                }
                let st = &s.gdes_state;
                let mut txt = format!("Camber: {}\n", fixed(&s.xf.name, 48));
                for k in 1..=st.ncm {
                    txt += &format!(
                        " {} {}\n",
                        f_fmt((st.xcm[k] - st.xcm[1]) / xbch, 12, 6),
                        f_fmt((st.ycm[k] - st.ycm[1]) / xbch, 12, 6)
                    );
                }
                if std::fs::write(&fname, txt).is_err() {
                    crate::outln!(" Error opening camber save file");
                }
            }
            "RDAC" => {
                let fname = s.input.ask("Enter added camber filename^").trim().to_string();
                if s.input.eof {
                    return;
                }
                let text = match std::fs::read_to_string(&fname) {
                    Ok(t) => t,
                    Err(_) => {
                        crate::outln!();
                        crate::outln!(" Error opening added camber file");
                        continue;
                    }
                };
                let st = &mut s.gdes_state;
                st.ncadd = 0;
                for line in text.lines().skip(1).take(NTX) {
                    match getflt(line, 2) {
                        Some(v) if v.len() >= 2 => {
                            st.ncadd += 1;
                            st.xcadd[st.ncadd] = v[0];
                            st.ycadd[st.ncadd] = v[1];
                        }
                        _ => break,
                    }
                }
                if st.ncadd <= 1 || (st.xcadd[st.ncadd] - st.xcadd[1]) == 0.0 {
                    st.ncadd = 2;
                    st.xcadd[1] = st.xcam[1];
                    st.xcadd[2] = st.xcam[st.ncam];
                    st.ycadd[1] = 0.0;
                    st.ycadd[2] = 0.0;
                    crate::outln!(" No added camber points found");
                    stage = 100;
                    continue;
                }
                let xcorg = st.xcadd[1];
                let xcscl = st.xcadd[st.ncadd] - xcorg;
                for k in 1..=st.ncadd {
                    st.xcadd[k] = (st.xcadd[k] - xcorg) / xcscl;
                    st.ycadd[k] /= xcscl;
                }
                for k in 1..=st.ncadd {
                    let (xc, yc) = (st.xcadd[k], st.ycadd[k]);
                    st.xcadd[k] = st.xcam[1] + xc * xbch - yc * ybch;
                    st.ycadd[k] = st.xcadd[k] * ybch + yc * xbch;
                }
                segspl(&st.ycadd, &mut st.ycaddp, &st.xcadd, st.ncadd);
                camb_interp_c(st);
                st.ldcplot = false;
                stage = 200;
            }
            "SETC" => {
                let st = &mut s.gdes_state;
                st.ncadd = st.ncm;
                for k in 1..=st.ncm {
                    st.xcadd[k] = st.xcm[k];
                    st.ycadd[k] = st.ycm[k];
                }
                segspl(&st.ycadd, &mut st.ycaddp, &st.xcadd, st.ncadd);
                camb_interp_c(st);
                st.ldcplot = false;
                stage = 200;
            }
            "INPC" => {
                crate::outln!("\n Manual input of camber x/c,y/c:\n\n Input x/c, y/c pairs from  x/c = 0  to  x/c = 1\n <cr> ends input");
                let pairs = read_pairs(&mut s.input);
                let st = &mut s.gdes_state;
                st.ncadd = 0;
                for (xx, yy) in pairs {
                    st.ncadd += 1;
                    st.xcadd[st.ncadd] = st.xcam[1] + xx * xbch - yy * ybch;
                    st.ycadd[st.ncadd] = xx * ybch + yy * xbch;
                }
                if st.ncadd >= 2 {
                    sortdup(st.ncadd, &mut st.xcadd, &mut st.ycadd);
                    fixdup(&mut st.ncadd, &mut st.xcadd, &mut st.ycadd);
                }
                if st.ncadd < 2 {
                    crate::outln!(" Too few added camber points.  Added camber cleared.");
                    st.ncadd = 2;
                    st.xcadd[1] = st.xcam[1];
                    st.xcadd[2] = st.xcam[st.ncam];
                    st.ycadd[1] = 0.0;
                    st.ycadd[2] = 0.0;
                }
                segspl(&st.ycadd, &mut st.ycaddp, &st.xcadd, st.ncadd);
                camb_interp_c(st);
                st.ldcplot = false;
                stage = 200;
            }
            "INPP" => {
                crate::outln!("\n Manual input of loading x/c, DCp:\n\n Input x/c, DCp pairs from  x/c = 0  to  x/c = 1\n <cr> ends input");
                let pairs = read_pairs(&mut s.input);
                let xcmref = s.xf.xcmref;
                let st = &mut s.gdes_state;
                for (xx, yy) in pairs {
                    st.npadd += 1;
                    st.xpadd[st.npadd] = st.xcam[1] + xx * xbch;
                    st.ypadd[st.npadd] = yy;
                }
                sortdup(st.npadd, &mut st.xpadd, &mut st.ypadd);
                fixdup(&mut st.npadd, &mut st.xpadd, &mut st.ypadd);
                segspl(&st.ypadd, &mut st.ypaddp, &st.xpadd, st.npadd);
                for k in 1..=st.ncam {
                    st.pcam[k] = seval(st.xcam[k], &st.ypadd, &st.ypaddp, &st.xpadd, st.npadd);
                    st.pcamp[k] = deval(st.xcam[k], &st.ypadd, &st.ypaddp, &st.xpadd, st.npadd);
                }
                {
                    let GdesState { ncam, xcam, ycam, ycamp, pcam, pcamp, .. } = st;
                    cpcam(*ncam, xcam, ycam, ycamp, pcam, pcamp);
                }
                let (clx, cmx) = camb_clcm(st, xcmref);
                print_clcm(clx, cmx);
                st.ldcplot = true;
                stage = 200;
            }
            "MODC" | "MODP" => cursor_note(&comand),
            "CLR" => {
                lclear = true;
                s.gdes_state.ldcplot = false;
                stage = 100;
            }
            "SCAL" => {
                let scal = if ninput >= 1 { rinput[1] } else { s.input.askr("Enter camber scaling factor^", 1.0) };
                let xcmref = s.xf.xcmref;
                let st = &mut s.gdes_state;
                for i in 1..=st.ncadd {
                    st.ycadd[i] *= scal;
                    st.ypadd[i] *= scal;
                    st.ycaddp[i] *= scal;
                    st.ypaddp[i] *= scal;
                }
                for i in 1..=st.ncam {
                    st.ycam[i] *= scal;
                    st.ycamp[i] *= scal;
                    st.pcam[i] *= scal;
                    st.pcamp[i] *= scal;
                }
                let (clx, cmx) = camb_clcm(st, xcmref);
                if clx != 0.0 && cmx != 0.0 {
                    print_clcm(clx, cmx);
                }
                comold = comand.clone();
                argold = comarg.clone();
                stage = 200;
            }
            "ADD" => {
                let st = &mut s.gdes_state;
                segspl(&st.ycam, &mut st.ycamp, &st.xcam, st.ncam);
                let xf = &mut s.xf;
                let nb = xf.nb;
                let mut w1 = vec![0.0; nb + 2];
                {
                    let (xb, xbp, yb, ybp, sb) = (&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb);
                    for i in 1..=nb {
                        let sbopp = sopps(sb[i], xb, xbp, yb, ybp, sb, nb, sbl);
                        let xbopp = seval(sbopp, xb, xbp, sb, nb);
                        let ybopp = seval(sbopp, yb, ybp, sb, nb);
                        let oldcam = 0.5 * (yb[i] + ybopp) * xbch / sbch - 0.5 * (xb[i] + xbopp) * ybch / sbch;
                        let cam = oldcam + seval(xb[i], &st.ycam, &st.ycamp, &st.xcam, st.ncam);
                        w1[i] = cam + 0.5 * (yb[i] - ybopp);
                    }
                }
                for i in 1..=nb {
                    xf.yb[i] = w1[i];
                }
                bupdate(xf);
                st.ldcplot = false;
                lclear = true;
                stage = 100;
            }
            "SLOP" => {
                s.gdes_state.lcslop = !s.gdes_state.lcslop;
                if s.gdes_state.lcslop {
                    crate::outln!(" Modified segment will be made tangent at endpoints");
                } else {
                    crate::outln!(" Modified segment will not be made tangent at endpoints");
                }
            }
            "BLOW" | "B" => cursor_note(&comand),
            "RESE" | "R" => {
                gofini(s);
                plotg(s);
                stage = 300;
            }
            "DCPL" => {
                s.gdes_state.ldcplot = !s.gdes_state.ldcplot;
                stage = 200;
            }
            "CPLI" => {
                let (mut ypmn, mut ypmx) = (0.0, 0.0);
                if ninput >= 2 {
                    ypmn = rinput[1];
                    ypmx = rinput[2];
                } else if ninput >= 1 {
                    ypmx = s.input.askr("Enter max DCp^", ypmx);
                } else {
                    ypmn = s.input.askr("Enter min DCp^", ypmn);
                    ypmx = s.input.askr("Enter max DCp^", ypmx);
                }
                if ypmx - ypmn > 0.0 {
                    gofini(s);
                }
                stage = 300;
            }
            "SIZE" => {
                if ninput >= 1 {
                    s.gdes_state.size = rinput[1];
                } else {
                    crate::outln!(" Current plot-object size ={}", lreal(s.gdes_state.size));
                    s.gdes_state.size = s.input.askr("Enter new plot-object size^", s.gdes_state.size);
                }
                gofini(s);
                stage = 300;
            }
            "ANNO" => {
                cursor_note(&comand);
                stage = 300;
            }
            "HARD" => stage = 300,
            _ => {
                crate::outln!(" {} command not recognized.  Type a \"?\" for list", fixed(&comand, 4));
            }
        }
        if s.input.eof {
            s.gdes_state.lgparm = lgparsave;
            return;
        }
    }
}
