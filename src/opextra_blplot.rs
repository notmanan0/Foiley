//! OPER .VPLO sub-menu (BLPLOT, blplot.f) without graphics: keeps the
//! plot-type / axis-limit state so that the menu reads exactly the same
//! input as XFoil, and implements the data commands DUMP and OVER.

use crate::cli::{args, Session};
use crate::fmt::{f_fmt, fixed, g_fmt};
use crate::polar::{eol, i_fmt};
use crate::xfoil::{DTOR, IVX};

const FILSUF: [&str; 13] = ["", "hk", "dT", "dB", "ue", "cf", "cd", "nc", "ct", "rt", "rl", "G2", "be"];

/// SCALIT for a single value
fn scalit1(y: f64) -> f64 {
    let ag2 = 2.0f64.log10();
    let ag5 = 5.0f64.log10();
    let mut ymax = y.abs();
    if ymax == 0.0 {
        ymax = 1.0e-8;
    }
    let ylog = ymax.log10();
    let ylog1 = (ylog + 100.0).trunc() - 99.0;
    let mut ylog2 = ylog1 + ag2;
    if ylog2 - 1.0 > ylog {
        ylog2 -= 1.0;
    }
    let mut ylog5 = ylog1 + ag5;
    if ylog5 - 1.0 > ylog {
        ylog5 -= 1.0;
    }
    let gmin = ylog1.min(ylog2).min(ylog5);
    10f64.powf(-gmin)
}

/// SCALIT over an array y[1..=n] (1-based slice start at index `i0`)
fn scalit(n: usize, y: &[f64], i0: usize) -> f64 {
    let mut ymax: f64 = if n >= 1 { y[i0].abs() } else { 0.0 };
    for i in 1..n {
        ymax = ymax.max(y[i0 + i].abs());
    }
    scalit1(ymax)
}

fn sign(a: f64, b: f64) -> f64 {
    if b >= 0.0 {
        a.abs()
    } else {
        -a.abs()
    }
}

struct Lim {
    xmin: f64,
    xmax: f64,
    xdel: f64,
    y: [[f64; 3]; 13],
}

fn defaults(s: &Session, nside: [usize; 3], numbl: [usize; 3]) -> Lim {
    let xf = &s.xf;
    let chx = 0.99 * (xf.xte - xf.xle);
    let xtew = xf.xte + 0.4 * chx;
    let xfac = scalit1(chx);
    let xdel = 1.0 / (5.0 * xfac);
    let xmax = (xtew.abs() / xdel + 0.05).trunc() * sign(xdel, xtew);
    let xmin = (xf.xle.abs() / xdel + 0.05).trunc() * sign(xdel, xf.xle);
    let mut y = [[0.0; 3]; 13];
    y[1] = [0.0, 6.0, 1.0];
    y[11] = [0.0, 5.0, 1.0];
    y[12] = [-1.0, 5.0, 1.0];
    let yfac = scalit(numbl[1], &xf.dstr[1], 2);
    y[2] = [0.0, 1.0 / yfac, (1.0 / yfac) / 5.0];
    let yfac = scalit(numbl[2], &xf.dstr[2], 2);
    y[3] = [0.0, 1.0 / yfac, (1.0 / yfac) / 5.0];
    y[4] = [0.0, 1.6, 0.2];
    let mut taumax: f64 = 0.0;
    let mut dismax: f64 = 0.0;
    for is in 1..=2 {
        for ibl in 2..=nside[is] {
            taumax = taumax.max(xf.tau[is][ibl]);
            dismax = dismax.max(xf.dis[is][ibl]);
        }
    }
    let que = 0.5 * xf.qinf * xf.qinf;
    let yfac = scalit1(taumax / que);
    y[5] = [0.0, 0.5 / yfac, (0.5 / yfac) / 5.0];
    let qrf = xf.qinf;
    let yfac = scalit1(dismax / qrf.powi(3));
    y[6] = [0.0, 0.5 / yfac, (0.5 / yfac) / 5.0];
    let acr1 = (xf.acrit + 1.5).max(1.0);
    let yfac = scalit1(acr1);
    let andel = 1.0 / (5.0 * yfac);
    y[7] = [0.0, andel * (acr1 / andel + 0.6).trunc(), andel];
    let mut cmax: f64 = 0.0;
    for is in 1..=2 {
        for ibl in xf.itran[is].max(1)..=nside[is] {
            cmax = cmax.max(xf.ctau[is][ibl].abs());
        }
    }
    let yfac = scalit1(cmax);
    y[8] = [0.0, 1.0 / yfac, (1.0 / yfac) / 5.0];
    let mut rmax: f64 = 0.0;
    for is in 1..=2 {
        for ibl in 2..=nside[is] {
            rmax = rmax.max(xf.reinf * xf.uedg[is][ibl] * xf.thet[is][ibl]);
        }
    }
    let yfac = scalit1(rmax);
    y[9] = [0.0, 1.0 / yfac, (1.0 / yfac) / 5.0];
    y[10] = [1.0, 5.0, 1.0];
    Lim { xmin, xmax, xdel, y }
}

/// list-directed READ(*,*) of n reals (may span lines). None on EOF,
/// Some(None) on a conversion error.
pub fn read_n(s: &mut Session, n: usize) -> Option<Option<Vec<f64>>> {
    let mut vals: Vec<f64> = Vec::new();
    loop {
        let l = s.input.read_line();
        if s.input.eof {
            return None;
        }
        let mut slash = false;
        for tok in l.split(|c: char| c == ' ' || c == ',' || c == '\t') {
            if tok.is_empty() {
                continue;
            }
            if vals.len() >= n {
                break;
            }
            let tok = match tok.find('/') {
                Some(k) => {
                    slash = true;
                    &tok[..k]
                }
                None => tok,
            };
            if !tok.is_empty() {
                match crate::polar::parse_fortran_real(tok) {
                    Some(v) => vals.push(v),
                    None => return Some(None),
                }
            }
            if slash {
                break;
            }
        }
        if vals.len() >= n || slash {
            return Some(Some(vals));
        }
    }
}

/// BLPLOT menu
pub fn blplot(s: &mut Session) {
    let mut kplot = 0usize;
    let (nside, numbl) = {
        let xf = &s.xf;
        let ns1 = xf.nbl[2] as i64 + xf.iblte[1] as i64 - xf.iblte[2] as i64;
        let ns1 = ns1.clamp(0, IVX as i64) as usize;
        let ns2 = xf.nbl[2];
        ([0, ns1, ns2], [0, ns1.saturating_sub(1), ns2.saturating_sub(1)])
    };
    // Cartesian BL x-arrays
    let mut xxbl = [vec![0.0; IVX + 2], vec![0.0; IVX + 2], vec![0.0; IVX + 2]];
    {
        let xf = &s.xf;
        for is in 1..=2 {
            for ibl in 2..=xf.nbl[is] {
                xxbl[is][ibl] = xf.x[xf.ipan[is][ibl]];
            }
        }
        for iblw in 1..=xf.nbl[2].saturating_sub(xf.iblte[2]) {
            xxbl[1][xf.iblte[1] + iblw] = xxbl[2][xf.iblte[2] + iblw];
        }
    }
    let mut lim = defaults(s, nside, numbl);
    let (mut ymin, mut ymax, mut ydel) = (0.0, 0.0, 0.0);
    loop {
        let (com, arg) = s.input.askc("..VPLO^");
        if s.input.eof && com == "QUIT" {
            return;
        }
        let (r, _) = args(&arg);
        let set_plot = |k: usize, lim: &Lim| (lim.y[k][0], lim.y[k][1], lim.y[k][2]);
        let plot_k = match com.as_str() {
            "" => return,
            "?" => {
                crate::outln!("\n   <cr>   Return to OPER menu\n\n   H      Plot kinematic shape parameter\n   DT     Plot top    side Dstar and Theta\n   DB     Plot bottom side Dstar and Theta\n   UE     Plot edge velocity\n   CF     Plot skin friction coefficient\n   CD     Plot dissipation coefficient\n   N      Plot amplification ratio\n   CT     Plot max shear coefficient\n   RT     Plot Re_theta\n   RTL    Plot log(Re_theta)\n\n   DUMP f Write current plot variable to file\n   OVER f Overlay current plot variable from file\n\n   X  rrr Change x-axis limits\n   Y  rrr Change y-axis limits on current plot\n\n   BLOW   Cursor blowup of current plot\n   RESE   Reset to default x,y-axis limits\n   SIZE r Change absolute plot-object size\n  .ANNO   Annotate plot\n   HARD   Hardcopy current plot\n\n   GRID   Toggle grid plotting\n   SYMB   Toggle node-symbol plotting\n   LABE   Toggle label plotting\n   CLIP   Toggle line-plot clipping\n   FRPL   Toggle TS frequency plotting");
                0
            }
            "H" | "HK" => 1,
            "DS" | "DT" => 2,
            "DP" | "DB" => 3,
            "UE" => 4,
            "CF" => 5,
            "CD" => 6,
            "N" => 7,
            "CT" => 8,
            "RT" => 9,
            "RTL" => 10,
            "G" => 11,
            "BE" => 12,
            "XLIM" | "X" => {
                let v = if r.len() >= 3 {
                    [r[0], r[1], r[2]]
                } else {
                    loop {
                        crate::out!(
                            "\n Currently, Xmin,Xmax,Xdel ={}{}{}\n Enter new  Xmin,Xmax,Xdel: ",
                            f_fmt(lim.xmin, 11, 4),
                            f_fmt(lim.xmax, 11, 4),
                            f_fmt(lim.xdel, 11, 4)
                        );
                        crate::flush_if_tty();
                        match read_n(s, 3) {
                            None => return,
                            Some(Some(v)) => break [v.first().copied().unwrap_or(lim.xmin), v.get(1).copied().unwrap_or(lim.xmax), v.get(2).copied().unwrap_or(lim.xdel)],
                            Some(None) => {}
                        }
                    }
                };
                lim.xmin = v[0];
                lim.xmax = v[1];
                lim.xdel = v[2];
                kplot
            }
            "YLIM" | "Y" => {
                if kplot == 0 {
                    crate::outln!(" No current plot");
                    0
                } else {
                    let v = if r.len() >= 3 {
                        Some([r[0], r[1], r[2]])
                    } else {
                        loop {
                            crate::out!(
                                "\n Currently, Ymin,Ymax,Ydel ={}{}{}\n Enter new  Ymin,Ymax,Ydel :    ",
                                f_fmt(ymin, 11, 4),
                                f_fmt(ymax, 11, 4),
                                f_fmt(ydel, 11, 4)
                            );
                            crate::flush_if_tty();
                            match read_n(s, 3) {
                                None => return,
                                Some(Some(v)) => break Some([v.first().copied().unwrap_or(ymin), v.get(1).copied().unwrap_or(ymax), v.get(2).copied().unwrap_or(ydel)]),
                                // READ error: the Fortran jumps to the DUMP/OVER code (ERR=140)
                                Some(None) => break None,
                            }
                        }
                    };
                    match v {
                        Some(v) => {
                            lim.y[kplot][0] = v[0];
                            lim.y[kplot][1] = v[1];
                            if kplot != 10 {
                                lim.y[kplot][2] = v[2];
                            }
                            // (the Fortran replots G for KPLOT=12)
                            if kplot == 12 {
                                11
                            } else {
                                kplot
                            }
                        }
                        None => {
                            let fname = dump_fname(s, &arg, kplot);
                            dump(s, kplot, &fname, &xxbl, nside);
                            0
                        }
                    }
                }
            }
            "BLOW" | "B" => {
                if kplot == 0 {
                    crate::outln!(" No current plot");
                } else if s.extra_state.lcurs {
                    crate::outln!();
                    crate::outln!(" Mark off corners of blowup area");
                    crate::outln!(" (2 identical points default to current area)");
                    crate::outln!(" (cursor input not available in this build)");
                } else {
                    crate::outln!();
                    crate::outln!(" Enter x,y coordinates of blowup area corners");
                    crate::outln!(" (2 identical points default to current area)");
                    crate::outln!();
                    for k in 1..=2 {
                        loop {
                            crate::outln!(" Point {}: ", k);
                            match read_n(s, 2) {
                                None => return,
                                Some(Some(_)) => break,
                                Some(None) => {}
                            }
                        }
                    }
                }
                0
            }
            "RESE" | "R" => {
                lim = defaults(s, nside, numbl);
                0
            }
            "GRID" | "SYMB" | "LABE" | "CLIP" | "FRPL" => kplot,
            "HARD" | "Z" | "U" => 0,
            "SIZE" => {
                if !r.is_empty() {
                    s.extra_state.size = r[0];
                } else {
                    crate::outln!(" Current plot-object size ={}", crate::fmt::list_real(s.extra_state.size));
                    let v = s.input.askr("Enter new plot-object size^", s.extra_state.size);
                    s.extra_state.size = v;
                }
                kplot
            }
            "ANNO" => {
                if s.extra_state.lplot {
                    crate::opextra::annot(s);
                } else {
                    crate::outln!(" No active plot to annotate");
                }
                0
            }
            "DUMP" | "OVER" => {
                if kplot == 0 {
                    crate::outln!(" No current plot");
                } else {
                    let fname = dump_fname(s, &arg, kplot);
                    if com == "DUMP" {
                        dump(s, kplot, &fname, &xxbl, nside);
                    } else {
                        over(s, &fname);
                    }
                }
                0
            }
            _ => {
                crate::outln!(" {:<4} command not recognized.  Type a \"?\" for list", com);
                0
            }
        };
        if plot_k > 0 {
            kplot = plot_k;
            let (a, b, c) = set_plot(kplot, &lim);
            ymin = a;
            ymax = b;
            ydel = c;
            s.extra_state.lplot = true;
        }
    }
}

fn dump_fname(s: &mut Session, arg: &str, kplot: usize) -> String {
    if !arg.is_empty() {
        arg.chars().take(64).collect::<String>().trim().to_string()
    } else {
        let prefix = s.xf.prefix.trim().to_string();
        if !prefix.is_empty() {
            let fildef = format!("{}.{}", prefix, FILSUF[kplot]);
            crate::outln!("\n Enter filename:  {}", fixed(&fildef, 80));
            let l = s.input.read_line();
            let l = l.chars().take(64).collect::<String>().trim().to_string();
            if l.is_empty() {
                fildef
            } else {
                l
            }
        } else {
            s.input.ask("Enter filename^").chars().take(64).collect::<String>().trim().to_string()
        }
    }
}

fn dump(s: &mut Session, kplot: usize, fname: &str, xxbl: &[Vec<f64>; 3], nside: [usize; 3]) {
    let xf = &s.xf;
    let gamm1 = xf.gamm1;
    let qinf = xf.qinf;
    let tklam = xf.tklam;
    let hstinv = gamm1 * (xf.minf / qinf).powi(2) / (1.0 + 0.5 * gamm1 * xf.minf * xf.minf);
    let hvrat = 0.35;
    let que = 0.5 * qinf * qinf;
    let mut out = String::new();
    out += &format!("#  {}\n", fixed(&xf.name, 48));
    out += &format!("#  alpha ={}\n", f_fmt(xf.alfa / DTOR, 9, 4));
    out += &format!("#  Mach  ={}\n", f_fmt(xf.minf, 9, 4));
    out += &format!("#  Reyn  ={}\n", i_fmt((xf.reinf + 0.5).trunc() as i64, 9));
    out += &format!("#  Ncrit ={}\n", f_fmt(xf.acrit, 9, 4));
    out += "#\n";
    out += &format!("#    x           {}\n", FILSUF[kplot]);
    let mut ibl1 = [0usize, 2, 2];
    let mut ibl2 = [0usize, nside[1], nside[2]];
    let mut ws = [vec![0.0; IVX + 2], vec![0.0; IVX + 2], vec![0.0; IVX + 2]];
    let mut xs = ws.clone();
    let uc_of = |ue: f64| ue * (1.0 - tklam) / (1.0 - tklam * (ue / qinf).powi(2));
    for is in 1..=2 {
        for ibl in 2..=nside[is] {
            match kplot {
                1 => {
                    let thi = xf.thet[is][ibl];
                    let dsi = xf.dstr[is][ibl];
                    let uc = uc_of(xf.uedg[is][ibl]);
                    let amsq = uc * uc * hstinv / (gamm1 * (1.0 - 0.5 * uc * uc * hstinv));
                    ws[is][ibl] = crate::blsys::hkin(dsi / thi, amsq).0;
                    xs[is][ibl] = xxbl[is][ibl];
                }
                2 if is == 1 => {
                    ibl1 = [0, 2, 2];
                    ibl2 = [0, nside[is], nside[is]];
                    ws[1][ibl] = xf.dstr[is][ibl];
                    ws[2][ibl] = xf.thet[is][ibl];
                    xs[1][ibl] = xxbl[is][ibl];
                    xs[2][ibl] = xxbl[is][ibl];
                }
                3 if is == 2 => {
                    ibl1 = [0, 2, 2];
                    ibl2 = [0, nside[is], nside[is]];
                    ws[1][ibl] = xf.dstr[is][ibl];
                    ws[2][ibl] = xf.thet[is][ibl];
                    xs[1][ibl] = xxbl[is][ibl];
                    xs[2][ibl] = xxbl[is][ibl];
                }
                4 => {
                    ws[is][ibl] = uc_of(xf.uedg[is][ibl]);
                    xs[is][ibl] = xxbl[is][ibl];
                }
                5 => {
                    ws[is][ibl] = xf.tau[is][ibl] / que;
                    xs[is][ibl] = xxbl[is][ibl];
                }
                6 => {
                    ws[is][ibl] = xf.dis[is][ibl] / qinf.powi(3);
                    xs[is][ibl] = xxbl[is][ibl];
                }
                7 => {
                    ibl1[is] = 2;
                    ibl2[is] = xf.itran[is].saturating_sub(1);
                    ws[is][ibl] = xf.ctau[is][ibl];
                    xs[is][ibl] = xxbl[is][ibl];
                }
                8 => {
                    ibl1[is] = xf.itran[is];
                    ibl2[is] = nside[is];
                    ws[is][ibl] = xf.ctau[is][ibl];
                    xs[is][ibl] = xxbl[is][ibl];
                }
                9 | 10 => {
                    let ue = uc_of(xf.uedg[is][ibl]);
                    let herat = (1.0 - 0.5 * hstinv * ue * ue) / (1.0 - 0.5 * hstinv * qinf * qinf);
                    let rhoe = herat.powf(1.0 / gamm1);
                    let amue = (herat.powi(3)).sqrt() * (1.0 + hvrat) / (herat + hvrat);
                    let rtheta = xf.reinf * rhoe * ue * xf.thet[is][ibl] / amue;
                    ws[is][ibl] = if kplot == 9 { rtheta } else { rtheta.max(1.0).log10() };
                    xs[is][ibl] = xxbl[is][ibl];
                }
                11 => {
                    let thi = xf.thet[is][ibl];
                    let dsi = xf.dstr[is][ibl];
                    let uc = uc_of(xf.uedg[is][ibl]);
                    let amsq = uc * uc * hstinv / (gamm1 * (1.0 - 0.5 * uc * uc * hstinv));
                    let hki = crate::blsys::hkin(dsi / thi, amsq).0;
                    let qloc = 0.5 * uc * uc / (1.0 + 0.5 * gamm1 * amsq).powf(1.0 / gamm1);
                    let cf = xf.tau[is][ibl] / qloc;
                    let cflim = cf.max(0.0001);
                    let gacon = xf.bl.par.gacon;
                    ws[is][ibl] = ((hki - 1.0) / (gacon * hki)).powi(2) / (0.5 * cflim);
                    xs[is][ibl] = xxbl[is][ibl];
                }
                12 => {
                    let dso = xf.dstr[is][ibl];
                    let ueo = xf.uedg[is][ibl];
                    let (uem, uep, xim, xip);
                    if ibl == xf.iblte[is] || ibl == nside[is] {
                        uem = xf.uedg[is][ibl - 1];
                        uep = xf.uedg[is][ibl];
                        xim = xf.xssi[is][ibl - 1];
                        xip = xf.xssi[is][ibl];
                    } else if ibl == xf.iblte[is] + 1 {
                        uem = xf.uedg[is][ibl];
                        uep = xf.uedg[is][ibl + 1];
                        xim = xf.xssi[is][ibl];
                        xip = xf.xssi[is][ibl + 1];
                    } else {
                        uem = xf.uedg[is][ibl - 1];
                        uep = xf.uedg[is][ibl + 1];
                        xim = xf.xssi[is][ibl - 1];
                        xip = xf.xssi[is][ibl + 1];
                    }
                    let uco = uc_of(ueo);
                    let ucm = uc_of(uem);
                    let ucp = uc_of(uep);
                    let duds = (ucp - ucm) / (xip - xim);
                    let amsq = uco * uco * hstinv / (gamm1 * (1.0 - 0.5 * uco * uco * hstinv));
                    let qloc = 0.5 * uco * uco / (1.0 + 0.5 * gamm1 * amsq).powf(1.0 / gamm1);
                    let cf = xf.tau[is][ibl] / qloc;
                    let cflim = cf.max(0.0001);
                    ws[is][ibl] = -dso * duds / (uco * 0.5 * cflim);
                    xs[is][ibl] = xxbl[is][ibl];
                }
                _ => {}
            }
        }
    }
    for is in 1..=2 {
        for ibl in ibl1[is]..=ibl2[is] {
            if ibl > IVX {
                break;
            }
            out += &format!(" {}{}\n", g_fmt(xs[is][ibl], 14, 6), g_fmt(ws[is][ibl], 14, 6));
        }
        out += "\n";
    }
    if std::fs::write(fname, eol(&out)).is_err() {
        crate::outln!(" File OPEN error");
    }
}

/// OVER: read an overlay file (echo comment lines), set the default prefix.
fn over(s: &mut Session, fname: &str) {
    let txt = match std::fs::read(fname) {
        Ok(b) => String::from_utf8_lossy(&b).to_string(),
        Err(_) => {
            crate::outln!(" File OPEN error");
            return;
        }
    };
    let mut lines = txt.lines();
    let mut line = lines.next();
    while let Some(l) = line {
        if let Some(rest) = l.strip_prefix('#') {
            let r: String = rest.chars().take(79).collect();
            crate::outln!(" {}", fixed(&r, 79));
            line = lines.next();
        } else {
            break;
        }
    }
    let mut is = 1;
    while let Some(l) = line {
        let first10: String = l.chars().take(10).collect();
        if first10.trim().is_empty() {
            if is == 2 {
                break;
            }
            is += 1;
        } else if crate::geom::getflt(l, 2).map(|v| v.len() < 2).unwrap_or(true) {
            crate::outln!(" File READ error");
            return;
        }
        line = lines.next();
    }
    let kdot = fname.find('.');
    s.xf.prefix = match kdot {
        Some(k) => fname[..k].trim().to_string(),
        None => fname.trim().to_string(),
    };
}
