//! Extra OPER-level and top-level commands: BL/Cp dumps (BLDUMP, CPDUMP),
//! flap hinge moments (MHINGE, GETXYF), VELS, NAMMOD, the polar
//! bookkeeping of xpol.f (PLRSET, PLRINI, PLXINI, PLRADD, PLXADD, PLRSUM,
//! PLRCOP, APCOPY, ...) and polar file commands (PGET, PWRT, PLIS, PDEL,
//! PSOR, PREM, PNAM, ASET, RGET, RDEL), and top-level INTE, INTX, HALF,
//! BEND/BENP, MSAV, WDEF/RDEF (WRTDEF/GETDEF), NAME/NINC, PLOP (OPLSET).

use crate::cli::{args, Session};
use crate::fmt::{f_fmt, fixed, g_fmt};
use crate::polar::*;
use crate::xfoil::{DTOR, IBX};

/// Global polar / plotting state that lives in XFOIL.INC commons.
#[derive(Clone, Debug)]
pub struct ExtraState {
    /// IPOL(1..IPTOT), NIPOL, NIPOL0
    pub ipol: [usize; IPTOT + 1],
    pub nipol: usize,
    pub nipol0: usize,
    /// JPOL(1..JPTOT), NJPOL
    pub jpol: [usize; JPTOT + 1],
    pub njpol: usize,
    pub lcminp: bool,
    pub lhmomp: bool,
    pub lpfile: bool,
    pub lpfilx: bool,
    /// reference polars (NPOLREF = len)
    pub refs: Vec<RefPolar>,
    // ---- plotting parameters (kept only for WDEF/RDEF/PLOP fidelity)
    pub idev: i32,
    pub idevrp: i32,
    pub size: f64,
    pub plotar: f64,
    pub ch: f64,
    pub scrnfr: f64,
    pub xpage: f64,
    pub ypage: f64,
    pub xmarg: f64,
    pub ymarg: f64,
    pub lcurs: bool,
    pub lland: bool,
    pub cpmax: f64,
    pub cpmin: f64,
    pub cpdel: f64,
    pub uemax: f64,
    pub uemin: f64,
    pub xofair: f64,
    pub facair: f64,
    pub uprwt: f64,
    /// CPOLPLF(j, kv): [kv][j], kv = IAL,ICL,ICD,ICM ; j = 1..3
    pub cpolplf: [[f64; 4]; 5],
    pub lppsho: bool,
    pub lcpgrd: bool,
    pub lcpref: bool,
    pub lforef: bool,
    /// VELS point (XXX,YYY)
    pub vels_xy: (f64, f64),
    /// LPLOT (a plot is "active"; only used to decide whether ANNO reads input)
    pub lplot: bool,
}

impl Default for ExtraState {
    fn default() -> Self {
        let mut ipol = [0usize; IPTOT + 1];
        ipol[1] = IAL;
        ipol[2] = ICL;
        ipol[3] = ICD;
        ipol[4] = ICP;
        ipol[5] = ICM;
        let mut jpol = [0usize; JPTOT + 1];
        jpol[1] = JTN;
        let mut cpolplf = [[0.0; 4]; 5];
        cpolplf[ICD] = [0.0, 0.0, 0.04, 0.01];
        cpolplf[ICL] = [0.0, 0.0, 1.5, 0.5];
        cpolplf[ICM] = [0.0, -0.25, 0.0, 0.05];
        cpolplf[IAL] = [0.0, -4.0, 10.0, 2.0];
        ExtraState {
            ipol,
            nipol: 5,
            nipol0: 5,
            jpol,
            njpol: 1,
            lcminp: false,
            lhmomp: false,
            lpfile: false,
            lpfilx: false,
            refs: Vec::new(),
            idev: 1,
            idevrp: 2,
            size: 10.0,
            plotar: 0.55,
            ch: 0.015,
            scrnfr: 0.80,
            xpage: 11.0,
            ypage: 8.5,
            xmarg: 0.0,
            ymarg: 0.0,
            lcurs: true,
            lland: true,
            cpmax: 1.0,
            cpmin: -2.0,
            cpdel: -0.5,
            uemax: 1.8,
            uemin: -1.0,
            xofair: 0.09,
            facair: 0.70,
            uprwt: 0.02,
            cpolplf,
            lppsho: false,
            lcpgrd: false,
            lcpref: false,
            lforef: false,
            vels_xy: (0.0, 0.0),
            lplot: false,
        }
    }
}

// ------------------------------------------------------------------------
// small helpers

/// CHARACTER*48 assignment followed by STRIP.
fn name48(s: &str) -> String {
    s.chars().take(48).collect::<String>().trim_matches(' ').to_string()
}

/// CHARACTER*64 FNAME assignment, stripped.
fn fname64(s: &str) -> String {
    s.chars().take(64).collect::<String>().trim_matches(' ').to_string()
}

/// gfortran list-directed REAL output of a single value (WRITE(*,*) x).
fn lreal(v: f64) -> String {
    crate::fmt::list_real(v)
}

/// Fortran Fw.0 (gfortran prints the decimal point)
fn f0(v: f64, w: usize) -> String {
    let s = format!("{:.0}.", v);
    if s.len() > w {
        "*".repeat(w)
    } else {
        format!("{:>w$}", s, w = w)
    }
}

fn lint(v: i64) -> String {
    format!("{:12}", v)
}

fn lbool(b: bool) -> &'static str {
    if b {
        "T"
    } else {
        "F"
    }
}

/// ASKS-style prompt + raw line (CHARACTER string input).
fn asks(s: &mut Session, prompt: &str) -> String {
    s.input.ask(prompt)
}

/// READ(*,'(A)') without prompt.
fn readln(s: &mut Session) -> String {
    s.input.read_line()
}

fn npol(s: &Session) -> usize {
    s.polars.len()
}

// ------------------------------------------------------------------------
// OPER-level dispatcher

/// OPER-level commands not handled in cli.rs. Return true if handled.
pub fn oper_cmd(s: &mut Session, com: &str, arg: &str) -> bool {
    let (r, iinp) = args(arg);
    let ninput = r.len();
    match com {
        "DUMP" => bldump(s, arg),
        "CPWR" => cpdump(s, arg),
        "FMOM" => {
            mhinge(s);
            let xf = &s.xf;
            crate::outln!(
                "\n Flap hinge x,y :{}{}\n                                           2  2\n Hinge moment/span = {}  x  1/2 rho V  c \n                                           2   \n x-Force     /span = {}  x  1/2 rho V  c \n                                           2   \n y-Force     /span = {}  x  1/2 rho V  c \n",
                f_fmt(xf.xof, 8, 4),
                f_fmt(xf.yof, 8, 4),
                f_fmt(xf.hmom, 8, 6),
                f_fmt(xf.hfx, 8, 6),
                f_fmt(xf.hfy, 8, 6)
            );
        }
        "FNEW" => {
            if ninput >= 2 {
                s.xf.xof = r[0];
                s.xf.yof = r[1];
            } else if ninput >= 1 {
                s.xf.xof = r[0];
                s.xf.yof = -999.0;
            } else {
                s.xf.xof = -999.0;
                s.xf.yof = -999.0;
            }
            s.xf.lflap = false;
        }
        "VELS" => vels(s, &r),
        "CPMN" => {
            let xf = &s.xf;
            if xf.lvisc {
                crate::outln!(
                    "  Minimum Inviscid Cp ={}   at x ={}\n  Minimum Viscous  Cp ={}   at x ={}",
                    f_fmt(xf.cpmni, 8, 4),
                    f_fmt(xf.xcpmni, 8, 4),
                    f_fmt(xf.cpmnv, 8, 4),
                    f_fmt(xf.xcpmnv, 8, 4)
                );
            } else {
                crate::outln!("  Minimum Inviscid Cp ={}   at x ={}", f_fmt(xf.cpmni, 8, 4), f_fmt(xf.xcpmni, 8, 4));
            }
        }
        "CINC" => {
            let e = &mut s.extra_state;
            e.lcminp = !e.lcminp;
            if e.lcminp {
                crate::outln!(" Min Cp will be written to polar save file");
            } else {
                crate::outln!(" Min Cp won't be written to polar save file");
            }
        }
        "HINC" => {
            let e = &mut s.extra_state;
            e.lhmomp = !e.lhmomp;
            if e.lhmomp {
                crate::outln!(" Hinge moment will be written to polar save file");
                if !s.xf.lflap {
                    crate::outln!();
                    crate::outln!(" Note: Flap hinge location not defined");
                    crate::outln!("       Set it with FNEW command");
                }
            } else {
                crate::outln!(" Hinge moment won't be written to polar save file");
            }
        }
        "NAME" => {
            if arg.trim().is_empty() {
                let nm = s.xf.name.clone();
                s.xf.name = nammod(s, &nm, 0, -1);
            } else {
                s.xf.name = name48(arg);
            }
        }
        "NINC" => {
            let nm = s.xf.name.clone();
            s.xf.name = nammod(s, &nm, 1, 1);
        }
        "NDEC" => {
            let nm = s.xf.name.clone();
            s.xf.name = nammod(s, &nm, -1, 1);
        }
        "DAMP" => {
            if s.xf.idamp == 0 {
                s.xf.idamp = 1;
                crate::outln!(" Modified amplification used");
            } else {
                s.xf.idamp = 0;
                crate::outln!(" Original amplification used");
            }
        }
        "SEQP" => {
            let e = &mut s.extra_state;
            e.lppsho = !e.lppsho;
            if e.lppsho {
                crate::outln!(" Polar will be plotted during point sequence");
            } else {
                crate::outln!(" Cp(x) will be plotted during point sequence");
            }
        }
        "GRID" => {
            let e = &mut s.extra_state;
            e.lcpgrd = !e.lcpgrd;
            if e.lcpgrd {
                crate::outln!(" Cp grid overlay enabled");
            } else {
                crate::outln!(" Cp grid overlay disabled");
            }
        }
        "CREF" => {
            let e = &mut s.extra_state;
            e.lcpref = !e.lcpref;
            if e.lcpref {
                crate::outln!(" Reference Cp plotting enabled");
            } else {
                crate::outln!(" Reference Cp plotting disabled");
            }
        }
        "FREF" => {
            let e = &mut s.extra_state;
            e.lforef = !e.lforef;
            if e.lforef {
                crate::outln!(" Reference force plotting enabled");
            } else {
                crate::outln!(" Reference force plotting disabled");
            }
        }
        "SIZE" => {
            if ninput >= 1 {
                s.extra_state.size = r[0];
            } else {
                crate::outln!(" Current plot-object size ={}", lreal(s.extra_state.size));
                let v = s.input.askr("Enter new plot-object size^", s.extra_state.size);
                s.extra_state.size = v;
            }
        }
        "CPMI" => {
            if ninput >= 1 {
                s.extra_state.cpmin = r[0];
            } else {
                crate::outln!(" Current CPmin ={}", lreal(s.extra_state.cpmin));
                let v = s.input.askr("Enter new CPmin^", s.extra_state.cpmin);
                s.extra_state.cpmin = v;
            }
            s.extra_state.cpdel = if s.extra_state.cpmin < -4.01 { -1.0 } else { -0.5 };
        }
        "UEMA" => {
            if ninput >= 1 {
                s.extra_state.uemax = r[0];
            } else {
                crate::outln!(" Current UEmax ={}", lreal(s.extra_state.uemax));
                let v = s.input.askr("Enter new UEMAX^", s.extra_state.uemax);
                s.extra_state.uemax = v;
            }
        }
        "UEMI" => {
            if ninput >= 1 {
                s.extra_state.uemin = r[0];
            } else {
                crate::outln!(" Current UEmin ={}", lreal(s.extra_state.uemin));
                let v = s.input.askr("Enter new UEMIN^", s.extra_state.uemin);
                s.extra_state.uemin = v;
            }
        }
        "BLWT" => {
            if ninput >= 1 {
                s.extra_state.uprwt = r[0];
            } else {
                crate::outln!(" Current u/Qinf profile plot weight ={}", lreal(s.extra_state.uprwt));
                let v = s.input.askr("Enter new plot weight^", s.extra_state.uprwt);
                s.extra_state.uprwt = v;
            }
            s.extra_state.lplot = true;
        }
        "BL" => {
            if !s.xf.lvconv {
                crate::outln!(" Compute valid viscous solution first");
            } else {
                if ninput == 0 {
                    crate::outln!(" Using default number of profiles:{}", lint(21));
                }
                s.extra_state.lplot = true;
            }
        }
        "BLC" => {
            if !s.xf.lvconv {
                crate::outln!(" Compute valid viscous solution first");
            } else {
                s.extra_state.lplot = true;
                crate::outln!(" BLC: cursor input not available in this build");
            }
        }
        "BLF" => crate::outln!(" BLF: not available in this build (plotting/debug command)"),
        "ANNO" => {
            if s.extra_state.lplot {
                annot(s);
            } else {
                crate::outln!(" No active plot to annotate");
            }
        }
        "VPLO" | "VP" => crate::opextra_blplot::blplot(s),
        "PPAX" => polaxi(s),
        "CPX" | "CP" | "UEX" | "UE" | "CPV" | "PPLO" | "PP" => s.extra_state.lplot = true,
        "HARD" => s.extra_state.lplot = false,
        "Z" | "U" => {}
        "APLO" | "AP" => {
            if npol(s) == 0 {
                crate::outln!();
                crate::outln!(" No polars are stored");
            } else {
                s.extra_state.lplot = true;
            }
        }
        // ---------------- polar commands
        "PGET" => pget(s, arg),
        "PWRT" => pwrt(s, &iinp),
        "RGET" => rget(s, arg),
        "RDEL" => rdel(s, &iinp),
        "PSUM" => {
            if npol(s) == 0 {
                crate::outln!();
                crate::outln!(" No polars are stored");
            } else {
                let n = npol(s);
                plrsum(s, 1, n);
            }
        }
        "PLIS" => plis(s, &iinp),
        "PDEL" => pdel(s, &iinp),
        "PSOR" => psor(s, &iinp),
        "ASET" => aset(s, &iinp),
        "PREM" => prem(s, &r, &iinp),
        "PNAM" => pnam(s, &iinp),
        "?" => oper_help(),
        _ => return false,
    }
    true
}

fn oper_help() {
    crate::outln!(
        "\n   <cr>     Return to Top Level\n   !        Redo last ALFA,CLI,CL,ASEQ,CSEQ,VELS\n\n   Visc r   Toggle Inviscid/Viscous mode\n  .VPAR     Change BL parameter(s)\n   Re   r   Change Reynolds number\n   Mach r   Change Mach number\n   Type i   Change type of Mach,Re variation with CL\n   ITER     Change viscous-solution iteration limit\n   INIT     Toggle BL initialization flag\n\n   Alfa r   Prescribe alpha\n   CLI  r   Prescribe inviscid CL\n   Cl   r   Prescribe CL\n   ASeq rrr Prescribe a sequence of alphas\n   CSeq rrr Prescribe a sequence of CLs\n\n   SEQP     Toggle polar/Cp(x) sequence plot display\n   CINC     Toggle  minimum Cp  inclusion in polar\n   HINC     Toggle hinge moment inclusion in polar\n   Pacc i   Toggle auto point accumulation to active polar\n   PGET f   Read new polar from save file\n   PWRT i   Write polar to save file\n   PSUM     Show summary of stored polars\n   PLIS i   List stored polar(s)\n   PDEL i   Delete stored polar\n   PSOR i   Sort stored polar\n   PPlo ii. Plot stored polar(s)\n   APlo ii. Plot stored airfoil(s) for each polar\n   ASET i   Copy stored airfoil into current airfoil\n   PREM ir. Remove point(s) from stored polar\n   PNAM i   Change airfoil name of stored polar\n   PPAX     Change polar plot axis limits\n\n   RGET f   Read new reference polar from file\n   RDEL i   Delete stored reference polar\n\n   GRID     Toggle Cp vs x grid overlay\n   CREF     Toggle reference Cp data overlay\n   FREF     Toggle reference CL,CD.. data display\n\n   CPx      Plot Cp vs x\n   CPV      Plot airfoil with pressure vectors (gee wiz)\n  .VPlo     BL variable plots\n  .ANNO     Annotate current plot\n   HARD     Hardcopy current plot\n   SIZE r   Change plot-object size\n   CPMI r   Change minimum Cp axis annotation\n\n   BL   i   Plot boundary layer velocity profiles\n   BLC      Plot boundary layer velocity profiles at cursor\n   BLWT r   Change velocity profile scale weight\n\n   FMOM     Calculate flap hinge moment and forces\n   FNEW rr  Set new flap hinge point\n   VELS rr  Calculate velocity components at a point\n   DUMP f   Output Ue,Dstar,Theta,Cf vs s,x,y to file\n   CPWR f   Output x vs Cp to file\n   CPMN     Report minimum surface Cp\n   NAME s   Specify new airfoil name\n   NINC     Increment name version number"
    );
}

/// ANNOT (plt_util.f): annotation menu; cursor input is not available,
/// but the menu reads the same input lines as XFoil.
pub fn annot(s: &mut Session) {
    loop {
        crate::outln!("
   C haracters           |       
   S lant characters     |       
   M ath characters      |  plot 
   P oint symbol         |       
   L ine                 |       
   A rrow                |       
 
   W idth of characters  | modify
   T ype of point symbol |       ");
        loop {
            crate::outln!();
            crate::out!(" Select option or <return>:  ");
            crate::flush_if_tty();
            let l = readln(s);
            let opt = l.chars().next().unwrap_or(' ');
            if opt == ' ' || s.input.eof {
                return;
            }
            match opt {
                'C' | 'c' | 'S' | 's' | 'M' | 'm' => {
                    crate::outln!(" Click on lower left point of character string...");
                    crate::out!(" Enter character string:  ");
                    let _ = readln(s);
                }
                'P' | 'p' => crate::outln!(" Click on symbol locations ..."),
                'L' | 'l' | 'A' | 'a' => crate::outln!(" Click on line points, twice on last point..."),
                'W' | 'w' => loop {
                    crate::out!(" Enter new character width factor (currently =  1.000):  ");
                    let l = readln(s);
                    if s.input.eof || crate::polar::list_real1(&l).is_some() {
                        break;
                    }
                },
                'T' | 't' => loop {
                    crate::outln!();
                    for l in [
                        "  0  square         7  Y         ",
                        "  1  circle         8  flipped Y ",
                        "  2  triangle       9  *         ",
                        "  3  +             10  flipped * ",
                        "  4  x             11  hourglass ",
                        "  5  diamond       12  bowtie    ",
                        "  6  yield sign    13  star      ",
                    ] {
                        crate::outln!(" {}", l);
                    }
                    crate::outln!();
                    crate::out!(" Enter new symbol type (currently = 0):  ");
                    let l = readln(s);
                    if s.input.eof || list_int1(&l).is_some() {
                        break;
                    }
                },
                _ => break,
            }
        }
    }
}

/// POLAXI: polar plot axis limits.
fn polaxi(s: &mut Session) {
    const CVAR: [&str; 5] = ["", "Alpha", "  CL ", "  CD ", " -CM "];
    crate::outln!(" Enter new axis annotations, or <return> to leave unchanged...");
    crate::outln!();
    for kv in 1..=4 {
        loop {
            let c = s.extra_state.cpolplf[kv];
            crate::outln!("   {}  min, max, delta:{}{}{}", CVAR[kv], f_fmt(c[1], 11, 5), f_fmt(c[2], 11, 5), f_fmt(c[3], 11, 5));
            let l = readln(s);
            if s.input.eof {
                return;
            }
            match crate::geom::getflt(&l, 3) {
                Some(v) => {
                    for (j, x) in v.iter().enumerate() {
                        s.extra_state.cpolplf[kv][j + 1] = *x;
                    }
                    break;
                }
                None => crate::outln!(" READ error.  Enter again."),
            }
        }
    }
}

// ------------------------------------------------------------------------
// NAMMOD

/// NAMMOD: returns new name.
pub fn nammod(s: &mut Session, name: &str, kdel: i64, kmod0: i64) -> String {
    let mut nm = name48(name);
    let kb1 = index(&nm, "[");
    let kb2 = index(&nm, "]");
    let mut namdef = nm.clone();
    let mut kmodp = 0i64;
    let mut ok = true;
    if kb1 != 0 && kb2 != 0 && kb2 as i64 - kb1 as i64 > 1 {
        match list_int1(&substr(&nm, kb1 + 1, kb2 - 1)) {
            Some(k) => {
                let kmod = k.abs();
                kmodp = (kmod + kdel) % 100;
                nm = if kb1 >= 2 { substr(&nm, 1, kb1 - 1).trim_matches(' ').to_string() } else { String::new() };
            }
            None => ok = false,
        }
    } else if kmod0 > 0 {
        kmodp = kmod0 % 100;
    } else {
        kmodp = 0;
    }
    if ok {
        if kmodp >= 10 {
            namdef = fixed(&format!("{} [{:2}]", nm, kmodp), 48).trim_end().to_string();
        } else if kmodp >= 1 {
            namdef = fixed(&format!("{} [{:1}]", nm, kmodp), 48).trim_end().to_string();
        }
    }
    crate::outln!("\n Enter airfoil name or <return> for default:  {}", fixed(&namdef, 48));
    let l = readln(s);
    let l = name48(&l);
    if l.is_empty() {
        namdef
    } else {
        l
    }
}

// ------------------------------------------------------------------------
// BLDUMP / CPDUMP

fn dump_fname(s: &mut Session, arg: &str, suffix: &str) -> String {
    if !arg.starts_with(' ') && !arg.is_empty() {
        return fname64(arg);
    }
    let prefix = s.xf.prefix.trim().to_string();
    if !prefix.is_empty() {
        let fildef = format!("{}{}", prefix, suffix);
        let fildef80 = fixed(&fildef, 80);
        crate::outln!("\n Enter filename:  {}", fildef80);
        let l = readln(s);
        let l = fname64(&l);
        if l.is_empty() {
            fildef
        } else {
            l
        }
    } else {
        fname64(&asks(s, "Enter filename^"))
    }
}

fn delim_of(kdelim: i32) -> char {
    match kdelim {
        0 => ' ',
        1 => ',',
        2 => '\t',
        _ => {
            crate::outln!(" ? Illegal delimiter.  Using blank.");
            ' '
        }
    }
}

/// BSTRIP: remove all blanks.
fn bstrip(l: &str) -> String {
    l.chars().filter(|&c| c != ' ').collect()
}

/// BLDUMP
fn bldump(s: &mut Session, arg: &str) {
    let delim = delim_of(s.xf.kdelim);
    let fname = dump_fname(s, arg, ".bl");
    let xf = &mut s.xf;
    let mut out = String::new();
    if xf.kdelim == 0 {
        out += "#    s        x        y     Ue/Vinf    Dstar     Theta      Cf       H\n";
    } else {
        out += &format!("#s{d}x{d}y{d}Ue/Vinf{d}Dstar{d}Theta{d}Cf{d}H\n", d = delim);
    }
    xf.comset();
    let gamm1 = xf.gamm1;
    let hstinv = gamm1 * (xf.minf / xf.qinf).powi(2) / (1.0 + 0.5 * gamm1 * xf.minf * xf.minf);
    let n = xf.n;
    let line = |s_: f64, x: f64, y: f64, ue: f64, ds: f64, th: f64, cf: f64, hk: f64| -> String {
        if xf.kdelim == 0 {
            format!(
                " {}{}{}{}{}{}{}{}\n",
                f_fmt(s_, 9, 5),
                f_fmt(x, 9, 5),
                f_fmt(y, 9, 5),
                f_fmt(ue, 9, 5),
                f_fmt(ds, 10, 6),
                f_fmt(th, 10, 6),
                f_fmt(cf, 10, 6),
                f_fmt(hk, 10, 3)
            )
        } else {
            let l = format!(
                " {}{d}{}{d}{}{d}{}{d}{}{d}{}{d}{}{d}{}",
                f_fmt(s_, 9, 5),
                f_fmt(x, 9, 5),
                f_fmt(y, 9, 5),
                f_fmt(ue, 9, 5),
                f_fmt(ds, 10, 6),
                f_fmt(th, 10, 6),
                f_fmt(cf, 10, 6),
                f_fmt(hk, 10, 3),
                d = delim
            );
            format!("{}\n", bstrip(&l))
        }
    };
    let tklam = xf.tklam;
    let qinf = xf.qinf;
    for i in 1..=n {
        let is = if xf.gam[i] < 0.0 { 2 } else { 1 };
        let (ds, th, cf, h);
        if xf.lipan && xf.lvisc {
            let ibl = if is == 1 { xf.iblte[is] + 1 - i } else { xf.iblte[is] + i - n };
            ds = xf.dstr[is][ibl];
            th = xf.thet[is][ibl];
            cf = xf.tau[is][ibl] / (0.5 * qinf * qinf);
            h = if th == 0.0 { 1.0 } else { ds / th };
        } else {
            ds = 0.0;
            th = 0.0;
            cf = 0.0;
            h = 1.0;
        }
        let g = xf.gam[i] / qinf;
        let ue = g * (1.0 - tklam) / (1.0 - tklam * g * g);
        let amsq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
        let (hk, _, _) = crate::blsys::hkin(h, amsq);
        out += &line(xf.s[i], xf.x[i], xf.y[i], ue, ds, th, cf, hk);
    }
    if xf.lwake {
        let is = 2;
        for i in n + 1..=n + xf.nw {
            let ibl = xf.iblte[is] + i - n;
            let ds = xf.dstr[is][ibl];
            let th = xf.thet[is][ibl];
            let h = ds / th;
            let cf = 0.0;
            let ui = xf.uedg[is][ibl];
            let g = ui / qinf;
            let ue = g * (1.0 - tklam) / (1.0 - tklam * g * g);
            let amsq = ue * ue * hstinv / (gamm1 * (1.0 - 0.5 * ue * ue * hstinv));
            let (hk, _, _) = crate::blsys::hkin(h, amsq);
            out += &line(xf.s[i], xf.x[i], xf.y[i], ue, ds, th, cf, hk);
        }
    }
    write_file(&fname, eol(&out).as_bytes());
}

/// CPDUMP
fn cpdump(s: &mut Session, arg: &str) {
    let delim = delim_of(s.xf.kdelim);
    let fname = dump_fname(s, arg, ".cp");
    let xf = &mut s.xf;
    let mut out = String::new();
    if xf.kdelim == 0 {
        out += "#      x          Cp  \n";
    } else {
        out += &format!("#x{}Cp\n", delim);
    }
    xf.comset();
    let beta = (1.0 - xf.minf * xf.minf).sqrt();
    let bfac = 0.5 * xf.minf * xf.minf / (1.0 + beta);
    for i in 1..=xf.n {
        let g = xf.gam[i] / xf.qinf;
        let cpinc = 1.0 - g * g;
        let den = beta + bfac * cpinc;
        let cpcom = cpinc / den;
        if xf.kdelim == 0 {
            out += &format!(" {}{}\n", f_fmt(xf.x[i], 11, 5), f_fmt(cpcom, 11, 5));
        } else {
            let l = format!(" {}{}{}", f_fmt(xf.x[i], 11, 5), delim, f_fmt(cpcom, 11, 5));
            out += &bstrip(&l);
            out.push('\n');
        }
    }
    write_file(&fname, eol(&out).as_bytes());
}

fn write_file(fname: &str, data: &[u8]) -> bool {
    let f = if fname.is_empty() { "fort.19" } else { fname };
    match std::fs::write(f, data) {
        Ok(_) => true,
        Err(_) => {
            crate::outln!(" File OPEN error: {}", f);
            false
        }
    }
}

// ------------------------------------------------------------------------
// MHINGE / GETXYF / VELS

/// GETXYF: returns (tops, bots, xf, yf)
fn getxyf(s: &mut Session, mut xff: f64, mut yff: f64) -> (f64, f64, f64, f64) {
    if xff == -999.0 {
        xff = s.input.askr("Enter flap hinge x location^", xff);
    }
    let xf = &s.xf;
    let n = xf.n;
    let mut tops = xf.s[1] + (xf.x[1] - xff);
    let mut bots = xf.s[n] - (xf.x[n] - xff);
    crate::spline::sinvrt(&mut tops, xff, &xf.x, &xf.xp, &xf.s, n);
    crate::spline::sinvrt(&mut bots, xff, &xf.x, &xf.xp, &xf.s, n);
    let topy = crate::spline::seval(tops, &xf.y, &xf.yp, &xf.s, n);
    let boty = crate::spline::seval(bots, &xf.y, &xf.yp, &xf.s, n);
    crate::outln!(
        "\n  Top    surface:  y ={}     y/t = 1.0\n  Bottom surface:  y ={}     y/t = 0.0",
        f_fmt(topy, 8, 4),
        f_fmt(boty, 8, 4)
    );
    if yff == -999.0 {
        yff = s.input.askr("Enter flap hinge y location (or 999 to specify y/t)^", yff);
    }
    if yff == 999.0 {
        let yrel = s.input.askr("Enter flap hinge relative y/t location^", 0.0);
        yff = topy * yrel + boty * (1.0 - yrel);
    }
    (tops, bots, xff, yff)
}

/// MHINGE: flap hinge moment about (XOF,YOF).
pub fn mhinge(s: &mut Session) {
    let (tops, bots);
    if !s.xf.lflap {
        let (t, b, xo, yo) = getxyf(s, s.xf.xof, s.xf.yof);
        tops = t;
        bots = b;
        s.xf.xof = xo;
        s.xf.yof = yo;
        s.xf.lflap = true;
    } else {
        let xf = &s.xf;
        let n = xf.n;
        let mut t = xf.xof;
        let mut b = xf.s[n] - xf.xof;
        crate::spline::sinvrt(&mut t, xf.xof, &xf.x, &xf.xp, &xf.s, n);
        crate::spline::sinvrt(&mut b, xf.xof, &xf.x, &xf.xp, &xf.s, n);
        tops = t;
        bots = b;
    }
    let xf = &mut s.xf;
    let n = xf.n;
    let (x, y, sarr) = (&xf.x, &xf.y, &xf.s);
    let topx = crate::spline::seval(tops, x, &xf.xp, sarr, n);
    let topy = crate::spline::seval(tops, y, &xf.yp, sarr, n);
    let botx = crate::spline::seval(bots, x, &xf.xp, sarr, n);
    let boty = crate::spline::seval(bots, y, &xf.yp, sarr, n);
    let cp: &Vec<f64> = if xf.lvisc { &xf.cpv } else { &xf.cpi };
    let (xof, yof) = (xf.xof, xf.yof);
    let mut hmom = 0.0;
    let mut hfx = 0.0;
    let mut hfy = 0.0;
    for i in 2..=n {
        if sarr[i - 1] >= tops && sarr[i] <= bots {
            continue;
        }
        let dx = x[i] - x[i - 1];
        let dy = y[i] - y[i - 1];
        let xmid = 0.5 * (x[i] + x[i - 1]) - xof;
        let ymid = 0.5 * (y[i] + y[i - 1]) - yof;
        let pmid = 0.5 * (cp[i] + cp[i - 1]);
        hmom += pmid * (xmid * dx + ymid * dy);
        hfx -= pmid * dy;
        hfy += pmid * dx;
    }
    let mut i = 2;
    while i <= n {
        if sarr[i] > tops {
            break;
        }
        i += 1;
    }
    let dx = topx - x[i - 1];
    let dy = topy - y[i - 1];
    let xmid = 0.5 * (topx + x[i - 1]) - xof;
    let ymid = 0.5 * (topy + y[i - 1]) - yof;
    let frac = if sarr[i] != sarr[i - 1] { (tops - sarr[i - 1]) / (sarr[i] - sarr[i - 1]) } else { 0.0 };
    let topp = cp[i] * frac + cp[i - 1] * (1.0 - frac);
    let mut pmid = 0.5 * (topp + cp[i - 1]);
    hmom += pmid * (xmid * dx + ymid * dy);
    hfx -= pmid * dy;
    hfy += pmid * dx;

    let dx = xof - topx;
    let dy = yof - topy;
    let xmid = 0.5 * (topx + xof) - xof;
    let ymid = 0.5 * (topy + yof) - yof;
    hmom += pmid * (xmid * dx + ymid * dy);
    hfx -= pmid * dy;
    hfy += pmid * dx;

    let mut i = n;
    while i >= 2 {
        if sarr[i - 1] < bots {
            break;
        }
        i -= 1;
    }
    let dx = x[i] - botx;
    let dy = y[i] - boty;
    let xmid = 0.5 * (botx + x[i]) - xof;
    let ymid = 0.5 * (boty + y[i]) - yof;
    let frac = if sarr[i] != sarr[i - 1] { (bots - sarr[i - 1]) / (sarr[i] - sarr[i - 1]) } else { 0.0 };
    let botp = cp[i] * frac + cp[i - 1] * (1.0 - frac);
    pmid = 0.5 * (botp + cp[i]);
    hmom += pmid * (xmid * dx + ymid * dy);
    hfx -= pmid * dy;
    hfy += pmid * dx;

    let dx = botx - xof;
    let dy = boty - yof;
    let xmid = 0.5 * (botx + xof) - xof;
    let ymid = 0.5 * (boty + yof) - yof;
    hmom += pmid * (xmid * dx + ymid * dy);
    hfx -= pmid * dy;
    hfy += pmid * dx;

    let dx = x[1] - x[n];
    let dy = y[1] - y[n];
    let xmid = 0.5 * (x[1] + x[n]) - xof;
    let ymid = 0.5 * (y[1] + y[n]) - yof;
    let pmid = 0.5 * (cp[1] + cp[n]);
    hmom += pmid * (xmid * dx + ymid * dy);
    hfx -= pmid * dy;
    hfy += pmid * dx;

    xf.hmom = hmom;
    xf.hfx = hfx;
    xf.hfy = hfy;
}

/// VELS command
fn vels(s: &mut Session, r: &[f64]) {
    let (xxx, yyy);
    if r.len() >= 2 {
        xxx = r[0];
        yyy = r[1];
    } else if r.len() == 1 {
        xxx = r[0];
        yyy = s.input.askr("Enter y^", s.extra_state.vels_xy.1);
    } else {
        xxx = s.input.askr("Enter x^", s.extra_state.vels_xy.0);
        yyy = s.input.askr("Enter y^", s.extra_state.vels_xy.1);
    }
    s.extra_state.vels_xy = (xxx, yyy);
    let mut o = crate::xfoil::PsiOut::new();
    s.xf.psilin(0, xxx, yyy, -1.0, 0.0, false, true, &mut o);
    let vvv = o.psi_ni;
    s.xf.psilin(0, xxx, yyy, 0.0, 1.0, false, true, &mut o);
    let uuu = o.psi_ni;
    let qqq = (uuu * uuu + vvv * vvv).sqrt();
    let cpp = 1.0 - (uuu * uuu + vvv * vvv);
    crate::outln!(
        "\n u/Uinf = {}   v/Uinf = {}\n q/Uinf = {}   Cp     = {}\n",
        f_fmt(uuu, 8, 4),
        f_fmt(vvv, 8, 4),
        f_fmt(qqq, 8, 4),
        f_fmt(cpp, 8, 4)
    );
}

// ------------------------------------------------------------------------
// Polar bookkeeping (xpol.f)

fn polwrit_cur(s: &Session, ip: usize, lhead: bool, ia1: i64, ia2: i64, code: &str, vers: f64) -> String {
    let e = &s.extra_state;
    polwrit(&s.polars[ip - 1], lhead, ia1, ia2, &e.ipol[1..], e.nipol, &e.jpol[1..], e.njpol, code, vers)
}

/// set NIPOL with IPOL(NIPOL) = IMC/ICH  (PLRSET/PLRINI/PLRADD convention)
fn set_nipol_a(s: &mut Session) {
    let e = &mut s.extra_state;
    e.nipol = e.nipol0;
    if e.lcminp {
        e.nipol += 1;
        e.ipol[e.nipol] = IMC;
    }
    if e.lhmomp {
        e.nipol += 1;
        e.ipol[e.nipol] = ICH;
    }
}

/// set NIPOL with IPOL(IMC/ICH) = NIPOL  (the OPER PWRT/PLIS/PREM/PNAM convention)
fn set_nipol_b(s: &mut Session) {
    let e = &mut s.extra_state;
    e.nipol = e.nipol0;
    if e.lcminp {
        e.nipol += 1;
        e.ipol[IMC] = e.nipol;
    }
    if e.lhmomp {
        e.nipol += 1;
        e.ipol[ICH] = e.nipol;
    }
}

/// PLRSUM
pub fn plrsum(s: &Session, ip1: usize, ip2: usize) {
    crate::outln!();
    crate::outln!(" {}{}", "       airfoil                    Re           Mach     ", "  Ncrit  XtripT  XtripB     file");
    crate::outln!(" {}{}", "      ------------------------  ------------  ----------", "  -----  ------  ------    -------------------");
    const CLTYP: [&str; 4] = ["     ", "     ", "/sqCL", "/CL  "];
    // IEXP keeps its previous value when Re = 0 (undefined at first: '*')
    let mut iexp: i64 = 99;
    for ip in ip1..=ip2 {
        if ip < 1 || ip > s.polars.len() {
            continue;
        }
        let p = &s.polars[ip - 1];
        let (cacc, cfil) = if ip == s.ipact { (">", if s.extra_state.lpfile { ">" } else { " " }) } else { (" ", " ") };
        let rman = if p.reyn > 0.0 {
            iexp = (p.reyn.log10().trunc() as i64).clamp(0, 9);
            p.reyn / 10f64.powi(iexp as i32)
        } else {
            0.0
        };
        let ct = |t: i32| if (1..=3).contains(&t) { CLTYP[t as usize] } else { "     " };
        crate::outln!(
            " {}{}  {}{}e{}{}{}{}{}{}{}  {} {}",
            cacc,
            i_fmt(ip as i64, 3),
            fixed(&p.name, 24),
            f_fmt(rman, 7, 3),
            i_fmt(iexp, 1),
            ct(p.retyp),
            f_fmt(p.mach, 7, 3),
            ct(p.matyp),
            f_fmt(p.acrit, 7, 2),
            f_fmt(p.xstrip[1], 8, 3),
            f_fmt(p.xstrip[2], 8, 3),
            cfil,
            p.pfname.trim()
        );
    }
}

/// PRFSUM
fn prfsum(s: &Session) {
    crate::outln!();
    crate::outln!(" {}", "       reference polar                          ");
    crate::outln!(" {}", "      ------------------------------------------");
    for (k, r) in s.extra_state.refs.iter().enumerate() {
        crate::outln!("  {}  {}", i_fmt(k as i64 + 1, 3), fixed(&r.name, 48));
    }
}

/// APCOPY
fn apcopy(s: &mut Session, ip: usize) {
    let p = s.polars[ip - 1].clone();
    if p.nxy < 2 {
        crate::outln!(" No airfoil is stored with this polar");
        return;
    }
    let xf = &mut s.xf;
    xf.n = p.nxy;
    for i in 1..=p.nxy {
        xf.x[i] = p.xyx[i];
        xf.y[i] = p.xyy[i];
    }
    xf.name = p.name.clone();
    xf.finish_current_airfoil();
    xf.invalidate_geometry();
}

/// PLRSET: returns the (possibly zeroed) polar index.
fn plrset(s: &mut Session, ip: usize) -> usize {
    if ip == 0 {
        return 0;
    }
    if ip <= npol(s) {
        crate::outln!();
        crate::outln!(" Existing stored polar is chosen for appending...");
        set_nipol_a(s);
        let na = s.polars[ip - 1].na() as i64;
        let t = polwrit_cur(s, ip, true, 1, na, "XFOIL", VERSION);
        crate::out!("{}", t);
        s.extra_state.nipol = s.extra_state.nipol0;
        // check if geometries differ
        let mut differ = false;
        {
            let xf = &s.xf;
            let p = &s.polars[ip - 1];
            if xf.n != p.nxy {
                differ = true;
            } else {
                let sizref = xf.s[xf.n] - xf.s[1];
                for i in 1..=xf.n {
                    let dsq = (xf.x[i] - p.xyx[i]).powi(2) + (xf.y[i] - p.xyy[i]).powi(2);
                    let dsfrac = dsq.sqrt() / sizref;
                    if dsfrac > 0.00001 {
                        differ = true;
                        break;
                    }
                }
            }
        }
        if differ {
            crate::outln!(" Current airfoil differs from airfoil of stored polar");
            crate::outln!(
                "\n   - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n    0  abort polar accumulation\n    1  compute with current airfoil\n    2  compute with stored  airfoil (overwrite current airfoil)"
            );
            let iopt = s.input.aski("   Select action^", 1);
            if iopt == 0 || s.input.eof {
                return 0;
            } else if iopt == 2 {
                apcopy(s, ip);
            }
        }
        crate::outln!();
        crate::outln!(" Setting current parameters to those of stored polar");
        let p = s.polars[ip - 1].clone();
        let xf = &mut s.xf;
        xf.name = p.name.clone();
        xf.retyp = p.retyp;
        xf.matyp = p.matyp;
        xf.minf1 = p.mach;
        xf.reinf1 = p.reyn;
        xf.acrit = p.acrit;
        xf.xstrip[1] = p.xstrip[1];
        xf.xstrip[2] = p.xstrip[2];
        ip
    } else {
        let xf = &s.xf;
        let mut p = Polar {
            name: name48(&xf.name),
            retyp: xf.retyp,
            matyp: xf.matyp,
            reyn: if xf.lvisc { xf.reinf1 } else { 0.0 },
            mach: xf.minf1,
            acrit: xf.acrit,
            nxy: xf.n,
            ..Default::default()
        };
        p.xstrip[1] = xf.xstrip[1];
        p.xstrip[2] = xf.xstrip[2];
        for i in 1..=xf.n {
            p.xyx[i] = xf.x[i];
            p.xyy[i] = xf.y[i];
        }
        crate::outln!("\n Polar{} newly created for accumulation\n Airfoil archived with polar: {}", i_fmt(ip as i64, 3), fixed(&p.name, 48));
        s.polars.push(p);
        s.polars.len()
    }
}

/// PLRINI: polar save file set-up.
fn plrini(s: &mut Session, ip: usize) {
    let pfname = s.polars[ip - 1].pfname.trim().to_string();
    let npf = pfname.len();
    let prompt = if npf == 0 {
        "Enter  polar save filename  OR  <return> for no file^".to_string()
    } else {
        crate::outln!(" Default polar save filename:  {}", pfname);
        "Enter  new filename  OR  \"none\"  OR  <return> for default^".to_string()
    };
    let mut fname = fname64(&asks(s, &prompt));
    let mut nfn = fname.len();
    if nfn == 0 {
        fname = pfname.clone();
        nfn = npf;
    } else if "NONEnone".contains(substr(&fname, 1, 4).as_str()) {
        nfn = 0;
    }
    if nfn == 0 {
        s.extra_state.lpfile = false;
        crate::outln!();
        crate::outln!(" Polar save file will NOT be written");
        return;
    }
    s.extra_state.lpfile = false;

    if let Ok(bytes) = std::fs::read(&fname) {
        // ---- existing file: read it in
        let txt = String::from_utf8_lossy(&bytes).to_string();
        let err = polread(&txt, &mut s.polars[ip - 1]);
        if err {
            crate::outln!(" Old polar save file READ error");
            return;
        }
        s.polars[ip - 1].pfname = fname.clone();
        let p = s.polars[ip - 1].clone();
        let xf = &s.xf;
        let namdif = name48(&xf.name) != name48(&p.name);
        if namdif
            || p.reyn != xf.reinf1
            || p.mach != xf.minf1
            || p.retyp != xf.retyp
            || p.matyp != xf.matyp
            || p.acrit != xf.acrit
            || p.xstrip[1] != xf.xstrip[1]
            || p.xstrip[2] != xf.xstrip[2]
        {
            crate::outln!(
                "\n              Current                         Save file\n          ------------------              ------------------\n name :   {}{}\n Re   :   {}{:20}{}\n Mach :   {}{:20}{}\n Retyp:   {}{:25}{}\n Matyp:   {}{:25}{}\n Ncrit:   {}{:20}{}\n xtr T:   {}{:20}{}\n xtr B:   {}{:20}{}",
                fixed(&xf.name, 48),
                fixed(&p.name, 48),
                f0(xf.reinf1, 12),
                "",
                f0(p.reyn, 12),
                f_fmt(xf.minf1, 12, 4),
                "",
                f_fmt(p.mach, 12, 4),
                i_fmt(xf.retyp as i64, 7),
                "",
                i_fmt(p.retyp as i64, 7),
                i_fmt(xf.matyp as i64, 7),
                "",
                i_fmt(p.matyp as i64, 7),
                f_fmt(xf.acrit, 12, 4),
                "",
                f_fmt(p.acrit, 12, 4),
                f_fmt(xf.xstrip[1], 12, 4),
                "",
                f_fmt(p.xstrip[1], 12, 4),
                f_fmt(xf.xstrip[2], 12, 4),
                "",
                f_fmt(p.xstrip[2], 12, 4)
            );
            crate::outln!();
            crate::outln!(" Current parameters different from old save file values.");
            let ok = s.input.askl("Set current parameters to old save file values ?^");
            if ok {
                let xf = &mut s.xf;
                xf.name = p.name.clone();
                xf.reinf1 = p.reyn;
                xf.minf1 = p.mach;
                xf.retyp = p.retyp;
                xf.matyp = p.matyp;
                xf.acrit = p.acrit;
                xf.xstrip[1] = p.xstrip[1];
                xf.xstrip[2] = p.xstrip[2];
            } else {
                crate::outln!();
                crate::outln!(" Old polar save file NOT available for appending");
                return;
            }
        }
        crate::outln!();
        crate::outln!(" Old polar save file read in ...");
        let na = p.na() as i64;
        let t = polwrit_cur(s, ip, true, 1, na, &p.code, p.vers);
        crate::out!("{}", t);
        s.extra_state.lpfile = true;
        crate::outln!();
        crate::outln!(" Old polar save file available for appending");
        return;
    }

    // ---- new file: write header
    set_nipol_a(s);
    let t = polwrit_cur(s, ip, true, 0, -1, "XFOIL", VERSION);
    let created = std::fs::OpenOptions::new().write(true).create_new(true).open(&fname).and_then(|mut f| {
        use std::io::Write;
        f.write_all(eol(&t).as_bytes())
    });
    if created.is_err() {
        // (NIPOL is left modified, as in the Fortran OPEN error branch)
        crate::outln!(" New polar save file OPEN error");
        return;
    }
    s.polars[ip - 1].pfname = fname;
    s.extra_state.nipol = s.extra_state.nipol0;
    s.extra_state.lpfile = true;
    crate::outln!();
    crate::outln!(" New polar save file available");
}

/// Payload helpers for unformatted records (REAL*8, INTEGER*4)
fn pr(b: &mut Vec<u8>, v: f64) {
    b.extend_from_slice(&v.to_le_bytes());
}
fn pi(b: &mut Vec<u8>, v: i32) {
    b.extend_from_slice(&v.to_le_bytes());
}
fn rd_f64(b: &[u8], k: usize) -> Option<f64> {
    b.get(8 * k..8 * k + 8).map(|s| f64::from_le_bytes(s.try_into().unwrap()))
}
fn rd_i32(b: &[u8], k: usize) -> Option<i32> {
    b.get(4 * k..4 * k + 4).map(|s| i32::from_le_bytes(s.try_into().unwrap()))
}

/// PLXINI: polar dump file set-up (Fortran unformatted file).
fn plxini(s: &mut Session, ip: usize) {
    let pfnamx = s.polars[ip - 1].pfnamx.trim().to_string();
    let npf = pfnamx.len();
    let prompt = if npf == 0 {
        "Enter  polar dump filename  OR  <return> for no file^".to_string()
    } else {
        crate::outln!(" Default polar dump filename:  {}", pfnamx);
        "Enter  new filename  OR  \"none\"  OR  <return> for default^".to_string()
    };
    let fname = fname64(&asks(s, &prompt));
    let mut nfn = fname.len();
    if "NONEnone".contains(substr(&fname, 1, 4).as_str()) {
        nfn = 0;
    }
    if nfn == 0 {
        s.extra_state.lpfilx = false;
        crate::outln!();
        crate::outln!(" Polar dump file will NOT be written");
        return;
    }
    s.extra_state.lpfilx = false;

    let data = std::fs::read(&fname).unwrap_or_default();
    if !data.is_empty() {
        let (recs, _ok) = ufread(&data);
        let hdr = (|| {
            let r1 = recs.first()?;
            if r1.len() < 32 {
                return None;
            }
            let namex = String::from_utf8_lossy(&r1[..32]).trim_matches(' ').to_string();
            let r2 = recs.get(1)?;
            let machx = rd_f64(r2, 0)?;
            let reynx = rd_f64(r2, 1)?;
            let acritx = rd_f64(r2, 2)?;
            let r3 = recs.get(2)?;
            let matypx = rd_i32(r3, 0)?;
            let retypx = rd_i32(r3, 1)?;
            let r4 = recs.get(3)?;
            rd_i32(r4, 3)?;
            Some((namex, machx, reynx * 1.0e6, acritx, matypx, retypx))
        })();
        match hdr {
            None => {
                crate::outln!(" Polar dump file READ error");
                return;
            }
            Some((namex, machx, reynx, acritx, matypx, retypx)) => {
                s.polars[ip - 1].pfnamx = fname.clone();
                let xf = &s.xf;
                let namdif = name48(&xf.name) != namex;
                if namdif || reynx != xf.reinf1 || machx != xf.minf1 || acritx != xf.acrit || retypx != xf.retyp || matypx != xf.matyp {
                    crate::outln!(
                        "\n              Dump file                       Current\n            ------------                    ------------\n name :   {}{}\n Re   :   {}{:20}{}\n Mach :   {}{:20}{}\n Retyp:   {}{:25}{}\n Matyp:   {}{:25}{}\n Ncrit:   {}{:20}{}",
                        fixed(&namex, 32),
                        fixed(&xf.name, 48),
                        f0(reynx, 12),
                        "",
                        f0(xf.reinf1, 12),
                        f_fmt(machx, 12, 4),
                        "",
                        f_fmt(xf.minf1, 12, 4),
                        i_fmt(retypx as i64, 7),
                        "",
                        i_fmt(xf.retyp as i64, 7),
                        i_fmt(matypx as i64, 7),
                        "",
                        i_fmt(xf.matyp as i64, 7),
                        f_fmt(acritx, 12, 4),
                        "",
                        f_fmt(xf.acrit, 12, 4)
                    );
                    crate::outln!();
                    crate::outln!(" Current parameters different from old dump file values.");
                    let ok = s.input.askl("Set current parameters to old dump file values ?^");
                    if ok {
                        let xf = &mut s.xf;
                        xf.name = namex;
                        xf.minf1 = machx;
                        xf.reinf1 = reynx;
                        xf.acrit = acritx;
                        xf.retyp = retypx;
                        xf.matyp = matypx;
                    } else {
                        crate::outln!();
                        crate::outln!(" Old polar dump file NOT available for appending");
                        return;
                    }
                }
                s.extra_state.lpfilx = true;
                crate::outln!();
                crate::outln!(" Old polar dump file available for appending");
                return;
            }
        }
    }

    // ---- new dump file: header records
    let xf = &s.xf;
    let mut buf = Vec::new();
    let mut r = Vec::new();
    r.extend_from_slice(fixed(&xf.name, 48).as_bytes());
    r.extend_from_slice(b"XFOIL   ");
    pr(&mut r, VERSION);
    ufrec(&mut buf, &r);
    let mut r = Vec::new();
    pr(&mut r, xf.minf1);
    pr(&mut r, xf.reinf1 / 1.0e6);
    pr(&mut r, xf.acrit);
    ufrec(&mut buf, &r);
    let mut r = Vec::new();
    pi(&mut r, xf.matyp);
    pi(&mut r, xf.retyp);
    ufrec(&mut buf, &r);
    let mut r = Vec::new();
    pi(&mut r, 0);
    pi(&mut r, 0);
    pi(&mut r, 0);
    pi(&mut r, xf.n as i32);
    ufrec(&mut buf, &r);
    let mut r = Vec::new();
    for i in 1..=xf.n {
        pr(&mut r, xf.x[i]);
        pr(&mut r, xf.y[i]);
    }
    ufrec(&mut buf, &r);
    if std::fs::write(&fname, &buf).is_err() {
        crate::outln!(" OPEN error on polar dump file {}", fixed(&fname, 48));
        return;
    }
    s.polars[ip - 1].pfnamx = fname;
    s.extra_state.lpfilx = true;
    crate::outln!();
    crate::outln!(" New polar dump file available");
}

/// PACC command (toggle polar accumulation).
pub fn pacc(s: &mut Session, iinput: &[i64]) {
    s.xf.lpacc = !s.xf.lpacc;
    if s.xf.lpacc {
        let np = npol(s);
        let ip = if !iinput.is_empty() { (iinput[0].max(0) as usize).min(np + 1) } else { np + 1 };
        let mut ipact;
        if ip > np {
            if np == NPX {
                crate::outln!();
                crate::outln!(" Number of polars is at array limit");
                crate::outln!(" New polar will not be stored");
                ipact = 0;
            } else {
                ipact = np + 1;
            }
        } else {
            ipact = ip;
        }
        ipact = plrset(s, ipact);
        s.ipact = ipact;
        if ipact == 0 {
            s.xf.lpacc = false;
            return;
        }
        plrini(s, ipact);
        plxini(s, ipact);
        crate::outln!();
        crate::outln!(" Polar accumulation enabled");
    } else {
        crate::outln!();
        crate::outln!(" Polar accumulation disabled");
        s.ipact = 0;
    }
}

/// PLRADD: add current point to active polar and its save file.
pub fn plradd(s: &mut Session) {
    let ip = s.ipact;
    if ip == 0 || ip > npol(s) {
        crate::outln!(" No active polar is declared. Point not stored.");
    } else if s.polars[ip - 1].na() >= NAX {
        crate::outln!(" Polar storage arrays full. Point not stored");
    } else {
        if s.xf.lflap {
            mhinge(s);
        }
        let xf = &s.xf;
        let (cdtot, cdv, re) = if xf.lvisc { (xf.cd, xf.cd, xf.reinf) } else { (0.0, 0.0, 0.0) };
        let mut row = [0.0f64; IPTOT + 1];
        row[IAL] = xf.adeg;
        row[ICL] = xf.cl;
        row[ICD] = cdtot;
        row[ICM] = xf.cm;
        row[ICP] = xf.cdp;
        row[ICV] = cdv;
        row[IMA] = xf.minf;
        row[IRE] = re;
        row[INC] = xf.acrit;
        let mut sd = [[0.0f64; 3]; JPTOT + 1];
        for is in 1..=2 {
            let xoct = if xf.lvisc { xf.xoctr[is] } else { 0.0 };
            sd[JTP][is] = xf.xstrip[is];
            sd[JTN][is] = xoct;
        }
        row[ICH] = if xf.lflap { xf.hmom } else { 0.0 };
        row[IMC] = xf.cpmn;
        let p = &mut s.polars[ip - 1];
        p.cpol.push(row);
        p.cpolsd.push(sd);
        crate::outln!("\n Point added to stored polar{}", i_fmt(ip as i64, 3));
    }
    if s.extra_state.lpfile && ip >= 1 && ip <= npol(s) {
        set_nipol_a(s);
        let na = s.polars[ip - 1].na() as i64;
        let t = polwrit_cur(s, ip, false, na, na, "XFOIL", VERSION);
        let f = s.polars[ip - 1].pfname.clone();
        if let Ok(mut fh) = std::fs::OpenOptions::new().append(true).create(true).open(&f) {
            use std::io::Write;
            let _ = fh.write_all(eol(&t).as_bytes());
        }
        s.extra_state.nipol = s.extra_state.nipol0;
        crate::outln!(" Point written to save file  {}", fixed(&f, 48));
    } else {
        crate::outln!(" Save file unspecified or not available");
    }
}

/// PLXADD: append current point to the unformatted polar dump file.
pub fn plxadd(s: &mut Session) {
    if !s.extra_state.lpfilx || s.ipact == 0 || s.ipact > npol(s) {
        crate::outln!(" Dump file unspecified or not available");
        return;
    }
    let fname = s.polars[s.ipact - 1].pfnamx.clone();
    let xf = &mut s.xf;
    let beta = (1.0 - xf.minf * xf.minf).sqrt();
    let bfac = 0.5 * xf.minf * xf.minf / (1.0 + beta);
    let mut buf = Vec::new();
    let (cdtot, xt1, xt2) = if xf.lvisc { (xf.cd, xf.xoctr[1], xf.xoctr[2]) } else { (0.0, 0.0, 0.0) };
    let mut r = Vec::new();
    for v in [xf.alfa / DTOR, xf.cl, cdtot, 0.0, xf.cm, xt1, xt2] {
        pr(&mut r, v);
    }
    ufrec(&mut buf, &r);
    let mut nside = [0usize; 3];
    nside[1] = xf.iblte[1] + (xf.nbl[2].wrapping_sub(xf.iblte[2]));
    nside[2] = xf.nbl[2];
    nside[1] = nside[1].max(2);
    nside[2] = nside[2].max(2);
    let mut r = Vec::new();
    pi(&mut r, nside[1] as i32);
    pi(&mut r, nside[2] as i32);
    pi(&mut r, xf.iblte[1] as i32);
    pi(&mut r, xf.iblte[2] as i32);
    ufrec(&mut buf, &r);
    let que = 0.5 * xf.qinf * xf.qinf;
    let nmax = nside[1].max(nside[2]) + 1;
    let mut xx = [vec![0.0; nmax + 1], vec![0.0; nmax + 1], vec![0.0; nmax + 1]];
    let mut cp = xx.clone();
    let mut cf = xx.clone();
    let n = xf.n;
    xx[1][1] = crate::spline::seval(xf.sst, &xf.x, &xf.xp, &xf.s, n);
    cp[1][1] = 1.0 / (beta + bfac);
    cf[1][1] = 0.0;
    xf.thet[1][1] = 0.5 * (xf.thet[1][2] + xf.thet[2][2]);
    xf.dstr[1][1] = 0.5 * (xf.dstr[1][2] + xf.dstr[2][2]);
    xf.ctau[1][1] = 0.0;
    xx[2][1] = xx[1][1];
    cp[2][1] = cp[1][1];
    cf[2][1] = cf[1][1];
    xf.thet[2][1] = xf.thet[1][1];
    xf.dstr[2][1] = xf.dstr[1][1];
    xf.ctau[2][1] = xf.ctau[1][1];
    for is in 1..=2 {
        for ibl in 2..=nside[is] {
            let i = xf.ipan[is][ibl];
            xx[is][ibl] = xf.x[i];
            cp[is][ibl] = xf.cpv[i];
            cf[is][ibl] = xf.tau[is][ibl] / que;
        }
    }
    for is in 1..=2 {
        let mut r = Vec::new();
        for ibl in 1..=nside[is] {
            pr(&mut r, xx[is][ibl]);
            pr(&mut r, cp[is][ibl]);
            pr(&mut r, xf.thet[is][ibl]);
            pr(&mut r, xf.dstr[is][ibl]);
            pr(&mut r, cf[is][ibl]);
            pr(&mut r, xf.ctau[is][ibl]);
        }
        ufrec(&mut buf, &r);
    }
    if let Ok(mut fh) = std::fs::OpenOptions::new().append(true).create(true).open(&fname) {
        use std::io::Write;
        let _ = fh.write_all(&buf);
    }
    crate::outln!(" Point written to dump file {}", fixed(&fname, 48));
}

/// PLRCOP-style removal of polar ip (1-based), reporting moves.
fn delete_polar(s: &mut Session, ip: usize, track_ipact: bool) {
    let np = npol(s);
    for jp in ip + 1..=np {
        crate::outln!(" Polar{}  moved into polar{}", i_fmt(jp as i64, 3), i_fmt(jp as i64 - 1, 3));
        if track_ipact && s.ipact == jp {
            s.ipact = jp - 1;
        }
    }
    s.polars.remove(ip - 1);
}

// ------------------------------------------------------------------------
// OPER polar commands

fn pget(s: &mut Session, arg: &str) {
    if npol(s) >= NPX {
        crate::outln!();
        crate::outln!(" Number of polars is at array limit");
        crate::outln!(" Delete with PDEL if necessary");
        return;
    }
    let fname = if arg.trim().is_empty() { fname64(&asks(s, "Enter polar filename^")) } else { fname64(arg) };
    let mut p = Polar::default();
    let err = match std::fs::read(&fname) {
        Ok(b) => polread(&String::from_utf8_lossy(&b), &mut p),
        Err(_) => true,
    };
    if err {
        crate::outln!(" Polar file READ error");
        return;
    }
    p.nxy = 0;
    p.name = name48(&p.name);
    p.pfname = fname;
    s.polars.push(p);
    let ip = npol(s);
    let na = s.polars[ip - 1].na() as i64;
    let (code, vers) = (s.polars[ip - 1].code.clone(), s.polars[ip - 1].vers);
    let t = polwrit_cur(s, ip, true, 1, na, &code, vers);
    crate::out!("{}", t);
    crate::outln!("\n Stored as  Polar{}", i_fmt(ip as i64, 4));
}

fn pwrt(s: &mut Session, iinput: &[i64]) {
    let mut ninput = iinput.len();
    let mut ip: i64 = 0;
    let (ip1, ip2);
    loop {
        let np = npol(s) as i64;
        if np == 1 {
            ip = 1;
        } else if ninput == 0 {
            plrsum(s, 1, np as usize);
            ip = s.input.aski("Enter index of polar to write (0=all, -1=abort)^", ip);
            if ip == -1 || s.input.eof {
                return;
            }
        } else {
            ip = iinput[0];
        }
        if ip == 0 {
            ip1 = 1;
            ip2 = np;
            break;
        } else if ip >= 1 && ip <= np {
            ip1 = ip;
            ip2 = ip;
            break;
        } else {
            ninput = 0;
        }
    }
    for ip in ip1..=ip2 {
        let ipu = ip as usize;
        plrsum(s, ipu, ipu);
        let pf = s.polars[ipu - 1].pfname.trim().to_string();
        let prompt = if pf.is_empty() { "Enter polar output filename^".to_string() } else { format!("Enter polar output filename [{}]^", pf) };
        let mut fname = fname64(&asks(s, &prompt));
        if !pf.is_empty() && fname.is_empty() {
            fname = pf.clone();
        }
        set_nipol_b(s);
        let t = polwrit_cur(s, ipu, true, 1, s.polars[ipu - 1].na() as i64, "XFOIL", VERSION);
        let ok = polwrit_file(s, &fname, &t);
        if !ok {
            crate::outln!(" Polar{}  not written", i_fmt(ip, 3));
        } else {
            s.polars[ipu - 1].pfname = fname;
        }
    }
}

/// POLWRIT file handling with LQUERY=.TRUE.: returns false on ERROR.
fn polwrit_file(s: &mut Session, fname: &str, text: &str) -> bool {
    let f = if fname.is_empty() { "fort.19".to_string() } else { fname.to_string() };
    if !fname.is_empty() && std::path::Path::new(&f).exists() {
        crate::outln!();
        crate::outln!(" Output file exists.  Overwrite?  Y");
        let ans = readln(s);
        if ans.starts_with(['N', 'n']) {
            crate::outln!(" Polar file not saved");
            return true;
        }
    }
    std::fs::write(&f, eol(text)).is_ok()
}

fn plis(s: &mut Session, iinput: &[i64]) {
    let np = npol(s);
    if np == 0 {
        crate::outln!();
        crate::outln!(" No polars are stored");
        return;
    }
    let (ip1, ip2) = if iinput.is_empty() {
        (1, np)
    } else {
        let ip = iinput[0];
        if ip == 0 {
            (1, np)
        } else if ip >= 1 && ip as usize <= np {
            (ip as usize, ip as usize)
        } else {
            crate::outln!();
            crate::outln!(" Specified stored polar does not exist");
            return;
        }
    };
    set_nipol_b(s);
    for ip in ip1..=ip2 {
        crate::outln!("\n ==============================================================\n Polar{}", i_fmt(ip as i64, 3));
        let na = s.polars[ip - 1].na() as i64;
        let t = polwrit_cur(s, ip, true, 1, na, "XFOIL", VERSION);
        crate::out!("{}", t);
    }
    s.extra_state.nipol = s.extra_state.nipol0;
}

fn pdel(s: &mut Session, iinput: &[i64]) {
    if npol(s) == 0 {
        crate::outln!(" No polars are stored");
        return;
    }
    let mut ip: i64 = if !iinput.is_empty() { iinput[0] } else { npol(s) as i64 + 1 };
    loop {
        let np = npol(s) as i64;
        if ip == 0 {
            s.polars.clear();
            s.ipact = 0;
            s.xf.lpacc = false;
            break;
        } else if ip == -1 {
            return;
        } else if ip < -1 || ip > np {
            plrsum(s, 1, np as usize);
            ip = s.input.aski("Specify polar to delete (0 = all, -1 = abort)^", ip);
            if s.input.eof {
                return;
            }
            continue;
        } else {
            let ipu = ip as usize;
            if s.ipact == ipu {
                crate::outln!(" Active polar deleted.  Accumulation turned off");
                s.ipact = 0;
                s.xf.lpacc = false;
            }
            delete_polar(s, ipu, true);
            break;
        }
    }
    if s.ipact > 0 {
        crate::outln!(" Polar{}  is now active", i_fmt(s.ipact as i64, 3));
    }
}

fn plrsrt(p: &mut Polar, idsort: usize) {
    let n = p.na();
    let mut a = vec![0.0; n + 1];
    for i in 1..=n {
        a[i] = p.cpol[i - 1][idsort];
    }
    let indx = hsort(n, &a);
    let old = p.cpol.clone();
    let oldsd = p.cpolsd.clone();
    for i in 1..=n {
        p.cpol[i - 1] = old[indx[i] - 1];
        p.cpolsd[i - 1] = oldsd[indx[i] - 1];
    }
}

fn psor(s: &mut Session, iinput: &[i64]) {
    if npol(s) == 0 {
        crate::outln!(" No polars are stored");
        return;
    }
    let mut ip: i64 = if !iinput.is_empty() { iinput[0] } else { npol(s) as i64 + 1 };
    loop {
        let np = npol(s) as i64;
        if ip == -1 {
            return;
        } else if ip < -1 || ip > np {
            plrsum(s, 1, np as usize);
            ip = s.input.aski("Specify polar to sort (0 = all, -1 = abort)^", ip);
            if s.input.eof {
                return;
            }
            continue;
        } else {
            let (ip1, ip2) = if ip == 0 { (1, np as usize) } else { (ip as usize, ip as usize) };
            for jp in ip1..=ip2 {
                plrsrt(&mut s.polars[jp - 1], IAL);
            }
            return;
        }
    }
}

/// common "select polar" logic of ASET/PREM/PNAM. Returns None to abort.
fn select_polar(s: &mut Session, iinput: &[i64], prompt: &str) -> Option<usize> {
    let mut ninput = iinput.len();
    let mut ip: i64 = 0;
    loop {
        let np = npol(s) as i64;
        if ninput == 0 {
            if np == 1 {
                ip = 1;
            } else {
                plrsum(s, 1, np as usize);
                ip = s.input.aski(prompt, ip);
                if s.input.eof {
                    return None;
                }
            }
        } else {
            ip = iinput[0];
        }
        if ip == 0 {
            return None;
        } else if ip < 1 || ip > np {
            crate::outln!();
            crate::outln!(" Specified polar airfoil does not exist");
            ninput = 0;
        } else {
            return Some(ip as usize);
        }
    }
}

fn aset(s: &mut Session, iinput: &[i64]) {
    if npol(s) == 0 {
        crate::outln!();
        crate::outln!(" No polar airfoils are stored");
        return;
    }
    let ip = match select_polar(s, iinput, "Enter index of polar airfoil to set^") {
        Some(ip) => ip,
        None => return,
    };
    crate::outln!();
    crate::outln!(" Current airfoil will be overwritten.  Proceed?  Y");
    let ans = readln(s);
    if ans.starts_with(['N', 'n']) {
        crate::outln!(" No action taken");
    } else {
        apcopy(s, ip);
    }
}

fn prem(s: &mut Session, r: &[f64], iinput: &[i64]) {
    if npol(s) == 0 {
        crate::outln!();
        crate::outln!(" No polars are stored");
        return;
    }
    let ip = match select_polar(s, iinput, "Enter index of polar to modify^") {
        Some(ip) => ip,
        None => return,
    };
    let alphas: Vec<f64> = if r.len() >= 2 {
        r[1..].to_vec()
    } else {
        set_nipol_b(s);
        crate::outln!("\n ==============================================================\n Polar{}", i_fmt(ip as i64, 3));
        let na = s.polars[ip - 1].na() as i64;
        let t = polwrit_cur(s, ip, true, 1, na, "XFOIL", VERSION);
        crate::out!("{}", t);
        loop {
            crate::out!("\n Enter alpha(s) of points to be removed:  ");
            crate::flush_if_tty();
            let l = readln(s);
            if s.input.eof {
                return;
            }
            if let Some(v) = crate::geom::getflt(&l, 19) {
                break v;
            }
        }
    };
    for a in alphas {
        let p = &mut s.polars[ip - 1];
        let mut found = None;
        for ia in 0..p.na() {
            if (p.cpol[ia][IAL] - a).abs() < 0.0005 {
                found = Some(ia);
                break;
            }
        }
        if let Some(ia) = found {
            p.cpol.remove(ia);
            p.cpolsd.remove(ia);
            if p.na() == 0 {
                delete_polar(s, ip, true);
                if s.ipact > 0 {
                    crate::outln!(" Polar{}  is now active", i_fmt(s.ipact as i64, 3));
                }
                return;
            }
        }
    }
}

fn pnam(s: &mut Session, iinput: &[i64]) {
    if npol(s) == 0 {
        crate::outln!();
        crate::outln!(" No polars are stored");
        return;
    }
    let ip = match select_polar(s, iinput, "Enter index of polar to modify^") {
        Some(ip) => ip,
        None => return,
    };
    set_nipol_b(s);
    crate::outln!("\n ==============================================================\n Polar{}", i_fmt(ip as i64, 3));
    let t = polwrit_cur(s, ip, true, 0, -1, "XFOIL", VERSION);
    crate::out!("{}", t);
    s.extra_state.nipol = s.extra_state.nipol0;
    crate::out!("\n Enter new airfoil name of polar:  ");
    crate::flush_if_tty();
    let l = readln(s);
    s.polars[ip - 1].name = name48(&l);
}

fn rget(s: &mut Session, arg: &str) {
    if s.extra_state.refs.len() >= NPX {
        crate::outln!();
        crate::outln!(" Number of reference polars is at array limit");
        crate::outln!(" Delete with RDEL if necessary");
        return;
    }
    let fname = if arg.trim().is_empty() { fname64(&asks(s, "Enter reference polar filename^")) } else { fname64(arg) };
    let r = match std::fs::read(&fname) {
        Ok(b) => polref(&String::from_utf8_lossy(&b)),
        Err(_) => None,
    };
    let mut r = match r {
        Some(r) => r,
        None => {
            crate::outln!(" File OPEN error");
            return;
        }
    };
    r.name = name48(&r.name);
    if r.name.is_empty() {
        let l = asks(s, "Enter label for reference polar^");
        r.name = name48(&l);
    } else {
        crate::outln!();
        crate::outln!(" {}", fixed(&r.name, 48));
    }
    s.extra_state.refs.push(r);
}

fn rdel(s: &mut Session, iinput: &[i64]) {
    if s.extra_state.refs.is_empty() {
        crate::outln!(" No reference polars are stored");
        return;
    }
    let mut ir: i64 = if !iinput.is_empty() { iinput[0] } else { s.extra_state.refs.len() as i64 + 1 };
    loop {
        let nr = s.extra_state.refs.len() as i64;
        if ir == 0 {
            s.extra_state.refs.clear();
            return;
        } else if ir == -1 {
            return;
        } else if ir < -1 || ir > nr {
            prfsum(s);
            ir = s.input.aski("Specify ref. polar to delete (0 = all, -1 = abort)^", ir);
            if s.input.eof {
                return;
            }
        } else {
            for jr in ir + 1..=nr {
                crate::outln!(" Polar{}  moved into polar{}", i_fmt(jr, 3), i_fmt(jr - 1, 3));
            }
            s.extra_state.refs.remove(ir as usize - 1);
            return;
        }
    }
}

// ------------------------------------------------------------------------
// Top level

/// Top-level commands not handled in cli.rs. Return true if handled.
pub fn top_level(s: &mut Session, com: &str, arg: &str) -> bool {
    match com {
        "INTE" => inte(s, false),
        "INTX" => inte(s, true),
        "HALF" => {
            let xf = &mut s.xf;
            let nb = xf.nb;
            let nnew = half(&mut xf.xb, &mut xf.yb, &xf.sb.clone(), nb);
            xf.nb = nnew;
            buffer_geopar(s);
        }
        "BEND" => {
            if s.xf.n == 0 {
                crate::outln!();
                crate::outln!("  No airfoil available");
            } else {
                bendump(&s.xf.x, &s.xf.y, s.xf.n);
            }
        }
        "BENP" => {
            if s.xf.n == 0 {
                crate::outln!();
                crate::outln!("  No airfoil available");
            } else {
                let t = vec![1.0; s.xf.n + 1];
                bendump2(&s.xf.x, &s.xf.y, &t, s.xf.n);
            }
        }
        "MSAV" => msave(s, arg),
        "WDEF" => {
            let fname = if arg.is_empty() || arg.starts_with(' ') { "xfoil.def".to_string() } else { fname64(arg) };
            if std::path::Path::new(&fname).exists() {
                crate::outln!("\n  File  {}  exists.  Overwrite?  Y", fname);
                let ans = readln(s);
                if ans.starts_with(['N', 'n']) {
                    // (the Fortran falls through and writes the file anyway)
                    crate::outln!();
                    crate::outln!(" No action taken");
                }
            }
            let t = wrtdef(s);
            let _ = std::fs::write(&fname, eol(&t));
            crate::outln!("\n  File  {}  written", fname);
        }
        "RDEF" => {
            let fname = if arg.is_empty() || arg.starts_with(' ') { "xfoil.def".to_string() } else { fname64(arg) };
            getdef(s, &fixed(&fname, 64), false);
        }
        "NAME" => {
            if arg.trim().is_empty() {
                let nm = s.xf.name.clone();
                s.xf.name = nammod(s, &nm, 0, -1);
            } else {
                s.xf.name = name48(arg);
            }
        }
        "NINC" => {
            let nm = s.xf.name.clone();
            s.xf.name = nammod(s, &nm, 1, 1);
        }
        "PLOP" => oplset(s),
        _ => return false,
    }
    true
}

/// SAVE file text. iftyp: -1 integer, 0 plain, 1 labeled, 2 ISES.
pub fn save_text(xf: &crate::xfoil::XFoil, iftyp: i32) -> String {
    let delim = match xf.kdelim {
        0 => ' ',
        1 => ',',
        2 => '\t',
        _ => {
            crate::outln!(" ? Illegal delimiter.  Using blank.");
            ' '
        }
    };
    let mut out = String::new();
    if iftyp >= 1 {
        out += &name48(&xf.name);
        out.push('\n');
    }
    if iftyp >= 2 {
        let isp: String = xf.ispars.chars().take(80).collect();
        out += isp.trim_end_matches(' ');
        out.push('\n');
    }
    let n = xf.n;
    let idx: Vec<usize> = if xf.lclock { (1..=n).rev().collect() } else { (1..=n).collect() };
    for i in idx {
        if iftyp == -1 {
            let r = |v: f64| (v + 0.5f64.copysign(v)).trunc() as i64;
            out += &format!(" {}{}\n", i_fmt(r(xf.x[i]), 12), i_fmt(r(xf.y[i]), 12));
        } else if xf.kdelim == 0 {
            out += &format!(" {}{}\n", g_fmt(xf.x[i], 15, 7), g_fmt(xf.y[i], 15, 7));
        } else {
            let line = format!(" {}{}{}", f_fmt(xf.x[i], 10, 6), delim, f_fmt(xf.y[i], 10, 6));
            out += &bstrip(&line);
            out.push('\n');
        }
    }
    out
}

/// SCALC, SEGSPL, GEOPAR on buffer airfoil (prints thickness/camber).
fn buffer_geopar(s: &mut Session) {
    let xf = &mut s.xf;
    let nb = xf.nb;
    crate::spline::scalc(&xf.xb, &xf.yb, &mut xf.sb, nb);
    crate::spline::segspl(&xf.xb, &mut xf.xbp, &xf.sb, nb);
    crate::spline::segspl(&xf.yb, &mut xf.ybp, &xf.sb, nb);
    let g = crate::geom::geopar(&xf.xb, &xf.xbp, &xf.yb, &xf.ybp, &xf.sb, nb, true);
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

/// HALF: returns new number of points.
pub fn half(x: &mut [f64], y: &mut [f64], s: &[f64], n: usize) -> usize {
    let mut k = 1;
    let mut inext = 3;
    for i in 2..n {
        if s[i] == s[i + 1] {
            k += 1;
            x[k] = x[i];
            y[k] = y[i];
            k += 1;
            x[k] = x[i + 1];
            y[k] = y[i + 1];
            inext = i + 3;
        }
        if i == inext {
            k += 1;
            x[k] = x[i];
            y[k] = y[i];
            inext = i + 2;
        }
    }
    k += 1;
    x[k] = x[n];
    y[k] = y[n];
    k
}

/// INTE / INTX
fn inte(s: &mut Session, lx: bool) {
    let np = npol(s);
    crate::outln!("\n  F  disk file\n  C  current airfoil  {}", fixed(&s.xf.name, 48));
    for ip in 1..=np {
        if s.polars[ip - 1].nxy > 0 {
            crate::outln!(" {}  polar airfoil    {}", i_fmt(ip as i64, 2), fixed(&s.polars[ip - 1].name, 48));
        }
    }
    let promptn = match np {
        0 => "\" ( F C ):  ",
        1 => "\" ( F C 1 ):  ",
        2 => "\" ( F C 1 2 ):  ",
        _ => "\" ( F C 1 2.. ):  ",
    };
    let invalid = |lx: bool| {
        crate::outln!();
        if lx {
            crate::outln!(" Invalid response.  No action taken.");
        } else {
            crate::outln!(" Invalid response");
        }
    };
    let mut xint: [Vec<f64>; 2] = [vec![0.0; IBX + 1], vec![0.0; IBX + 1]];
    let mut yint = xint.clone();
    let mut xpint = xint.clone();
    let mut ypint = xint.clone();
    let mut sint = xint.clone();
    let mut nint = [0usize; 2];
    let mut sleint = [0.0f64; 2];
    let mut nameint = [String::new(), String::new()];
    for k in 0..2 {
        loop {
            crate::out!("\n  Select source of airfoil \"{}{}", k, promptn);
            crate::flush_if_tty();
            let line = readln(s);
            if s.input.eof && line.is_empty() {
                invalid(lx);
                return;
            }
            let cair: String = line.chars().take(2).collect();
            let c1 = cair.chars().next().unwrap_or(' ');
            if cair.trim().is_empty() {
                invalid(lx);
                return;
            } else if c1 == 'F' || c1 == 'f' {
                let fname = asks(s, "Enter filename^");
                let r = crate::geom::aread(fname.trim(), IBX);
                if r.itype == 0 {
                    return;
                }
                nint[k] = r.n;
                for i in 1..=r.n {
                    xint[k][i] = r.x[i];
                    yint[k][i] = r.y[i];
                }
                nameint[k] = name48(&r.name);
                break;
            } else if c1 == 'C' || c1 == 'c' {
                if s.xf.n <= 1 {
                    crate::outln!(" No current airfoil available");
                    continue;
                }
                nint[k] = s.xf.n;
                for i in 1..=s.xf.n {
                    xint[k][i] = s.xf.x[i];
                    yint[k][i] = s.xf.y[i];
                }
                nameint[k] = s.xf.name.clone();
                break;
            } else {
                let ip = match list_int1(&cair) {
                    Some(v) => v,
                    None => {
                        invalid(lx);
                        return;
                    }
                };
                if ip < 1 || ip as usize > np || s.polars[ip as usize - 1].nxy == 0 {
                    invalid(lx);
                    return;
                }
                let p = &s.polars[ip as usize - 1];
                nint[k] = p.nxy;
                // (the Fortran copies N points of the stored airfoil)
                for i in 1..=s.xf.n.min(p.xyx.len() - 1) {
                    xint[k][i] = p.xyx[i];
                    yint[k][i] = p.xyy[i];
                }
                nameint[k] = p.name.clone();
                break;
            }
        }
        let n = nint[k];
        crate::spline::scalc(&xint[k], &yint[k], &mut sint[k], n);
        crate::spline::segspld(&xint[k], &mut xpint[k], &sint[k], n, -999.0, -999.0);
        crate::spline::segspld(&yint[k], &mut ypint[k], &sint[k], n, -999.0, -999.0);
        sleint[k] = crate::geom::lefind(&xint[k], &xpint[k], &yint[k], &ypint[k], &sint[k], n);
    }
    crate::outln!();
    crate::outln!(" airfoil \"0\":  {}", fixed(&nameint[0], 48));
    crate::outln!(" airfoil \"1\":  {}", fixed(&nameint[1], 48));
    let frac = s.input.askr("Specify interpolating fraction  0...1^", 0.5);
    {
        let xf = &mut s.xf;
        let (n0, n1) = (nint[0], nint[1]);
        let (x0, xp0, y0, yp0, s0, sle0) = (&xint[0], &xpint[0], &yint[0], &ypint[0], &sint[0], sleint[0]);
        let (x1, xp1, y1, yp1, s1, sle1) = (&xint[1], &xpint[1], &yint[1], &ypint[1], &sint[1], sleint[1]);
        let n = n0;
        let f0 = 1.0 - frac;
        let f1 = frac;
        use crate::spline::{seval, sinvrt};
        if !lx {
            // INTER
            let tops0 = s0[1] - sle0;
            let tops1 = s1[1] - sle1;
            let bots0 = s0[n0] - sle0;
            let bots1 = s1[n1] - sle1;
            let mut sn = 0.0;
            for i in 1..=n {
                if s0[i] < sle0 {
                    sn = (s0[i] - sle0) / tops0;
                }
                if s0[i] >= sle0 {
                    sn = (s0[i] - sle0) / bots0;
                }
                let st0 = s0[i];
                let st1 = if st0 < sle0 { sle1 + tops1 * sn } else { sle1 + bots1 * sn };
                let xt0 = seval(st0, x0, xp0, s0, n0);
                let yt0 = seval(st0, y0, yp0, s0, n0);
                let xt1 = seval(st1, x1, xp1, s1, n1);
                let yt1 = seval(st1, y1, yp1, s1, n1);
                xf.xb[i] = f0 * xt0 + f1 * xt1;
                xf.yb[i] = f0 * yt0 + f1 * yt1;
            }
        } else {
            // INTERX
            let xle0 = seval(sle0, x0, xp0, s0, n0);
            let xle1 = seval(sle1, x1, xp1, s1, n1);
            let mut xn = 0.0;
            for i in 1..=n {
                if s0[i] < sle0 {
                    xn = (x0[i] - xle0) / (x0[1] - xle0);
                }
                if s0[i] >= sle0 {
                    xn = (x0[i] - xle0) / (x0[n0] - xle0);
                }
                let mut xt0 = x0[i];
                let mut st0 = s0[i];
                let (mut xt1, mut st1);
                if st0 < sle0 {
                    xt1 = xle1 + (x1[1] - xle1) * xn;
                    st1 = sle1 + (s1[1] - sle1) * xn;
                } else {
                    xt1 = xle1 + (x1[n1] - xle1) * xn;
                    st1 = sle1 + (s1[n1] - sle1) * xn;
                }
                sinvrt(&mut st0, xt0, x0, xp0, s0, n0);
                sinvrt(&mut st1, xt1, x1, xp1, s1, n1);
                xt0 = seval(st0, x0, xp0, s0, n0);
                let yt0 = seval(st0, y0, yp0, s0, n0);
                xt1 = seval(st1, x1, xp1, s1, n1);
                let yt1 = seval(st1, y1, yp1, s1, n1);
                xf.xb[i] = f0 * xt0 + f1 * xt1;
                xf.yb[i] = f0 * yt0 + f1 * yt1;
            }
        }
        xf.nb = n;
    }
    buffer_geopar(s);
    let nm = asks(s, "Enter new airfoil name^");
    s.xf.name = name48(&nm);
    crate::outln!();
    crate::outln!(" Result has been placed in buffer airfoil");
    crate::outln!(" Execute PCOP or PANE to set new current airfoil");
}

/// IJSECT: section properties.
#[allow(clippy::type_complexity)]
fn ijsect(n: usize, x: &[f64], y: &[f64]) -> [f64; 16] {
    let mut xmin = x[1];
    let mut xmax = x[1];
    let mut ymin = y[1];
    let mut ymax = y[1];
    let mut dx = x[1] - x[n];
    let mut dy = y[1] - y[n];
    let mut ds = (dx * dx + dy * dy).sqrt();
    let mut xavg = 0.5 * (x[1] + x[n]);
    let mut yavg = 0.5 * (y[1] + y[n]);
    let mut x_dy = dy * xavg;
    let mut xx_dy = dy * xavg.powi(2);
    let mut xxx_dy = dy * xavg.powi(3);
    let mut x_ds = ds * xavg;
    let mut xx_ds = ds * xavg.powi(2);
    let mut y_dx = dx * yavg;
    let mut yy_dx = dx * yavg.powi(2);
    let mut yyy_dx = dx * yavg.powi(3);
    let mut y_ds = ds * yavg;
    let mut yy_ds = ds * yavg.powi(2);
    let mut c_ds = ds;
    for i in 2..=n {
        dx = x[i] - x[i - 1];
        dy = y[i] - y[i - 1];
        ds = (dx * dx + dy * dy).sqrt();
        xavg = 0.5 * (x[i] + x[i - 1]);
        yavg = 0.5 * (y[i] + y[i - 1]);
        x_dy += dy * xavg;
        xx_dy += dy * xavg.powi(2);
        xxx_dy += dy * xavg.powi(3);
        x_ds += ds * xavg;
        xx_ds += ds * xavg.powi(2);
        y_dx += dx * yavg;
        yy_dx += dx * yavg.powi(2);
        yyy_dx += dx * yavg.powi(3);
        y_ds += ds * yavg;
        yy_ds += ds * yavg.powi(2);
        c_ds += ds;
        xmin = xmin.min(x[i]);
        xmax = xmax.max(x[i]);
        ymin = ymin.min(y[i]);
        ymax = ymax.max(y[i]);
    }
    let area = -y_dx;
    let slen = c_ds;
    let mut o = [0.0f64; 16];
    o[0] = area;
    o[1] = slen;
    o[2] = xmin;
    o[3] = xmax;
    o[4] = ymin;
    o[5] = ymax;
    if area == 0.0 {
        return o;
    }
    let xc = xx_dy / (2.0 * x_dy);
    let xct = x_ds / c_ds;
    let aiyy = xxx_dy / 3.0 - xx_dy * xc + x_dy * xc.powi(2);
    let aiyyt = xx_ds - x_ds * xct * 2.0 + c_ds * xct.powi(2);
    let yc = yy_dx / (2.0 * y_dx);
    let yct = y_ds / c_ds;
    let aixx = -yyy_dx / 3.0 + yy_dx * yc - y_dx * yc.powi(2);
    let aixxt = yy_ds - y_ds * yct * 2.0 + c_ds * yct.powi(2);

    // IMID: first i with x(i+1) >= x(i), else n
    let mut imid = n;
    for i in 1..n {
        if x[i + 1] >= x[i] {
            imid = i;
            break;
        }
    }
    let mut aj = 0.0;
    for i in 2..=imid {
        let xavg = 0.5 * (x[i] + x[i - 1]);
        let yavg = 0.5 * (y[i] + y[i - 1]);
        let dx = x[i - 1] - x[i];
        let mut yopp = 0.0;
        if xavg > x[n] {
            yopp = y[n];
        } else if xavg <= x[imid] {
            yopp = y[imid];
        } else {
            let mut j = n;
            while j >= imid && j >= 1 {
                if xavg > x[j - 1] && xavg <= x[j] {
                    let frac = (xavg - x[j - 1]) / (x[j] - x[j - 1]);
                    yopp = y[j - 1] + (y[j] - y[j - 1]) * frac;
                    break;
                }
                j -= 1;
            }
        }
        aj += (yavg - yopp).abs().powi(3) * dx / 3.0;
    }
    let ajt = 4.0 * area * area / slen;
    o[6] = xc;
    o[7] = yc;
    o[8] = xct;
    o[9] = yct;
    o[10] = aixx;
    o[11] = aiyy;
    o[12] = aixxt;
    o[13] = aiyyt;
    o[14] = aj;
    o[15] = ajt;
    o
}

fn g14(v: f64) -> String {
    g_fmt(v, 14, 6)
}

fn bendump(x: &[f64], y: &[f64], n: usize) {
    let o = ijsect(n, x, y);
    let (area, slen, xmin, xmax, ymin, ymax) = (o[0], o[1], o[2], o[3], o[4], o[5]);
    let (xc, yc, xct, yct, aixx, aiyy, aixxt, aiyyt, aj, ajt) = (o[6], o[7], o[8], o[9], o[10], o[11], o[12], o[13], o[14], o[15]);
    let w = |a: &str, v: f64| crate::outln!(" {}{}", a, g14(v));
    crate::outln!();
    w("Area =", area);
    w("Slen =", slen);
    crate::outln!();
    crate::outln!(" X-bending parameters(solid):");
    w("        Xc =", xc);
    w("  max X-Xc =", xmax - xc);
    w("  min X-Xc =", xmin - xc);
    w("       Iyy =", aiyy);
    let xbar = (xmax - xc).abs().max((xmin - xc).abs());
    w(" Iyy/(X-Xc)=", aiyy / xbar);
    crate::outln!();
    crate::outln!(" Y-bending parameters(solid):");
    w("        Yc =", yc);
    w("  max Y-Yc =", ymax - yc);
    w("  min Y-Yc =", ymin - yc);
    w("       Ixx =", aixx);
    let ybar = (ymax - yc).abs().max((ymin - yc).abs());
    w(" Ixx/(Y-Yc)=", aixx / ybar);
    crate::outln!();
    w("       J   =", aj);
    crate::outln!();
    crate::outln!();
    crate::outln!(" X-bending parameters(skin):");
    w("         Xc =", xct);
    w("   max X-Xc =", xmax - xct);
    w("   min X-Xc =", xmin - xct);
    w("      Iyy/t =", aiyyt);
    let xbart = (xmax - xct).abs().max((xmin - xct).abs());
    w(" Iyy/t(X-Xc)=", aiyyt / xbart);
    crate::outln!();
    crate::outln!(" Y-bending parameters(skin):");
    w("         Yc =", yct);
    w("   max Y-Yc =", ymax - yct);
    w("   min Y-Yc =", ymin - yct);
    w("      Ixx/t =", aixxt);
    let ybart = (ymax - yct).abs().max((ymin - yct).abs());
    w(" Ixx/t(Y-Yc)=", aixxt / ybart);
    crate::outln!();
    w("      J/t   =", ajt);
}

fn bendump2(x: &[f64], y: &[f64], t: &[f64], n: usize) {
    let dtr = 1.0f64.atan() / 45.0;
    let o = ijsect(n, x, y);
    let (xmin, xmax, ymin, ymax) = (o[2], o[3], o[4], o[5]);
    let (aj, ajt) = (o[14], o[15]);
    let (area, xcena, ycena, ei11a, ei22a, apx1a, _apx2a) = crate::geom::aecalc(n, x, y, t, 1);
    let (slen, xcent, ycent, ei11t, ei22t, apx1t, _apx2t) = crate::geom::aecalc(n, x, y, t, 2);
    let w = |a: &str, v: f64| crate::outln!(" {}{}", a, g14(v));
    crate::outln!();
    w("Area =", area);
    w("Slen =", slen);
    crate::outln!();
    crate::outln!(" X-bending parameters:");
    w("solid centroid Xc=", xcena);
    w("skin  centroid Xc=", xcent);
    w(" solid max X-Xc  =", xmax - xcena);
    w(" solid min X-Xc  =", xmin - xcena);
    w(" skin  max X-Xc  =", xmax - xcent);
    w(" skin  min X-Xc  =", xmin - xcent);
    w("     solid Iyy   =", ei22a);
    w("     skin  Iyy/t =", ei22t);
    let xbara = (xmax - xcena).abs().max((xmin - xcena).abs());
    let xbart = (xmax - xcent).abs().max((xmin - xcent).abs());
    w(" solid Iyy/(X-Xc)=", ei22a / xbara);
    w(" skin Iyy/t(X-Xc)=", ei22t / xbart);
    crate::outln!();
    crate::outln!(" Y-bending parameters:");
    w("solid centroid Yc=", ycena);
    w("skin  centroid Yc=", ycent);
    w(" solid max Y-Yc  =", ymax - ycena);
    w(" solid min Y-Yc  =", ymin - ycena);
    w(" skin  max Y-Yc  =", ymax - ycent);
    w(" skin  min Y-Yc  =", ymin - ycent);
    w("     solid Ixx   =", ei11a);
    w("     skin  Ixx/t =", ei11t);
    let ybara = (ymax - ycena).abs().max((ymin - ycena).abs());
    let ybart = (ymax - ycent).abs().max((ymin - ycent).abs());
    w(" solid Ixx/(Y-Yc)=", ei11a / ybara);
    w(" skin Ixx/t(Y-Yc)=", ei11t / ybart);
    crate::outln!();
    w(" solid principal axis angle (deg ccw) =", apx1a / dtr);
    w(" skin  principal axis angle (deg ccw) =", apx1t / dtr);
    crate::outln!();
    w("    solid J     =", aj);
    w("    skin  J/t   =", ajt);
}

/// MSAVE: replace/add element in an MSES multi-element coordinate file.
fn msave(s: &mut Session, arg: &str) {
    let fname = if arg.is_empty() || arg.starts_with(' ') {
        fname64(&asks(s, "Enter output filename for element replacement^"))
    } else {
        fname64(arg)
    };
    let txt = match std::fs::read(&fname) {
        Ok(b) => String::from_utf8_lossy(&b).to_string(),
        Err(_) => {
            crate::outln!(" Old file OPEN error.  Airfoil not saved.");
            return;
        }
    };
    let lines: Vec<String> = txt.lines().map(|l| l.to_string()).collect();
    if lines.len() < 2 {
        crate::outln!(" Old file READ error.  Airfoil not saved.");
        return;
    }
    let name1: String = substr(&lines[0], 1, 80);
    let ispars1: String = substr(&lines[1], 1, 80);
    let lastnb = |t: &str| {
        let b = t.as_bytes();
        let mut k = 80;
        while k >= 2 {
            if b[k - 1] != b' ' {
                break;
            }
            k -= 1;
        }
        k
    };
    let nn1 = lastnb(&name1);
    let ni1 = lastnb(&ispars1);
    const NEX: usize = 5;
    let mut xt: Vec<Vec<f64>> = vec![Vec::new(); NEX + 2];
    let mut yt: Vec<Vec<f64>> = vec![Vec::new(); NEX + 2];
    let mut rd = ListReader::new(lines[2..].to_vec());
    let mut nel;
    let mut iel = 1usize;
    let mut i;
    'outer: loop {
        // iel = 1..NEX
        i = 1usize;
        loop {
            match rd.read_reals(2) {
                ListRead::Ok(v) => {
                    if v.len() < 2 {
                        i = 1;
                        break 'outer;
                    }
                    if v[0] == 999.0 {
                        break;
                    }
                    xt[iel].push(v[0]);
                    yt[iel].push(v[1]);
                    i += 1;
                }
                ListRead::Eof => break 'outer,
                ListRead::Err => {
                    crate::outln!(" Old file READ error.  Airfoil not saved.");
                    return;
                }
            }
        }
        iel += 1;
        if iel > NEX {
            break;
        }
    }
    if iel > NEX {
        nel = NEX;
    } else if i == 1 {
        nel = iel - 1;
    } else {
        nel = iel;
    }
    crate::outln!("\n Specified multielement airfoil has{} elements.", i_fmt(nel as i64, 2));
    let ieln = s.input.aski("Enter element to be replaced by current airfoil^", iel as i64);
    if ieln < 1 || ieln as usize > nel + 1 {
        crate::outln!(" Element number inappropriate.  Airfoil not written.");
        return;
    }
    let ieln = ieln as usize;
    if ieln == nel + 1 {
        nel += 1;
    }
    let n = s.xf.n;
    xt[ieln].clear();
    yt[ieln].clear();
    for i in 1..=n {
        let idir = if s.xf.lclock { n - i + 1 } else { i };
        xt[ieln].push(s.xf.x[idir]);
        yt[ieln].push(s.xf.y[idir]);
    }
    let mut out = String::new();
    out += &substr(&name1, 1, nn1);
    out.push('\n');
    out += &substr(&ispars1, 1, ni1);
    out.push('\n');
    for ie in 1..=nel {
        for k in 0..xt[ie].len() {
            out += &format!(" {}{}\n", g_fmt(xt[ie][k], 15, 7), g_fmt(yt[ie][k], 15, 7));
        }
        if ie < nel {
            out += "  999.0  999.0\n";
        }
    }
    let _ = std::fs::write(&fname, eol(&out));
}

/// WRTDEF: settings file text.
pub fn wrtdef(s: &Session) -> String {
    let e = &s.extra_state;
    let xf = &s.xf;
    let f = |v: f64| f_fmt(v, 9, 4);
    let b9 = "         ";
    let lcolor = e.idevrp == 4;
    let l2 = |b: bool| format!(" {}", lbool(b));
    let mut o = String::new();
    o += &format!(" {}    {}{}{} | Npan    PPanel  TErat  REFrat\n", i_fmt(xf.npan as i64, 5), f(xf.cvpar), f(xf.cterat), f(xf.ctrrat));
    o += &format!(" {}{}{}{} | XrefS1  XrefS2  XrefP1 XrefP2\n", f(xf.xsref1), f(xf.xsref2), f(xf.xpref1), f(xf.xpref2));
    o += &format!(" {}{}{}{} | Size    plotAR  CHsize ScrnFr\n", f(e.size), f(e.plotar), f(e.ch), f(e.scrnfr));
    o += &format!(" {}{}{}{} | Xpage   Ypage   Xmargn Ymargn\n", f(e.xpage), f(e.ypage), f(e.xmarg), f(e.ymarg));
    o += &format!(" {}       {}       {}{} | Lcolor  Lcursor\n", l2(lcolor), l2(e.lcurs), b9, b9);
    o += &format!(" {}{}{}{} | CPmax   CPmin   CPdel\n", f(e.cpmax), f(e.cpmin), f(e.cpdel), b9);
    o += &format!(" {}{}{}{} | XoffAir ScalAir BLUwt\n", f(e.xofair), f(e.facair), f(e.uprwt), b9);
    let c = &e.cpolplf;
    o += &format!(" {}{}{}{} | CLmin   CLmax   CLdel\n", f(c[ICL][1]), f(c[ICL][2]), f(c[ICL][3]), b9);
    o += &format!(" {}{}{}{} | CDmin   CDmax   CDdel\n", f(c[ICD][1]), f(c[ICD][2]), f(c[ICD][3]), b9);
    o += &format!(" {}{}{}{} | ALmin   ALmax   ALdel\n", f(c[IAL][1]), f(c[IAL][2]), f(c[IAL][3]), b9);
    o += &format!(" {}{}{}{} | CMmin   CMmax   CMdel\n", f(c[ICM][1]), f(c[ICM][2]), f(c[ICM][3]), b9);
    o += &format!(" {}      {}{}{} | MAtype  Mach    Vaccel\n", i_fmt(xf.matyp as i64, 3), f(xf.minf1), f(xf.vaccel), b9);
    o += &format!(" {}      {}{}{} | REtype  Re/10^6 Ncrit\n", i_fmt(xf.retyp as i64, 3), f(xf.reinf1 / 1.0e6), f(xf.acrit), b9);
    o += &format!(" {}{}{}{} | XtripT  XtripB\n", f(xf.xstrip[1]), f(xf.xstrip[2]), b9, b9);
    o
}

fn parse_logical(t: &str) -> Option<bool> {
    let u = t.trim_start_matches('.').to_ascii_uppercase();
    match u.chars().next() {
        Some('T') => Some(true),
        Some('F') => Some(false),
        _ => None,
    }
}

/// GETDEF: read settings file. `lask` asks first (start-up behaviour).
pub fn getdef(s: &mut Session, filnam: &str, lask: bool) {
    let txt = match std::fs::read(filnam.trim_end()) {
        Ok(b) => String::from_utf8_lossy(&b).to_string(),
        Err(_) => {
            if !lask {
                crate::outln!("\n File  {}  not found", filnam);
            } else {
                crate::outln!("\n File  {}  not found", filnam);
            }
            return;
        }
    };
    if lask {
        crate::outln!("\n  Read settings from file  {} ?  Y", filnam);
        let ans = readln(s);
        if ans.starts_with(['N', 'n']) {
            return;
        }
    }
    let mut rd = ListReader::new(txt.lines().map(|l| l.to_string()).collect());
    let ok = (|| -> Option<()> {
        let rr = |rd: &mut ListReader, n: usize| -> Option<Vec<f64>> {
            match rd.read_reals(n) {
                ListRead::Ok(v) if v.len() == n => Some(v),
                _ => None,
            }
        };
        let v = match rd.read_tokens(4) {
            ListRead::Ok(t) if t.len() == 4 => t,
            _ => return None,
        };
        let npan = list_int1(&v[0])?;
        let cv: Vec<f64> = v[1..].iter().map(|t| parse_fortran_real(t)).collect::<Option<Vec<f64>>>()?;
        s.xf.npan = npan.max(0) as usize;
        s.xf.cvpar = cv[0];
        s.xf.cterat = cv[1];
        s.xf.ctrrat = cv[2];
        let v = rr(&mut rd, 4)?;
        s.xf.xsref1 = v[0];
        s.xf.xsref2 = v[1];
        s.xf.xpref1 = v[2];
        s.xf.xpref2 = v[3];
        let v = rr(&mut rd, 4)?;
        let e = &mut s.extra_state;
        e.size = v[0];
        e.plotar = v[1];
        e.ch = v[2];
        e.scrnfr = v[3];
        let v = rr(&mut rd, 4)?;
        e.xpage = v[0];
        e.ypage = v[1];
        e.xmarg = v[2];
        e.ymarg = v[3];
        let t = match rd.read_tokens(2) {
            ListRead::Ok(t) if t.len() == 2 => t,
            _ => return None,
        };
        let lcolor = parse_logical(&t[0])?;
        e.lcurs = parse_logical(&t[1])?;
        let v = rr(&mut rd, 3)?;
        let e = &mut s.extra_state;
        e.cpmax = v[0];
        e.cpmin = v[1];
        e.cpdel = v[2];
        let v = rr(&mut rd, 3)?;
        let e = &mut s.extra_state;
        e.xofair = v[0];
        e.facair = v[1];
        e.uprwt = v[2];
        for kv in [ICL, ICD, IAL, ICM] {
            let v = rr(&mut rd, 3)?;
            for j in 0..3 {
                s.extra_state.cpolplf[kv][j + 1] = v[j];
            }
        }
        let t = match rd.read_tokens(3) {
            ListRead::Ok(t) if t.len() == 3 => t,
            _ => return None,
        };
        s.xf.matyp = list_int1(&t[0])? as i32;
        s.xf.minf1 = parse_fortran_real(&t[1])?;
        s.xf.vaccel = parse_fortran_real(&t[2])?;
        let t = match rd.read_tokens(3) {
            ListRead::Ok(t) if t.len() == 3 => t,
            _ => return None,
        };
        s.xf.retyp = list_int1(&t[0])? as i32;
        let rmill = parse_fortran_real(&t[1])?;
        s.xf.acrit = parse_fortran_real(&t[2])?;
        let v = rr(&mut rd, 2)?;
        s.xf.xstrip[1] = v[0];
        s.xf.xstrip[2] = v[1];
        s.extra_state.idevrp = if lcolor { 4 } else { 2 };
        s.xf.reinf1 = rmill * 1.0e6;
        Some(())
    })();
    if ok.is_none() {
        crate::outln!("\n File  {}  read error\n Settings may have been changed", filnam);
        return;
    }
    let mut r = 0.0;
    s.xf.mrcl(1.0, &mut r);
    s.xf.comset();
    crate::outln!("\n Default parameters read in from file  {}:\n", filnam);
    let t = wrtdef(s);
    crate::out!("{}", t);
}

/// OPLSET (PLOP): plot option menu; values are kept for WDEF.
pub fn oplset(s: &mut Session) {
    loop {
        {
            let e = &s.extra_state;
            let lgraph = e.idev >= 1;
            let lcolor = e.idevrp == 4;
            crate::outln!(
                " ...............................................\n\n  G raphics-enable flag:        {}\n  S ize of plot object         {}\"\n  A spect ratio of plot object {}\n  P age dimensions             {} x{}\"\n  M argins from page edges     {}\",{}\"\n  F ont size (relative)        {}\n  W indow/screen size fraction {}\n  B lowup input method:        {}\n  O rientation of plot:        {}\n  C olor PostScript output?     {}",
                lbool(lgraph),
                f_fmt(e.size, 6, 2),
                f_fmt(e.plotar, 8, 4),
                f_fmt(e.xpage, 6, 2),
                f_fmt(e.ypage, 6, 2),
                f_fmt(e.xmarg, 6, 2),
                f_fmt(e.ymarg, 6, 2),
                f_fmt(e.ch, 8, 4),
                f_fmt(e.scrnfr, 8, 4),
                if e.lcurs { "Cursor    " } else { "Keyboard  " },
                if e.lland { "Landscape " } else { "Portrait  " },
                lbool(lcolor)
            );
        }
        let (com, arg) = s.input.askc("      Option, Value   (or <Return>) ^");
        if s.input.eof && com == "QUIT" {
            return;
        }
        let (r, _) = args(&arg);
        let var = com.chars().next().unwrap_or(' ');
        let e = &mut s.extra_state;
        match var.to_ascii_uppercase() {
            '0' | ' ' => return,
            'G' => e.idev = if e.idev == 0 { 1 } else { 0 },
            'S' => {
                if !r.is_empty() {
                    e.size = r[0];
                } else {
                    let v = s.input.askr("Enter size (in)^", s.extra_state.size);
                    s.extra_state.size = v;
                }
            }
            'A' => {
                if !r.is_empty() {
                    e.plotar = r[0];
                } else {
                    let v = s.input.askr("Enter aspect ratio^", s.extra_state.plotar);
                    s.extra_state.plotar = v;
                }
            }
            'P' => {
                if r.len() >= 2 {
                    e.xpage = r[0];
                    e.ypage = r[1];
                } else if r.len() == 1 {
                    e.xpage = r[0];
                    let v = s.input.askr("Enter page Y dimension (in)^", s.extra_state.ypage);
                    s.extra_state.ypage = v;
                } else {
                    let v = s.input.askr("Enter page X dimension (in)^", s.extra_state.xpage);
                    s.extra_state.xpage = v;
                    let v = s.input.askr("Enter page Y dimension (in)^", s.extra_state.ypage);
                    s.extra_state.ypage = v;
                }
            }
            'M' => {
                if r.len() >= 2 {
                    e.xmarg = r[0];
                    e.ymarg = r[1];
                } else if r.len() == 1 {
                    e.xmarg = r[0];
                    let v = s.input.askr("Enter page Y margin (in)^", s.extra_state.ymarg);
                    s.extra_state.ymarg = v;
                } else {
                    let v = s.input.askr("Enter page X margin (in)^", s.extra_state.xmarg);
                    s.extra_state.xmarg = v;
                    let v = s.input.askr("Enter page Y margin (in)^", s.extra_state.ymarg);
                    s.extra_state.ymarg = v;
                }
            }
            'F' => {
                if !r.is_empty() {
                    e.ch = r[0];
                } else {
                    let v = s.input.askr("Enter character font size^", s.extra_state.ch);
                    s.extra_state.ch = v;
                }
            }
            'W' => {
                if !r.is_empty() {
                    e.scrnfr = r[0];
                } else {
                    let v = s.input.askr("Enter window/screen size fraction^", s.extra_state.scrnfr);
                    s.extra_state.scrnfr = v;
                }
            }
            'B' => e.lcurs = !e.lcurs,
            'O' => {
                e.lland = !e.lland;
                crate::outln!();
                crate::outln!(" Swapping X,Y page dimensions");
                std::mem::swap(&mut e.xpage, &mut e.ypage);
            }
            'C' => {
                let lcolor = !(e.idevrp == 4);
                e.idevrp = if lcolor { 4 } else { 2 };
            }
            _ => crate::outln!(" *** Item not recognized ***"),
        }
    }
}
