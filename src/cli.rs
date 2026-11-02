//! XFoil-compatible command interpreter (top level, OPER, VPAR, PPAR).
//! Reads commands from stdin exactly like XFoil so existing batch scripts
//! work unchanged.

use crate::geom::getflt;
use crate::polar::Polar;
use crate::xfoil::{XFoil, DTOR, IQX};
use std::io::BufRead;

/// Line-oriented input source shared by all menus.
pub struct Input {
    lines: Box<dyn Iterator<Item = String>>,
    pub eof: bool,
}

impl Input {
    pub fn stdin() -> Self {
        let stdin = std::io::stdin();
        let it = stdin.lock().lines().map_while(Result::ok).collect::<Vec<_>>().into_iter();
        Input { lines: Box::new(it), eof: false }
    }
    pub fn from_string(s: &str) -> Self {
        let v: Vec<String> = s.lines().map(|l| l.to_string()).collect();
        Input { lines: Box::new(v.into_iter()), eof: false }
    }
    /// Read a raw line after printing a prompt (ASKS: FORMAT(/A,'   s>  ',$)).
    pub fn ask(&mut self, prompt: &str) -> String {
        self.prompt_line(prompt, "   s>  ")
    }
    fn prompt_line(&mut self, prompt: &str, tag: &str) -> String {
        let p = match prompt.find('^') {
            Some(k) if k > 0 => &prompt[..k],
            _ => prompt,
        };
        crate::out!("\n{}{}", p, tag);
        crate::flush_if_tty();
        match self.lines.next() {
            Some(l) => l,
            None => {
                self.eof = true;
                String::new()
            }
        }
    }
    /// Read command + argument (ASKC).
    pub fn askc(&mut self, prompt: &str) -> (String, String) {
        let p = match prompt.find('^') {
            Some(k) if k > 0 => &prompt[..k],
            _ => prompt,
        };
        crate::out!("\n{}   c>  ", p);
        crate::flush_if_tty();
        let line = match self.lines.next() {
            Some(l) => l,
            None => {
                self.eof = true;
                return ("QUIT".into(), String::new());
            }
        };
        split_command(&line)
    }
    /// Read a raw line without prompting (Fortran READ(*,'(A)')).
    pub fn read_line(&mut self) -> String {
        match self.lines.next() {
            Some(l) => l,
            None => {
                self.eof = true;
                String::new()
            }
        }
    }
    /// ASKR: blank input keeps the current value `def`.
    pub fn askr(&mut self, prompt: &str, def: f64) -> f64 {
        loop {
            let l = self.prompt_line(prompt, "   r>  ");
            if self.eof {
                return def;
            }
            if l.trim().is_empty() {
                return def;
            }
            if let Some(v) = getflt(&l, 1) {
                if !v.is_empty() {
                    return v[0];
                }
            }
        }
    }
    /// ASKI: blank input keeps the current value `def`.
    pub fn aski(&mut self, prompt: &str, def: i64) -> i64 {
        loop {
            let l = self.prompt_line(prompt, "   i>  ");
            if self.eof {
                return def;
            }
            if l.trim().is_empty() {
                return def;
            }
            if let Some(v) = getflt(&l, 1) {
                if !v.is_empty() {
                    return v[0] as i64;
                }
            }
        }
    }
    /// ASKL
    pub fn askl(&mut self, prompt: &str) -> bool {
        loop {
            let l = self.prompt_line(prompt, " y/n>  ");
            if self.eof {
                return false;
            }
            match l.chars().next() {
                Some('Y') | Some('y') => return true,
                Some('N') | Some('n') => return false,
                _ => {}
            }
        }
    }
}

/// Split "comand arguments" into (COMAND uppercased <=4 chars, COMARG).
pub fn split_command(line: &str) -> (String, String) {
    // ASKC: LINE is CHARACTER*128, leading blanks stripped; the command ends
    // at the first blank, "+", "-", ".", "," or numeral.
    let l: String = line.chars().take(128).collect::<String>().replace('\t', " ");
    let l = l.trim_start_matches(' ');
    let k = l.find(|c: char| c == ' ' || c == '-' || c == '+' || c == '.' || c == ',' || c.is_ascii_digit());
    let (c, rest): (String, &str) = match k {
        Some(0) => (l.chars().take(4).collect(), l),
        Some(k) => (l[..k].chars().take(4).collect(), &l[k..]),
        None => {
            // no blank: command is the first 4 characters (K=5)
            let c: String = l.chars().take(4).collect();
            let rest = if l.len() > 4 && l.is_char_boundary(4) { &l[4..] } else { "" };
            (c, rest)
        }
    };
    let c = c.trim_end_matches(' ').to_ascii_uppercase();
    // CARGS = LINE(K:128), STRIPped
    (c, rest.trim_matches(' ').to_string())
}

/// Parsed numeric arguments from a command argument string.
pub fn args(comarg: &str) -> (Vec<f64>, Vec<i64>) {
    let r = getflt(comarg, 20).unwrap_or_default();
    let i = r.iter().map(|&v| if v.abs() > 2.1e9 { 1 << 30 } else { v as i64 }).collect();
    (r, i)
}

/// Whole interactive session state.
pub struct Session {
    pub xf: XFoil,
    pub input: Input,
    pub polars: Vec<Polar>,
    pub ipact: usize,
    pub lclock_write: bool,
    pub gdes_state: crate::gdes::GdesState,
    pub qdes_state: crate::qdes::QdesState,
    pub mdes_state: crate::mdes::MdesState,
    pub extra_state: crate::opextra::ExtraState,
    comold: String,
    argold: String,
    aa: (f64, f64, f64),
    clseq: (f64, f64, f64),
}

const ANGTOL: f64 = 40.0;

impl Session {
    pub fn new(input: Input) -> Self {
        let mut xf = XFoil::new();
        xf.verbose = true;
        Session {
            xf,
            input,
            polars: Vec::new(),
            ipact: 0,
            lclock_write: false,
            gdes_state: Default::default(),
            qdes_state: Default::default(),
            mdes_state: Default::default(),
            extra_state: Default::default(),
            comold: "****".into(),
            argold: String::new(),
            aa: (0.0, 0.0, 0.0),
            clseq: (0.0, 0.0, 0.0),
        }
    }

    pub fn banner() {
        crate::outln!();
        crate::outln!(" ===================================================");
        crate::outln!("  FOILEY (XFOIL 6.97 compatible, Rust port)");
        crate::outln!("  Multithreaded + SIMD panel/BL solver");
        crate::outln!(" ===================================================");
    }

    /// Top-level menu loop (PROGRAM XFOIL).
    pub fn run(&mut self, initial_file: Option<String>) {
        Self::banner();
        crate::opextra::getdef(self, "xfoil.def", true);
        if let Some(f) = initial_file {
            self.cmd_load(&f);
        }
        loop {
            let (com, arg) = self.input.askc(" XFOIL");
            if self.input.eof && com == "QUIT" {
                break;
            }
            let (r, i) = args(&arg);
            match com.as_str() {
                "" => {}
                "?" => self.help_top(),
                "QUIT" | "Q" => break,
                "OPER" => self.oper(),
                "NACA" => {
                    let ides = if !i.is_empty() { i[0] } else { self.input.aski("Enter NACA 4 or 5-digit airfoil designation^", 0) };
                    self.xf.naca(ides as i32);
                }
                "LOAD" => self.cmd_load(&arg),
                "PANE" => self.xf.pangen(true),
                "PCOP" => self.xf.abcopy(true),
                "PPAR" => self.ppar(),
                "NORM" => {
                    self.xf.lnorm = !self.xf.lnorm;
                    if self.xf.lnorm {
                        crate::outln!(" Loaded airfoil will  be normalized");
                    } else {
                        crate::outln!(" Loaded airfoil won't be normalized");
                    }
                }
                "SAVE" => self.save(1, &arg),
                "PSAV" => self.save(0, &arg),
                "USAV" => self.save(-1, &arg),
                "ISAV" => self.save(2, &arg),
                "REVE" => {
                    self.xf.lclock = !self.xf.lclock;
                    if self.xf.lclock {
                        crate::outln!(" Airfoil will be written in clockwise order");
                    } else {
                        crate::outln!(" Airfoil will be written in counterclockwise order");
                    }
                }
                "DELI" => {
                    let mut ninput = i.len();
                    let mut kdnew = self.xf.kdelim as i64;
                    loop {
                        if ninput >= 1 {
                            kdnew = i[0];
                        } else {
                            crate::outln!(
                                "\n  --------------------------\n   0  blank\n   1  comma\n   2  tab\n\n  currently, delimiter ={:2}",
                                self.xf.kdelim
                            );
                            kdnew = self.input.aski("Enter new delimiter", kdnew);
                            if self.input.eof {
                                break;
                            }
                        }
                        if (0..=2).contains(&kdnew) {
                            self.xf.kdelim = kdnew as i32;
                            break;
                        }
                        ninput = 0;
                    }
                }
                "XYCM" => {
                    if r.len() >= 2 {
                        self.xf.xcmref = r[0];
                        self.xf.ycmref = r[1];
                    } else {
                        self.xf.xcmref = self.input.askr("Enter new CM reference X^", self.xf.xcmref);
                        self.xf.ycmref = self.input.askr("Enter new CM reference Y^", self.xf.ycmref);
                    }
                }
                "PLOP" => crate::opextra::oplset(self),
                "Z" | "U" => {}
                _ => {
                    if !crate::opextra::top_level(self, &com, &arg) {
                        crate::extra::top_level(self, &com, &arg);
                    }
                }
            }
        }
    }

    fn help_top(&self) {
        crate::outln!(
            "\n   QUIT    Exit program\n\n  .OPER    Direct operating point(s)\n  .MDES    Complex mapping design routine\n  .QDES    Surface speed design routine\n  .GDES    Geometry design routine\n\n   SAVE f  Write airfoil to labeled coordinate file\n   PSAV f  Write airfoil to plain coordinate file\n   ISAV f  Write airfoil to ISES coordinate file\n   MSAV f  Write airfoil to MSES coordinate file\n   REVE    Reverse written-airfoil node ordering\n   DELI i  Change written-airfoil file delimiters\n\n   LOAD f  Read buffer airfoil from coordinate file\n   NACA i  Set NACA 4,5-digit airfoil and buffer airfoil\n   INTE    Set buffer airfoil by interpolating two airfoils\n   NORM    Buffer airfoil normalization toggle\n   HALF    Halve the number of points in buffer airfoil\n   XYCM rr Change CM reference location, currently {:8.5}{:8.5}\n\n   BEND    Display structural properties of current airfoil\n\n   PCOP    Set current-airfoil panel nodes directly from buffer airfoil points\n   PANE    Set current-airfoil panel nodes ({:4} ) based on curvature\n  .PPAR    Show/change paneling\n\n  .PLOP    Plotting options\n\n   WDEF f  Write  current-settings file\n   RDEF f  Reread current-settings file\n   NAME s  Specify new airfoil name\n   NINC    Increment name version number",
            self.xf.xcmref, self.xf.ycmref, self.xf.npan
        );
    }

    pub fn cmd_load(&mut self, arg: &str) {
        let fname = if arg.trim().is_empty() { self.input.ask("Enter filename^") } else { arg.to_string() };
        let itype = self.xf.load(&fname);
        if itype == 1 {
            self.xf.name = self.input.ask("Enter airfoil name^");
        }
        if itype > 0 && self.xf.nb > 0 {
            self.xf.abcopy(true);
            let (imax, amax) = crate::geom::cang(&self.xf.x, &self.xf.y, self.xf.n);
            if amax.abs() > ANGTOL {
                crate::outln!(
                    "\n WARNING: Poor input coordinate distribution\n          Excessive panel angle{:7.1}  at i ={:4}\n          Repaneling with PANE and/or PPAR suggested\n           (doing GDES,CADD before repaneling _may_\n            improve excessively coarse LE spacing",
                    amax, imax
                );
            }
        }
    }

    /// SAVE (xfoil.f)
    fn save(&mut self, iftyp: i32, arg: &str) {
        let fname = if arg.is_empty() || arg.starts_with(' ') { self.input.ask("Enter output filename^") } else { arg.to_string() };
        let fname: String = fname.chars().take(64).collect::<String>().trim().to_string();
        if std::path::Path::new(&fname).exists() {
            crate::outln!();
            crate::outln!(" Output file exists.  Overwrite?  Y");
            let ans = self.input.read_line();
            if ans.starts_with(['N', 'n']) {
                crate::outln!(" Current airfoil not saved.");
                return;
            }
        }
        let s = crate::opextra::save_text(&self.xf, iftyp);
        if fname.is_empty() || std::fs::write(&fname, crate::polar::eol(&s)).is_err() {
            crate::outln!(" Bad filename.");
            crate::outln!(" Current airfoil not saved.");
        }
    }

    fn ppar(&mut self) {
        if self.xf.nb <= 1 {
            crate::outln!(" GETPAN: Buffer airfoil not available.");
            return;
        }
        let mut lchange = false;
        let show = |xf: &XFoil| {
            crate::outln!(
                "\n    Present paneling parameters...\n  N  i   Number of panel nodes      {:4}\n  P  r   Panel bunching parameter   {:6.3}\n  T  r   TE/LE panel density ratio  {:6.3}\n  R  r   Refined area/LE  panel density ratio  {:6.3}\n  XT rr  Top    side refined area x/c limits   {:6.3}{:6.3}\n  XB rr  Bottom side refined area x/c limits   {:6.3}{:6.3}\n  Z oom\n  U nzoom",
                xf.npan, xf.cvpar, xf.cterat, xf.ctrrat, xf.xsref1, xf.xsref2, xf.xpref1, xf.xpref2
            );
        };
        show(&self.xf);
        loop {
            let (var, arg) = self.input.askc("Change what ? (<cr> if nothing else)");
            if self.input.eof {
                return;
            }
            let (r, i) = args(&arg);
            match var.as_str() {
                "" => {
                    if lchange {
                        self.xf.pangen(false);
                        if self.xf.n > 0 {
                            let (imax, amax) = crate::geom::cang(&self.xf.x, &self.xf.y, self.xf.n);
                            crate::outln!(
                                "\n Maximum panel corner angle ={:7.3}   at  i,x,y  = {:3}{:9.4}{:9.4}",
                                amax, imax, self.xf.x[imax], self.xf.y[imax]
                            );
                        }
                        lchange = false;
                        show(&self.xf);
                        continue;
                    }
                    return;
                }
                "N" => {
                    self.xf.npan = if !i.is_empty() { i[0] as usize } else { self.input.aski("Enter number of panel nodes^", self.xf.npan as i64) as usize };
                    if self.xf.npan > IQX - 6 {
                        self.xf.npan = IQX - 6;
                        crate::outln!("  Number of panel nodes reduced to array limit:{:4}", self.xf.npan);
                    }
                    lchange = true;
                }
                "P" => {
                    self.xf.cvpar = if !r.is_empty() { r[0] } else { self.input.askr("Enter panel bunching parameter (0 to ~1)^", self.xf.cvpar) };
                    lchange = true;
                }
                "T" => {
                    self.xf.cterat = if !r.is_empty() { r[0] } else { self.input.askr("Enter TE/LE panel density ratio^", self.xf.cterat) };
                    lchange = true;
                }
                "R" => {
                    self.xf.ctrrat = if !r.is_empty() { r[0] } else { self.input.askr("Enter refined-area panel density ratio^", self.xf.ctrrat) };
                    lchange = true;
                }
                "XT" => {
                    if r.len() >= 2 {
                        self.xf.xsref1 = r[0];
                        self.xf.xsref2 = r[1];
                    } else {
                        self.xf.xsref1 = self.input.askr("Enter left   top   side refinement limit^", self.xf.xsref1);
                        self.xf.xsref2 = self.input.askr("Enter right  top   side refinement limit^", self.xf.xsref2);
                    }
                    lchange = true;
                }
                "XB" => {
                    if r.len() >= 2 {
                        self.xf.xpref1 = r[0];
                        self.xf.xpref2 = r[1];
                    } else {
                        self.xf.xpref1 = self.input.askr("Enter left  bottom side refinement limit^", self.xf.xpref1);
                        self.xf.xpref2 = self.input.askr("Enter right bottom side refinement limit^", self.xf.xpref2);
                    }
                    lchange = true;
                }
                "Z" | "U" => {}
                _ => {
                    crate::outln!("\n ***  Input not recognized  ***");
                    show(&self.xf);
                }
            }
        }
    }

    fn mrshow(&self, lm: bool, lr: bool) {
        let xf = &self.xf;
        if lm || lr {
            crate::outln!();
        }
        let suf = |t: i32| match t {
            2 => " / sqrt(CL)",
            3 => " / CL",
            _ => "",
        };
        if lm && (1..=3).contains(&xf.matyp) {
            crate::outln!(" M  ={}{}", crate::fmt::f_fmt(xf.minf1, 10, 4), suf(xf.matyp));
        }
        if lr && (1..=3).contains(&xf.retyp) {
            crate::outln!(" Re ={}{}", crate::fmt::g_fmt(xf.reinf1, 12, 4), suf(xf.retyp));
        }
    }

    fn after_point(&mut self) {
        if self.xf.lpacc && (self.xf.lvconv || !self.xf.lvisc) {
            crate::opextra::plradd(self);
            crate::opextra::plxadd(self);
        }
    }

    fn set_alpha_point(&mut self) {
        let xf = &mut self.xf;
        xf.lalfa = true;
        xf.qinf = 1.0;
        xf.specal();
        if (xf.alfa - xf.awake).abs() > 1.0e-5 {
            xf.lwake = false;
        }
        if (xf.alfa - xf.avisc).abs() > 1.0e-5 {
            xf.lvconv = false;
        }
        if (xf.minf - xf.mvisc).abs() > 1.0e-5 {
            xf.lvconv = false;
        }
    }

    /// OPER menu.
    pub fn oper(&mut self) {
        if self.xf.n == 0 {
            crate::outln!("\n ***  No airfoil available  ***");
            return;
        }
        if self.ipact != 0 {
            crate::outln!("\n  Polar{:3}  is active", self.ipact);
        }
        self.comold = "****".into();
        self.argold = String::new();
        loop {
            let prompt = match (self.xf.lvisc, self.xf.lpacc) {
                (true, true) => ".OPERva",
                (true, false) => ".OPERv",
                (false, true) => ".OPERia",
                (false, false) => ".OPERi",
            };
            let (mut com, mut arg) = self.input.askc(prompt);
            if self.input.eof && com == "QUIT" {
                // EOF: propagate quit to top level
                return;
            }
            let lrecalc;
            if com.starts_with('!') {
                if self.comold == "****" {
                    crate::outln!(" Previous .OPER command not valid");
                    continue;
                }
                com = self.comold.clone();
                arg = self.argold.clone();
                lrecalc = true;
            } else {
                lrecalc = false;
            }
            if com.is_empty() {
                self.extra_state.lplot = false;
                return;
            }
            let (r, i) = args(&arg);
            let ninput = r.len();
            match com.as_str() {
                "VISC" | "V" => {
                    if self.xf.lpacc {
                        crate::outln!("\n * Polar is being accumulated.\n * Cannot change its parameters in midstream.");
                        continue;
                    }
                    self.xf.lvisc = !self.xf.lvisc;
                    if self.xf.lvisc {
                        if ninput >= 1 {
                            self.xf.reinf1 = r[0];
                        } else if self.xf.reinf1 == 0.0 {
                            self.xf.reinf1 = self.input.askr("Enter Reynolds number^", self.xf.reinf1);
                        }
                        self.mrshow(true, true);
                    }
                    self.xf.lvconv = false;
                }
                "RE" | "R" => {
                    if self.xf.lpacc && self.xf.lvisc {
                        crate::outln!("\n * Polar is being accumulated.\n * Cannot change its parameters in midstream.");
                        continue;
                    }
                    if ninput >= 1 {
                        self.xf.reinf1 = r[0];
                    } else {
                        crate::outln!("\n Currently...");
                        self.mrshow(false, true);
                        self.xf.reinf1 = self.input.askr("Enter new Reynolds number^", self.xf.reinf1);
                    }
                    let mut rr = 0.0;
                    self.xf.mrcl(1.0, &mut rr);
                    self.xf.lvconv = false;
                }
                "MACH" | "M" => {
                    if self.xf.lpacc {
                        crate::outln!("\n * Polar is being accumulated.\n * Cannot change its parameters in midstream.");
                        continue;
                    }
                    let mut m = if ninput >= 1 { r[0] } else { self.input.askr("Enter Mach number^", self.xf.minf1) };
                    while m >= 1.0 {
                        crate::outln!(" Supersonic freestream not allowed");
                        m = self.input.askr("Enter Mach number^", 0.0);
                        if self.input.eof {
                            break;
                        }
                    }
                    self.xf.minf1 = m;
                    let mut rr = 0.0;
                    self.xf.mrcl(1.0, &mut rr);
                    self.xf.comset();
                    if self.xf.minf > 0.0 {
                        crate::outln!("\n Sonic Cp ={:10.2}      Sonic Q/Qinf ={:10.3}\n", self.xf.cpstar, self.xf.qstar / self.xf.qinf);
                    }
                    let n = self.xf.n;
                    crate::xfoil::cpcalc(n, &self.xf.qinv.clone(), self.xf.qinf, self.xf.minf, &mut self.xf.cpi);
                    if self.xf.lvisc {
                        let nn = n + self.xf.nw;
                        crate::xfoil::cpcalc(nn, &self.xf.qvis.clone(), self.xf.qinf, self.xf.minf, &mut self.xf.cpv);
                    }
                    self.xf.clcalc();
                    self.xf.cdcalc();
                    self.xf.lvconv = false;
                }
                "TYPE" | "T" => {
                    if self.xf.lpacc {
                        crate::outln!("\n * Polar is being accumulated.\n * Cannot change its parameters in midstream.");
                        continue;
                    }
                    let mut ityp = if ninput >= 1 { i[0] } else { 0 };
                    while !(1..=3).contains(&ityp) {
                        crate::outln!("\n Type   parameters held constant       varying      fixed   \n ----   ------------------------       -------   -----------\n   1    M          , Re            ..   lift     chord, vel.\n   2    M sqrt(CL) , Re sqrt(CL)   ..   vel.     chord, lift\n   3    M          , Re CL         ..   chord    lift , vel.");
                        ityp = self.input.aski("Enter type of Mach,Re variation with CL^", 1);
                        if self.input.eof {
                            break;
                        }
                    }
                    match ityp {
                        2 => {
                            self.xf.matyp = 2;
                            self.xf.retyp = 2;
                        }
                        3 => {
                            self.xf.matyp = 1;
                            self.xf.retyp = 3;
                        }
                        _ => {
                            self.xf.matyp = 1;
                            self.xf.retyp = 1;
                        }
                    }
                    self.mrshow(true, true);
                    let mut rr = 0.0;
                    self.xf.mrcl(1.0, &mut rr);
                    self.xf.comset();
                    self.xf.lvconv = false;
                }
                "ITER" => {
                    let mut it = if ninput >= 1 { i[0] } else { self.input.aski("Enter new iteration limit^", self.xf.itmax as i64) };
                    while it < 1 && !self.input.eof {
                        it = self.input.aski("Enter new iteration limit^", self.xf.itmax as i64);
                    }
                    self.xf.itmax = it.max(1) as usize;
                }
                "INIT" => {
                    self.xf.lblini = !self.xf.lblini;
                    if self.xf.lblini {
                        crate::outln!(" BLs are assumed to be initialized");
                    } else {
                        crate::outln!(" BLs will be initialized on next point");
                        self.xf.lipan = false;
                    }
                }
                "ALFA" | "A" => {
                    if !lrecalc {
                        let adeg = if ninput >= 1 { r[0] } else { self.input.askr("Enter angle of attack (deg)^", self.xf.alfa / DTOR) };
                        self.xf.adeg = adeg;
                        self.xf.alfa = DTOR * adeg;
                        self.set_alpha_point();
                    }
                    if self.xf.lvisc {
                        let it = self.xf.itmax;
                        self.xf.viscal(it);
                    }
                    self.xf.fcpmin();
                    self.extra_state.lplot = true;
                    self.after_point();
                    if self.xf.lvisc && !self.xf.lpacc && !self.xf.lvconv {
                        crate::outln!(" Type \"!\" to continue iterating");
                    }
                    self.comold = com.clone();
                    self.argold = arg.clone();
                }
                "CLI" | "CL" | "C" => {
                    let invisc = com == "CLI";
                    if !lrecalc {
                        self.xf.clspec = if ninput >= 1 {
                            r[0]
                        } else if invisc {
                            self.input.askr("Enter inviscid lift coefficient^", self.xf.clspec)
                        } else {
                            self.input.askr("Enter lift coefficient^", self.xf.clspec)
                        };
                        self.xf.lalfa = invisc;
                        self.xf.alfa = 0.0;
                        self.xf.qinf = 1.0;
                        self.xf.speccl();
                        self.xf.adeg = self.xf.alfa / DTOR;
                        let xf = &mut self.xf;
                        if (xf.alfa - xf.awake).abs() > 1.0e-5 {
                            xf.lwake = false;
                        }
                        if (xf.alfa - xf.avisc).abs() > 1.0e-5 {
                            xf.lvconv = false;
                        }
                        if (xf.minf - xf.mvisc).abs() > 1.0e-5 {
                            xf.lvconv = false;
                        }
                    }
                    if self.xf.lvisc {
                        let it = self.xf.itmax;
                        self.xf.viscal(it);
                        self.xf.adeg = self.xf.alfa / DTOR;
                    }
                    self.xf.fcpmin();
                    self.extra_state.lplot = true;
                    self.after_point();
                    if self.xf.lvisc && !self.xf.lpacc && !self.xf.lvconv {
                        crate::outln!(" Type \"!\" to continue iterating");
                    }
                    self.comold = com.clone();
                    self.argold = arg.clone();
                }
                "ASEQ" | "AS" | "CSEQ" | "CS" => {
                    let lalfa = com == "ASEQ" || com == "AS";
                    self.xf.lalfa = lalfa;
                    let npoint;
                    if lalfa {
                        let (mut a1, mut a2, mut da) = (self.aa.0, self.aa.1, self.aa.2);
                        if ninput >= 3 {
                            a1 = r[0];
                            a2 = r[1];
                            da = r[2];
                        } else if ninput >= 2 {
                            a1 = r[0];
                            a2 = r[1];
                            da = self.input.askr("Enter alfa increment   (deg)^", da);
                        } else if ninput >= 1 {
                            a1 = r[0];
                            a2 = self.input.askr("Enter last  alfa value (deg)^", a2);
                            da = self.input.askr("Enter alfa increment   (deg)^", da);
                        } else {
                            a1 = self.input.askr("Enter first alfa value (deg)^", a1);
                            a2 = self.input.askr("Enter last  alfa value (deg)^", a2);
                            da = self.input.askr("Enter alfa increment   (deg)^", da);
                        }
                        da = if a2 < a1 { -da.abs() } else { da.abs() };
                        self.aa = (a1, a2, da);
                        let (aa1, aa2, daa) = (a1 * DTOR, a2 * DTOR, da * DTOR);
                        npoint = if daa != 0.0 { ((aa2 - aa1) / daa + 0.5) as i64 + 1 } else { 1 };
                    } else {
                        let (mut c1, mut c2, mut dc) = self.clseq;
                        if ninput >= 3 {
                            c1 = r[0];
                            c2 = r[1];
                            dc = r[2];
                        } else if ninput >= 2 {
                            c1 = r[0];
                            c2 = r[1];
                            dc = self.input.askr("Enter CL increment  ^", dc);
                        } else if ninput >= 1 {
                            c1 = r[0];
                            c2 = self.input.askr("Enter last  CL value^", c2);
                            dc = self.input.askr("Enter CL increment  ^", dc);
                        } else {
                            c1 = self.input.askr("Enter first CL value^", c1);
                            c2 = self.input.askr("Enter last  CL value^", c2);
                            dc = self.input.askr("Enter CL increment  ^", dc);
                        }
                        dc = if c2 < c1 { -dc.abs() } else { dc.abs() };
                        self.clseq = (c1, c2, dc);
                        npoint = if dc != 0.0 { ((c2 - c1) / dc + 0.5) as i64 + 1 } else { 1 };
                    }
                    self.extra_state.lplot = true;
                    self.sequence(lalfa, npoint.max(0) as usize);
                    self.comold = com.clone();
                    self.argold = arg.clone();
                }
                "PACC" | "P" => crate::opextra::pacc(self, &i),
                "VPAR" => self.vpar(),
                "VELS" => {
                    crate::extra::oper_cmd(self, &com, &arg);
                    self.comold = com.clone();
                    self.argold = arg.clone();
                }
                _ => {
                    if !crate::extra::oper_cmd(self, &com, &arg) {
                        crate::outln!(" {:<4} command not recognized.  Type a \"?\" for list", com);
                    }
                }
            }
        }
    }

    fn sequence(&mut self, lalfa: bool, npoint: usize) {
        let specs: Vec<f64> = if lalfa {
            let (aa1, daa) = (self.aa.0 * DTOR, self.aa.2 * DTOR);
            (0..npoint).map(|i| aa1 + daa * i as f64).collect()
        } else {
            let (cl1, dcl) = (self.clseq.0, self.clseq.2);
            (0..npoint).map(|i| cl1 + dcl * i as f64).collect()
        };
        let itmaxs = self.xf.itmax + 5;
        let jobs = crate::sweep::default_jobs();
        // Parallel chunks only report per-point coefficients; the polar dump
        // file and flap hinge moments need the full per-point solver state,
        // so those always use the (bit-exact) sequential path.
        let parallel = jobs > 1 && self.xf.lvisc && !self.extra_state.lpfilx && !self.xf.lflap && npoint >= 8;
        let mut xf = std::mem::take(&mut self.xf);
        if parallel {
            let pts = crate::sweep::run_parallel(&mut xf, lalfa, &specs, itmaxs, jobs, |_, _| {});
            self.xf = xf;
            let keep = {
                let x = &self.xf;
                (x.adeg, x.alfa, x.cl, x.cd, x.cdp, x.cm, x.minf, x.reinf, x.xoctr, x.cpmn, x.lvconv)
            };
            for p in &pts {
                let q = &p.point;
                {
                    let x = &mut self.xf;
                    x.adeg = q.alpha;
                    x.alfa = q.alpha * DTOR;
                    x.cl = q.cl;
                    x.cd = q.cd;
                    x.cdp = q.cdp;
                    x.cm = q.cm;
                    x.minf = q.mach;
                    if q.re > 0.0 {
                        x.reinf = q.re;
                    }
                    x.xoctr[1] = q.xtr_top;
                    x.xoctr[2] = q.xtr_bot;
                    x.cpmn = q.cpmin;
                    x.lvconv = p.converged;
                }
                crate::outln!("\n       a ={:7.3}      CL ={:8.4}", q.alpha, q.cl);
                crate::outln!("      Cm ={:8.4}     CD ={:9.5}   =>   CDp ={:9.5}", q.cm, q.cd, q.cdp);
                self.after_point();
            }
            let x = &mut self.xf;
            (x.adeg, x.alfa, x.cl, x.cd, x.cdp, x.cm, x.minf, x.reinf, x.xoctr, x.cpmn, x.lvconv) = keep;
        } else {
            // Exact XFoil semantics (alpha sweeps precompute the next point's
            // wake on a helper thread; results are bit-identical).  The live
            // per-point state is swapped into the session for PLRADD/PLXADD.
            crate::sweep::run_sequential(&mut xf, lalfa, &specs, itmaxs, |st, _p| {
                std::mem::swap(&mut self.xf, st);
                self.after_point();
                std::mem::swap(&mut self.xf, st);
            });
            self.xf = xf;
        }
    }

    fn vpar(&mut self) {
        let show = |xf: &XFoil| {
            let turb = 100.0 * (-(xf.acrit + 8.43) / 2.4).exp();
            let p = xf.bl.par;
            crate::outln!(
                "\n Xtr/c     ={:8.4}    top    side\n Xtr/c     ={:8.4}    bottom side\n Ncrit     ={:8.2}   ({:6.3} % turb. level )\n Vacc      ={:8.4}\n\n Klag  ={:8.4}     Uxwt  ={:8.2}\n A     ={:8.4}     B     ={:8.4}       KCt ={:8.5}\n CtiniK={:8.4}     CtiniX={:8.4}",
                xf.xstrip[1], xf.xstrip[2], xf.acrit, turb, xf.vaccel, p.sccon, p.duxcon, p.gacon, p.gbcon, p.ctcon, p.ctrcon, p.ctrcex
            );
        };
        show(&self.xf);
        loop {
            let (com, arg) = self.input.askc("..VPAR");
            if com.is_empty() || (self.input.eof && com == "QUIT") {
                return;
            }
            let (r, _i) = args(&arg);
            let locked = self.xf.lpacc && self.xf.lvisc;
            let lockmsg = || crate::outln!("\n * Polar is being accumulated.\n * Cannot change its parameters in midstream.");
            match com.as_str() {
                "SHOW" => show(&self.xf),
                "XTR" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    if r.len() >= 2 {
                        self.xf.xstrip[1] = r[0];
                        self.xf.xstrip[2] = r[1];
                    } else {
                        self.xf.xstrip[1] = self.input.askr("Enter top    side Xtrip/c^", self.xf.xstrip[1]);
                        self.xf.xstrip[2] = self.input.askr("Enter bottom side Xtrip/c^", self.xf.xstrip[2]);
                    }
                    self.xf.lvconv = false;
                }
                "N" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    self.xf.acrit = if !r.is_empty() { r[0] } else { self.input.askr("Enter critical amplification ratio^", self.xf.acrit) };
                    self.xf.lvconv = false;
                }
                "VACC" => {
                    self.xf.vaccel = if !r.is_empty() { r[0] } else { self.input.askr("Enter viscous acceleration parameter^", self.xf.vaccel) };
                }
                "INIT" => {
                    self.xf.lblini = !self.xf.lblini;
                    if !self.xf.lblini {
                        crate::outln!(" BLs will be initialized on next point");
                        self.xf.lipan = false;
                    } else {
                        crate::outln!(" BLs are assumed to be initialized");
                    }
                }
                "LAG" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    let p = &mut self.xf.bl.par;
                    if r.len() >= 2 {
                        p.sccon = r[0];
                        p.duxcon = r[1];
                    } else {
                        let a = self.input.askr("Enter shear lag constant^", self.xf.bl.par.sccon);
                        let b = self.input.askr("Enter shear lag UxEQ weight^", self.xf.bl.par.duxcon);
                        self.xf.bl.par.sccon = a;
                        self.xf.bl.par.duxcon = b;
                    }
                }
                "GB" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    if r.len() >= 2 {
                        self.xf.bl.par.gacon = r[0];
                        self.xf.bl.par.gbcon = r[1];
                    } else {
                        let a = self.input.askr("Enter G-beta constant A^", self.xf.bl.par.gacon);
                        let b = self.input.askr("Enter G-beta constant B^", self.xf.bl.par.gbcon);
                        self.xf.bl.par.gacon = a;
                        self.xf.bl.par.gbcon = b;
                    }
                    let p = &mut self.xf.bl.par;
                    p.ctcon = 0.5 / (p.gacon * p.gacon * p.gbcon);
                }
                "CTR" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    if r.len() >= 2 {
                        self.xf.bl.par.ctrcon = r[0];
                        self.xf.bl.par.ctrcex = r[1];
                    } else {
                        let a = self.input.askr("Enter initial-Ctau constant^", self.xf.bl.par.ctrcon);
                        let b = self.input.askr("Enter initial-Ctau exponent^", self.xf.bl.par.ctrcex);
                        self.xf.bl.par.ctrcon = a;
                        self.xf.bl.par.ctrcex = b;
                    }
                }
                "CFAC" => {
                    self.xf.bl.par.cffac = if !r.is_empty() { r[0] } else { self.input.askr("Enter Cf scaling factor^", self.xf.bl.par.cffac) };
                }
                "REST" => {
                    if locked {
                        lockmsg();
                        continue;
                    }
                    self.xf.bl.par = crate::blsys::BlPar::default();
                }
                "?" => crate::outln!("\n   <cr>    Return to OPER menu\n   SHOW    Display viscous parameters\n   XTR  rr Change trip positions Xtr/c\n   N    r  Change critical amplification exponent Ncrit\n   VACC r  Change Newton solution acceleration parameter\n   INIT    BL initialization flag toggle\n\n   LAG     change lag equation constants\n   GB      change G-beta constants\n   CTR     change initial transition-Ctau constants\n   REST    restore BL calibration to baseline"),
                _ => crate::outln!(" {:<4} command not recognized.  Type a \"?\" for list", com),
            }
        }
    }
}
