//! Main XFoil solver state (port of XFOIL.INC) and top-level aerodynamic
//! routines (MRCL, COMSET, CPCALC, CLCALC, CDCALC, TECALC, SPECAL, SPECCL,
//! VISCAL).
//!
//! Arrays are 1-based (index 0 unused) to mirror the Fortran exactly.

use crate::blsys::BlSys;

pub const IQX: usize = 360;
pub const IWX: usize = IQX / 8 + 2;
pub const IBX: usize = 4 * IQX;
pub const IZX: usize = IQX + IWX;
pub const IVX: usize = IQX / 2 + IWX + 50;

pub const PI: f64 = std::f64::consts::PI;
pub const HOPI: f64 = 0.5 / PI;
pub const QOPI: f64 = 0.25 / PI;
pub const DTOR: f64 = PI / 180.0;

/// Row-major dense matrix with 1-based indexing helpers.
#[derive(Clone, Debug, Default)]
pub struct Mat {
    pub nr: usize,
    pub nc: usize,
    pub ld: usize,
    pub d: Vec<f64>,
}

impl Mat {
    pub fn new(nr: usize, nc: usize) -> Self {
        let ld = nc + 1;
        Mat { nr, nc, ld, d: vec![0.0; (nr + 1) * ld] }
    }
    #[inline(always)]
    pub fn at(&self, i: usize, j: usize) -> f64 {
        unsafe { *self.d.get_unchecked(i * self.ld + j) }
    }
    #[inline(always)]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        let ld = self.ld;
        unsafe { *self.d.get_unchecked_mut(i * ld + j) = v }
    }
    #[inline(always)]
    pub fn row(&self, i: usize) -> &[f64] {
        &self.d[i * self.ld..(i + 1) * self.ld]
    }
    #[inline(always)]
    pub fn row_mut(&mut self, i: usize) -> &mut [f64] {
        let ld = self.ld;
        &mut self.d[i * ld..(i + 1) * ld]
    }
    pub fn ensure(&mut self, nr: usize, nc: usize) {
        if self.nr != nr || self.nc != nc {
            *self = Mat::new(nr, nc);
        }
    }
}

/// Outputs of PSILIN / PSWLIN (were in COMMON in Fortran).
#[derive(Clone, Debug)]
pub struct PsiOut {
    pub psi: f64,
    pub psi_ni: f64,
    pub qtan1: f64,
    pub qtan2: f64,
    pub z_qinf: f64,
    pub z_alfa: f64,
    pub z_qdof0: f64,
    pub z_qdof1: f64,
    pub z_qdof2: f64,
    pub z_qdof3: f64,
    pub dzdg: Vec<f64>,
    pub dzdn: Vec<f64>,
    pub dqdg: Vec<f64>,
    pub dzdm: Vec<f64>,
    pub dqdm: Vec<f64>,
}

impl PsiOut {
    pub fn new() -> Self {
        PsiOut {
            psi: 0.0,
            psi_ni: 0.0,
            qtan1: 0.0,
            qtan2: 0.0,
            z_qinf: 0.0,
            z_alfa: 0.0,
            z_qdof0: 0.0,
            z_qdof1: 0.0,
            z_qdof2: 0.0,
            z_qdof3: 0.0,
            dzdg: vec![0.0; IZX + 2],
            dzdn: vec![0.0; IZX + 2],
            dqdg: vec![0.0; IZX + 2],
            dzdm: vec![0.0; IZX + 2],
            dqdm: vec![0.0; IZX + 2],
        }
    }
}

impl Default for PsiOut {
    fn default() -> Self {
        Self::new()
    }
}

fn v(n: usize) -> Vec<f64> {
    vec![0.0; n + 1]
}
fn v3(n: usize) -> [Vec<f64>; 3] {
    [v(n), v(n), v(n)]
}

/// Complete solver state.
#[derive(Clone, Debug)]
pub struct XFoil {
    pub name: String,
    pub prefix: String,
    pub ispars: String,

    // current (paneled) airfoil + wake
    pub n: usize,
    pub nw: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub xp: Vec<f64>,
    pub yp: Vec<f64>,
    pub s: Vec<f64>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub apanel: Vec<f64>,
    pub sle: f64,
    pub xle: f64,
    pub yle: f64,
    pub xte: f64,
    pub yte: f64,
    pub chord: f64,
    pub wgap: Vec<f64>,
    pub waklen: f64,

    // vorticity / sources
    pub gam: Vec<f64>,
    pub gamu: [Vec<f64>; 3],
    pub gam_a: Vec<f64>,
    pub sig: Vec<f64>,
    pub sst: f64,
    pub sst_go: f64,
    pub sst_gp: f64,
    pub gamte: f64,
    pub gamte_a: f64,
    pub sigte: f64,
    pub sigte_a: f64,
    pub dste: f64,
    pub ante: f64,
    pub aste: f64,
    pub sharp: bool,

    pub qinv: Vec<f64>,
    pub qvis: Vec<f64>,
    pub cpi: Vec<f64>,
    pub cpv: Vec<f64>,
    pub qinvu: [Vec<f64>; 3],
    pub qinv_a: Vec<f64>,

    // influence matrices
    pub aij: Mat,
    pub aijpiv: Vec<usize>,
    /// BIJ stored column-major: column j is bijc[j*(ldb)..], rows 1..=n+1
    pub bijc: Vec<f64>,
    pub ldb: usize,
    pub cij: Mat,
    pub dij: Mat,

    // operating point
    pub adeg: f64,
    pub alfa: f64,
    pub awake: f64,
    pub avisc: f64,
    pub mvisc: f64,
    pub xcmref: f64,
    pub ycmref: f64,
    pub cl: f64,
    pub cm: f64,
    pub cd: f64,
    pub cdp: f64,
    pub cdf: f64,
    pub cl_alf: f64,
    pub cl_msq: f64,
    pub psio: f64,
    pub cosa: f64,
    pub sina: f64,
    pub qinf: f64,
    pub gamma: f64,
    pub gamm1: f64,
    pub minf1: f64,
    pub minf: f64,
    pub minf_cl: f64,
    pub tklam: f64,
    pub tkl_msq: f64,
    pub cpstar: f64,
    pub qstar: f64,
    pub cpmn: f64,
    pub cpmni: f64,
    pub cpmnv: f64,
    pub xcpmni: f64,
    pub xcpmnv: f64,
    pub clspec: f64,

    // paneling parameters
    pub npan: usize,
    pub cvpar: f64,
    pub cterat: f64,
    pub ctrrat: f64,
    pub xsref1: f64,
    pub xsref2: f64,
    pub xpref1: f64,
    pub xpref2: f64,

    // flags
    pub lgamu: bool,
    pub lqinu: bool,
    pub lvisc: bool,
    pub lalfa: bool,
    pub lwake: bool,
    pub lpacc: bool,
    pub lblini: bool,
    pub lipan: bool,
    pub lqaij: bool,
    pub ladij: bool,
    pub lwdij: bool,
    pub lvconv: bool,
    /// Newton iterations taken by the last `viscal` call (its `iter` count)
    pub viscal_iters: usize,
    pub lclock: bool,
    pub lnorm: bool,
    pub lbflap: bool,
    pub lflap: bool,
    pub lgsame: bool,
    pub lscini: bool,
    pub lqspec: bool,

    pub retyp: i32,
    pub matyp: i32,
    pub itmax: usize,
    pub nseqex: usize,
    pub kdelim: i32,
    pub idamp: i32,

    // buffer airfoil
    pub nb: usize,
    pub xb: Vec<f64>,
    pub yb: Vec<f64>,
    pub xbp: Vec<f64>,
    pub ybp: Vec<f64>,
    pub sb: Vec<f64>,
    pub snew: Vec<f64>,
    pub xbf: f64,
    pub ybf: f64,
    pub xof: f64,
    pub yof: f64,
    pub hmom: f64,
    pub hfx: f64,
    pub hfy: f64,
    pub sble: f64,
    pub chordb: f64,
    pub areab: f64,
    pub radble: f64,
    pub angbte: f64,
    pub ei11ba: f64,
    pub ei22ba: f64,
    pub apx1ba: f64,
    pub apx2ba: f64,
    pub ei11bt: f64,
    pub ei22bt: f64,
    pub apx1bt: f64,
    pub apx2bt: f64,
    pub thickb: f64,
    pub cambrb: f64,

    // boundary layer
    pub xssi: [Vec<f64>; 3],
    pub uedg: [Vec<f64>; 3],
    pub uinv: [Vec<f64>; 3],
    pub mass: [Vec<f64>; 3],
    pub thet: [Vec<f64>; 3],
    pub dstr: [Vec<f64>; 3],
    pub ctau: [Vec<f64>; 3],
    pub delt: [Vec<f64>; 3],
    pub tstr: [Vec<f64>; 3],
    pub uslp: [Vec<f64>; 3],
    pub tau: [Vec<f64>; 3],
    pub dis: [Vec<f64>; 3],
    pub ctq: [Vec<f64>; 3],
    pub vti: [Vec<f64>; 3],
    pub uinv_a: [Vec<f64>; 3],
    pub reinf1: f64,
    pub reinf: f64,
    pub reinf_cl: f64,
    pub acrit: f64,
    pub xstrip: [f64; 3],
    pub xoctr: [f64; 3],
    pub yoctr: [f64; 3],
    pub xssitr: [f64; 3],
    pub iblte: [usize; 3],
    pub nbl: [usize; 3],
    pub ipan: [Vec<usize>; 3],
    pub isys: [Vec<usize>; 3],
    pub nsys: usize,
    pub itran: [usize; 3],
    pub ist: usize,
    pub tforce: [bool; 3],
    pub rmsbl: f64,
    pub rmxbl: f64,
    pub rlx: f64,
    pub vaccel: f64,
    pub imxbl: usize,
    pub ismxbl: usize,
    pub vmxbl: char,

    // Newton system
    pub va: Vec<[[f64; 2]; 3]>,
    pub vb: Vec<[[f64; 2]; 3]>,
    pub vdel: Vec<[[f64; 2]; 3]>,
    /// VM(k,l,iv) stored at ((iv*3 + k)*ldm + l), k 0-based, l/iv 1-based
    pub vm: Vec<f64>,
    pub ldm: usize,
    pub vz: [[f64; 2]; 3],
    /// mass influence matrix M[iv][jv] (ld = ldm)
    pub uem: Vec<f64>,
    pub stcoef: Vec<crate::bl::StCoef>,

    pub bl: BlSys,

    // QDES-related arrays used by PSILIN (geolin)
    pub qf0: Vec<f64>,
    pub qf1: Vec<f64>,
    pub qf2: Vec<f64>,
    pub qf3: Vec<f64>,

    /// print iteration info to stdout
    pub verbose: bool,
    /// allow rayon parallelism inside a single solve (QDCALC, GGCALC)
    pub par_inner: bool,
    /// wake QINVU values are current (set by XYWAKE/QWCALC, cleared by GGCALC)
    pub qw_valid: bool,
    /// scratch for psilin results
    pub pso: PsiOut,
}

impl Default for XFoil {
    fn default() -> Self {
        Self::new()
    }
}

impl XFoil {
    pub fn new() -> Self {
        let mut xf = XFoil {
            name: String::new(),
            prefix: String::new(),
            ispars: " -2.0  3.0  -2.5  3.5".to_string(),
            n: 0,
            nw: 0,
            x: v(IZX),
            y: v(IZX),
            xp: v(IZX),
            yp: v(IZX),
            s: v(IZX),
            nx: v(IZX),
            ny: v(IZX),
            apanel: v(IZX),
            sle: 0.0,
            xle: 0.0,
            yle: 0.0,
            xte: 0.0,
            yte: 0.0,
            chord: 0.0,
            wgap: v(IWX),
            waklen: 1.0,
            gam: v(IQX),
            gamu: v3(IQX),
            gam_a: v(IQX),
            sig: v(IZX),
            sst: 0.0,
            sst_go: 0.0,
            sst_gp: 0.0,
            gamte: 0.0,
            gamte_a: 0.0,
            sigte: 0.0,
            sigte_a: 0.0,
            dste: 0.0,
            ante: 0.0,
            aste: 0.0,
            sharp: false,
            qinv: v(IZX),
            qvis: v(IZX),
            cpi: v(IZX),
            cpv: v(IZX),
            qinvu: v3(IZX),
            qinv_a: v(IZX),
            aij: Mat::default(),
            aijpiv: vec![0; IQX + 2],
            bijc: Vec::new(),
            ldb: 0,
            cij: Mat::default(),
            dij: Mat::default(),
            adeg: 0.0,
            alfa: 0.0,
            awake: 0.0,
            avisc: 0.0,
            mvisc: 0.0,
            xcmref: 0.25,
            ycmref: 0.0,
            cl: 0.0,
            cm: 0.0,
            cd: 0.0,
            cdp: 0.0,
            cdf: 0.0,
            cl_alf: 0.0,
            cl_msq: 0.0,
            psio: 0.0,
            cosa: 1.0,
            sina: 0.0,
            qinf: 1.0,
            gamma: 1.4,
            gamm1: 0.4,
            minf1: 0.0,
            minf: 0.0,
            minf_cl: 0.0,
            tklam: 0.0,
            tkl_msq: 0.0,
            cpstar: 0.0,
            qstar: 0.0,
            cpmn: 0.0,
            cpmni: 0.0,
            cpmnv: 0.0,
            xcpmni: 0.0,
            xcpmnv: 0.0,
            clspec: 0.0,
            npan: 160,
            cvpar: 1.0,
            cterat: 0.15,
            ctrrat: 0.2,
            xsref1: 1.0,
            xsref2: 1.0,
            xpref1: 1.0,
            xpref2: 1.0,
            lgamu: false,
            lqinu: false,
            lvisc: false,
            lalfa: false,
            lwake: false,
            lpacc: false,
            lblini: false,
            lipan: false,
            lqaij: false,
            ladij: false,
            lwdij: false,
            lvconv: false,
            viscal_iters: 0,
            lclock: false,
            lnorm: false,
            lbflap: false,
            lflap: false,
            lgsame: false,
            lscini: false,
            lqspec: false,
            retyp: 1,
            matyp: 1,
            itmax: 20,
            nseqex: 4,
            kdelim: 0,
            idamp: 0,
            nb: 0,
            xb: v(IBX),
            yb: v(IBX),
            xbp: v(IBX),
            ybp: v(IBX),
            sb: v(IBX),
            snew: v(5 * IBX),
            xbf: 0.0,
            ybf: 0.0,
            xof: 0.0,
            yof: 0.0,
            hmom: 0.0,
            hfx: 0.0,
            hfy: 0.0,
            sble: 0.0,
            chordb: 0.0,
            areab: 0.0,
            radble: 0.0,
            angbte: 0.0,
            ei11ba: 0.0,
            ei22ba: 0.0,
            apx1ba: 0.0,
            apx2ba: 0.0,
            ei11bt: 0.0,
            ei22bt: 0.0,
            apx1bt: 0.0,
            apx2bt: 0.0,
            thickb: 0.0,
            cambrb: 0.0,
            xssi: v3(IVX),
            uedg: v3(IVX),
            uinv: v3(IVX),
            mass: v3(IVX),
            thet: v3(IVX),
            dstr: v3(IVX),
            ctau: v3(IVX),
            delt: v3(IVX),
            tstr: v3(IVX),
            uslp: v3(IVX),
            tau: v3(IVX),
            dis: v3(IVX),
            ctq: v3(IVX),
            vti: v3(IVX),
            uinv_a: v3(IVX),
            reinf1: 0.0,
            reinf: 0.0,
            reinf_cl: 0.0,
            acrit: 9.0,
            xstrip: [0.0, 1.0, 1.0],
            xoctr: [0.0, 1.0, 1.0],
            yoctr: [0.0; 3],
            xssitr: [0.0; 3],
            iblte: [0; 3],
            nbl: [0; 3],
            ipan: [vec![0; IVX + 1], vec![0; IVX + 1], vec![0; IVX + 1]],
            isys: [vec![0; IVX + 1], vec![0; IVX + 1], vec![0; IVX + 1]],
            nsys: 0,
            itran: [0; 3],
            ist: 0,
            tforce: [false; 3],
            rmsbl: 0.0,
            rmxbl: 0.0,
            rlx: 1.0,
            vaccel: 0.01,
            imxbl: 0,
            ismxbl: 0,
            vmxbl: ' ',
            va: Vec::new(),
            vb: Vec::new(),
            vdel: Vec::new(),
            vm: Vec::new(),
            ldm: 0,
            vz: [[0.0; 2]; 3],
            uem: Vec::new(),
            stcoef: Vec::new(),
            bl: BlSys::new(),
            qf0: v(IQX),
            qf1: v(IQX),
            qf2: v(IQX),
            qf3: v(IQX),
            verbose: false,
            par_inner: true,
            qw_valid: false,
            pso: PsiOut::new(),
        };
        let mut r = 0.0;
        xf.mrcl(1.0, &mut r);
        xf.comset();
        xf
    }

    #[inline(always)]
    pub fn vm_idx(&self, k: usize, l: usize, iv: usize) -> usize {
        (iv * 3 + k) * self.ldm + l
    }

    /// Mach and Reynolds numbers from unit-CL values (MRCL). Returns M_CLS, sets R_CLS.
    pub fn mrcl(&mut self, cls: f64, r_cls: &mut f64) -> f64 {
        let cla = cls.max(0.000001);
        if self.retyp < 1 || self.retyp > 3 {
            self.retyp = 1;
        }
        if self.matyp < 1 || self.matyp > 3 {
            self.matyp = 1;
        }
        let mut m_cls;
        match self.matyp {
            2 => {
                self.minf = self.minf1 / cla.sqrt();
                m_cls = -0.5 * self.minf / cla;
            }
            _ => {
                self.minf = self.minf1;
                m_cls = 0.0;
            }
        }
        match self.retyp {
            2 => {
                self.reinf = self.reinf1 / cla.sqrt();
                *r_cls = -0.5 * self.reinf / cla;
            }
            3 => {
                self.reinf = self.reinf1 / cla;
                *r_cls = -self.reinf / cla;
            }
            _ => {
                self.reinf = self.reinf1;
                *r_cls = 0.0;
            }
        }
        if self.minf >= 0.99 {
            if self.verbose {
                crate::outln!("\n MRCL: CL too low for chosen Mach(CL) dependence\n      Aritificially limiting Mach to  0.99");
            }
            self.minf = 0.99;
            m_cls = 0.0;
        }
        let mut rrat = 1.0;
        if self.reinf1 > 0.0 {
            rrat = self.reinf / self.reinf1;
        }
        if rrat > 100.0 {
            self.reinf = self.reinf1 * 100.0;
            *r_cls = 0.0;
        }
        m_cls
    }

    /// Karman-Tsien parameters (COMSET).
    pub fn comset(&mut self) {
        let beta = (1.0 - self.minf * self.minf).sqrt();
        let beta_msq = -0.5 / beta;
        self.tklam = self.minf * self.minf / ((1.0 + beta) * (1.0 + beta));
        self.tkl_msq = 1.0 / ((1.0 + beta) * (1.0 + beta)) - 2.0 * self.tklam / (1.0 + beta) * beta_msq;
        if self.minf == 0.0 {
            self.cpstar = -999.0;
            self.qstar = 999.0;
        } else {
            let g = self.gamma;
            let gm1 = self.gamm1;
            let m2 = self.minf * self.minf;
            self.cpstar = 2.0 / (g * m2) * (((1.0 + 0.5 * gm1 * m2) / (1.0 + 0.5 * gm1)).powf(g / gm1) - 1.0);
            self.qstar = self.qinf / self.minf * ((1.0 + 0.5 * gm1 * m2) / (1.0 + 0.5 * gm1)).sqrt();
        }
    }

    /// TE gap geometry and TE panel strengths (TECALC).
    pub fn tecalc(&mut self) {
        let n = self.n;
        let dxte = self.x[1] - self.x[n];
        let dyte = self.y[1] - self.y[n];
        let dxs = 0.5 * (-self.xp[1] + self.xp[n]);
        let dys = 0.5 * (-self.yp[1] + self.yp[n]);
        self.ante = dxs * dyte - dys * dxte;
        self.aste = dxs * dxte + dys * dyte;
        self.dste = (dxte * dxte + dyte * dyte).sqrt();
        self.sharp = self.dste < 0.0001 * self.chord;
        let (scs, sds) = if self.sharp { (1.0, 0.0) } else { (self.ante / self.dste, self.aste / self.dste) };
        self.sigte = 0.5 * (self.gam[1] - self.gam[n]) * scs;
        self.gamte = -0.5 * (self.gam[1] - self.gam[n]) * sds;
        self.sigte_a = 0.5 * (self.gam_a[1] - self.gam_a[n]) * scs;
        self.gamte_a = -0.5 * (self.gam_a[1] - self.gam_a[n]) * sds;
    }

    /// Integrate surface pressures for CL, CM, CDP (CLCALC) using current gam.
    pub fn clcalc(&mut self) {
        let r = clcalc(
            self.n, &self.x, &self.y, &self.gam, &self.gam_a, self.alfa, self.minf, self.qinf, self.xcmref, self.ycmref,
        );
        self.cl = r.0;
        self.cm = r.1;
        self.cdp = r.2;
        self.cl_alf = r.3;
        self.cl_msq = r.4;
    }

    /// Viscous drag via Squire-Young and friction drag (CDCALC).
    pub fn cdcalc(&mut self) {
        let sa = self.alfa.sin();
        let ca = self.alfa.cos();
        if self.lvisc && self.lblini {
            let nb2 = self.nbl[2];
            let thwake = self.thet[2][nb2];
            let urat = self.uedg[2][nb2] / self.qinf;
            let uewake = self.uedg[2][nb2] * (1.0 - self.tklam) / (1.0 - self.tklam * urat * urat);
            let shwake = self.dstr[2][nb2] / self.thet[2][nb2];
            self.cd = 2.0 * thwake * (uewake / self.qinf).powf(0.5 * (5.0 + shwake));
        } else {
            self.cd = 0.0;
        }
        let mut cdf = 0.0;
        for is in 1..=2 {
            for ibl in 3..=self.iblte[is] {
                let i = self.ipan[is][ibl];
                let im = self.ipan[is][ibl - 1];
                let dx = (self.x[i] - self.x[im]) * ca + (self.y[i] - self.y[im]) * sa;
                cdf += 0.5 * (self.tau[is][ibl] + self.tau[is][ibl - 1]) * dx * 2.0 / (self.qinf * self.qinf);
            }
        }
        self.cdf = cdf;
    }

    fn set_gam_from_alfa(&mut self) {
        let cosa = self.alfa.cos();
        let sina = self.alfa.sin();
        for i in 1..=self.n {
            self.gam[i] = cosa * self.gamu[1][i] + sina * self.gamu[2][i];
            self.gam_a[i] = -sina * self.gamu[1][i] + cosa * self.gamu[2][i];
        }
        self.psio = cosa * self.gamu[1][self.n + 1] + sina * self.gamu[2][self.n + 1];
    }

    /// Converge to specified alpha (SPECAL).
    pub fn specal(&mut self) {
        if !self.lgamu || !self.lqaij {
            self.ggcalc();
        }
        self.cosa = self.alfa.cos();
        self.sina = self.alfa.sin();
        self.set_gam_from_alfa();
        self.tecalc();
        self.qiset();

        let mut clm = 1.0;
        let mut reinf_clm = 0.0;
        let mut minf_clm = self.mrcl(clm, &mut reinf_clm);
        self.comset();
        self.clcalc();
        let mut ok = false;
        for _itcl in 0..20 {
            let msq_clm = 2.0 * self.minf * minf_clm;
            let dclm = (self.cl - clm) / (1.0 - self.cl_msq * msq_clm);
            let clm1 = clm;
            let mut rlx = 1.0;
            for _irlx in 0..12 {
                clm = clm1 + rlx * dclm;
                minf_clm = self.mrcl(clm, &mut reinf_clm);
                if self.matyp == 1 || self.minf == 0.0 || minf_clm != 0.0 {
                    break;
                }
                rlx *= 0.5;
            }
            self.comset();
            self.clcalc();
            if dclm.abs() <= 1.0e-6 {
                ok = true;
                break;
            }
        }
        if !ok && self.verbose {
            crate::outln!(" SPECAL:  Minf convergence failed");
        }
        let cl = self.cl;
        let mut rcl = 0.0;
        self.minf_cl = self.mrcl(cl, &mut rcl);
        self.reinf_cl = rcl;
        self.comset();
        self.clcalc();
        let n = self.n;
        let nw = self.nw;
        if self.lvisc {
            cpcalc(n + nw, &self.qvis, self.qinf, self.minf, &mut self.cpv);
            cpcalc(n + nw, &self.qinv, self.qinf, self.minf, &mut self.cpi);
        } else {
            cpcalc(n, &self.qinv, self.qinf, self.minf, &mut self.cpi);
        }
    }

    /// Converge to specified inviscid CL (SPECCL).
    pub fn speccl(&mut self) {
        if !self.lgamu || !self.lqaij {
            self.ggcalc();
        }
        let mut rcl = 0.0;
        self.minf_cl = self.mrcl(self.clspec, &mut rcl);
        self.reinf_cl = rcl;
        self.comset();
        self.set_gam_from_alfa();
        self.clcalc();
        let mut ok = false;
        for _ in 0..20 {
            let dalfa = (self.clspec - self.cl) / self.cl_alf;
            self.alfa += dalfa;
            self.set_gam_from_alfa();
            self.clcalc();
            if dalfa.abs() <= 1.0e-6 {
                ok = true;
                break;
            }
        }
        if !ok && self.verbose {
            crate::outln!(" SPECCL:  CL convergence failed");
        }
        self.tecalc();
        self.qiset();
        let n = self.n;
        let nw = self.nw;
        if self.lvisc {
            cpcalc(n + nw, &self.qvis, self.qinf, self.minf, &mut self.cpv);
            cpcalc(n + nw, &self.qinv, self.qinf, self.minf, &mut self.cpi);
        } else {
            cpcalc(n, &self.qinv, self.qinf, self.minf, &mut self.cpi);
        }
    }

    /// Converge viscous operating point (VISCAL). Returns true if converged.
    pub fn viscal(&mut self, niter: usize) -> bool {
        const EPS1: f64 = 1.0e-4;
        if !self.lwake {
            self.xywake();
        }
        self.qwcalc();
        self.qiset();
        if !self.lipan {
            if self.lblini {
                self.gamqv();
            }
            self.stfind();
            self.iblpan();
            self.xicalc();
            self.iblsys();
        }
        self.uicalc();
        if !self.lblini {
            for is in 1..=2 {
                for ibl in 1..=self.nbl[is] {
                    self.uedg[is][ibl] = self.uinv[is][ibl];
                }
            }
        }
        let n = self.n;
        let nw = self.nw;
        if self.lvconv {
            self.qvfue();
            if self.lvisc {
                cpcalc(n + nw, &self.qvis, self.qinf, self.minf, &mut self.cpv);
                cpcalc(n + nw, &self.qinv, self.qinf, self.minf, &mut self.cpi);
            } else {
                cpcalc(n, &self.qinv, self.qinf, self.minf, &mut self.cpi);
            }
            self.gamqv();
            self.clcalc();
            self.cdcalc();
        }
        if !self.lwdij || !self.ladij {
            self.qdcalc();
        }
        if self.verbose {
            crate::outln!("\n Solving BL system ...");
        }
        let mut conv = false;
        self.viscal_iters = 0;
        for iter in 1..=niter {
            self.viscal_iters = iter;
            self.setbl();
            self.blsolv();
            self.update();
            if self.lalfa {
                let cl = self.cl;
                let mut rcl = 0.0;
                self.minf_cl = self.mrcl(cl, &mut rcl);
                self.reinf_cl = rcl;
                self.comset();
            } else {
                self.qiset();
                self.uicalc();
            }
            self.qvfue();
            self.gamqv();
            self.stmove();
            self.clcalc();
            self.cdcalc();
            if self.verbose {
                if self.rlx < 1.0 {
                    crate::outln!(
                        "\n {:3}   rms: {}   max: {}   {} at {:4}{:3}   RLX:{:6.3}",
                        iter,
                        crate::fmt::e_fmt(self.rmsbl, 10, 4),
                        crate::fmt::e_fmt(self.rmxbl, 10, 4),
                        self.vmxbl,
                        self.imxbl,
                        self.ismxbl,
                        self.rlx
                    );
                } else {
                    crate::outln!(
                        "\n {:3}   rms: {}   max: {}   {} at {:4}{:3}",
                        iter,
                        crate::fmt::e_fmt(self.rmsbl, 10, 4),
                        crate::fmt::e_fmt(self.rmxbl, 10, 4),
                        self.vmxbl,
                        self.imxbl,
                        self.ismxbl
                    );
                }
            }
            self.cdp = self.cd - self.cdf;
            if self.verbose {
                crate::outln!("       a ={:7.3}      CL ={:8.4}", self.alfa / DTOR, self.cl);
                crate::outln!("      Cm ={:8.4}     CD ={:9.5}   =>   CDf ={:9.5}    CDp ={:9.5}", self.cm, self.cd, self.cdf, self.cdp);
            }
            if self.rmsbl < EPS1 {
                self.lvconv = true;
                self.avisc = self.alfa;
                self.mvisc = self.minf;
                conv = true;
                break;
            }
        }
        if !conv && self.verbose {
            crate::outln!(" VISCAL:  Convergence failed");
        }
        cpcalc(n + nw, &self.qinv, self.qinf, self.minf, &mut self.cpi);
        cpcalc(n + nw, &self.qvis, self.qinf, self.minf, &mut self.cpv);
        conv
    }

    /// Minimum Cp (FCPMIN).
    pub fn fcpmin(&mut self) {
        let n = self.n;
        self.xcpmni = self.x[1];
        self.xcpmnv = self.x[1];
        self.cpmni = self.cpi[1];
        self.cpmnv = self.cpv[1];
        for i in 2..=n + self.nw {
            if self.cpi[i] < self.cpmni {
                self.xcpmni = self.x[i];
                self.cpmni = self.cpi[i];
            }
            if self.cpv[i] < self.cpmnv {
                self.xcpmnv = self.x[i];
                self.cpmnv = self.cpv[i];
            }
        }
        if self.lvisc {
            self.cpmn = self.cpmnv;
        } else {
            self.cpmn = self.cpmni;
            self.cpmnv = self.cpmni;
            self.xcpmnv = self.xcpmni;
        }
    }
}

/// Compressible Cp from speed (CPCALC). Returns true if density went negative.
pub fn cpcalc(n: usize, q: &[f64], qinf: f64, minf: f64, cp: &mut [f64]) -> bool {
    let beta = (1.0 - minf * minf).sqrt();
    let bfac = 0.5 * minf * minf / (1.0 + beta);
    let mut denneg = false;
    for i in 1..=n {
        let cpinc = 1.0 - (q[i] / qinf) * (q[i] / qinf);
        let den = beta + bfac * cpinc;
        cp[i] = cpinc / den;
        if den <= 0.0 {
            denneg = true;
        }
    }
    denneg
}

/// CLCALC: returns (cl, cm, cdp, cl_alf, cl_msq).
#[allow(clippy::too_many_arguments)]
pub fn clcalc(
    n: usize,
    x: &[f64],
    y: &[f64],
    gam: &[f64],
    gam_a: &[f64],
    alfa: f64,
    minf: f64,
    qinf: f64,
    xref: f64,
    yref: f64,
) -> (f64, f64, f64, f64, f64) {
    let sa = alfa.sin();
    let ca = alfa.cos();
    let beta = (1.0 - minf * minf).sqrt();
    let beta_msq = -0.5 / beta;
    let bfac = 0.5 * minf * minf / (1.0 + beta);
    let bfac_msq = 0.5 / (1.0 + beta) - bfac / (1.0 + beta) * beta_msq;
    let mut cl = 0.0;
    let mut cm = 0.0;
    let mut cdp = 0.0;
    let mut cl_alf = 0.0;
    let mut cl_msq = 0.0;
    let i = 1;
    let cginc = 1.0 - (gam[i] / qinf) * (gam[i] / qinf);
    let mut cpg1 = cginc / (beta + bfac * cginc);
    let mut cpg1_msq = -cpg1 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
    let cpi_gam = -2.0 * gam[i] / (qinf * qinf);
    let cpc_cpi = (1.0 - bfac * cpg1) / (beta + bfac * cginc);
    let mut cpg1_alf = cpc_cpi * cpi_gam * gam_a[i];
    for i in 1..=n {
        let ip = if i == n { 1 } else { i + 1 };
        let cginc = 1.0 - (gam[ip] / qinf) * (gam[ip] / qinf);
        let cpg2 = cginc / (beta + bfac * cginc);
        let cpg2_msq = -cpg2 / (beta + bfac * cginc) * (beta_msq + bfac_msq * cginc);
        let cpi_gam = -2.0 * gam[ip] / (qinf * qinf);
        let cpc_cpi = (1.0 - bfac * cpg2) / (beta + bfac * cginc);
        let cpg2_alf = cpc_cpi * cpi_gam * gam_a[ip];

        let dx = (x[ip] - x[i]) * ca + (y[ip] - y[i]) * sa;
        let dy = (y[ip] - y[i]) * ca - (x[ip] - x[i]) * sa;
        let dg = cpg2 - cpg1;
        let ax = (0.5 * (x[ip] + x[i]) - xref) * ca + (0.5 * (y[ip] + y[i]) - yref) * sa;
        let ay = (0.5 * (y[ip] + y[i]) - yref) * ca - (0.5 * (x[ip] + x[i]) - xref) * sa;
        let ag = 0.5 * (cpg2 + cpg1);
        let dx_alf = -(x[ip] - x[i]) * sa + (y[ip] - y[i]) * ca;
        let ag_alf = 0.5 * (cpg2_alf + cpg1_alf);
        let ag_msq = 0.5 * (cpg2_msq + cpg1_msq);
        cl += dx * ag;
        cdp -= dy * ag;
        cm -= dx * (ag * ax + dg * dx / 12.0) + dy * (ag * ay + dg * dy / 12.0);
        cl_alf += dx * ag_alf + ag * dx_alf;
        cl_msq += dx * ag_msq;
        cpg1 = cpg2;
        cpg1_alf = cpg2_alf;
        cpg1_msq = cpg2_msq;
    }
    (cl, cm, cdp, cl_alf, cl_msq)
}

/// Wake geometry and wake source-influence data for one angle of attack
/// (everything XYWAKE + QWCALC + QDCALC produce that depends on alpha).
/// Used to precompute the next sweep point on a helper thread.
#[derive(Clone, Debug)]
pub struct WakePack {
    pub alfa: f64,
    pub nw: usize,
    pub n: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub s: Vec<f64>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub apanel: Vec<f64>,
    pub q1: Vec<f64>,
    pub q2: Vec<f64>,
    pub dij: Mat,
}

impl XFoil {
    /// Extract the wake data after XYWAKE/QWCALC/QDCALC have run.
    pub fn wake_pack(&self) -> WakePack {
        let (n, nw) = (self.n, self.nw);
        let r = n + 1..=n + nw;
        WakePack {
            alfa: self.awake,
            nw,
            n,
            x: self.x[r.clone()].to_vec(),
            y: self.y[r.clone()].to_vec(),
            s: self.s[r.clone()].to_vec(),
            nx: self.nx[r.clone()].to_vec(),
            ny: self.ny[r.clone()].to_vec(),
            apanel: self.apanel[r.clone()].to_vec(),
            q1: self.qinvu[1][r.clone()].to_vec(),
            q2: self.qinvu[2][r].to_vec(),
            dij: self.dij.clone(),
        }
    }

    /// Install precomputed wake data (bit-identical to computing it here).
    pub fn install_wake(&mut self, p: WakePack) {
        let (n, nw) = (p.n, p.nw);
        debug_assert_eq!(n, self.n);
        self.nw = nw;
        self.x[n + 1..=n + nw].copy_from_slice(&p.x);
        self.y[n + 1..=n + nw].copy_from_slice(&p.y);
        self.s[n + 1..=n + nw].copy_from_slice(&p.s);
        self.nx[n + 1..=n + nw].copy_from_slice(&p.nx);
        self.ny[n + 1..=n + nw].copy_from_slice(&p.ny);
        self.apanel[n + 1..=n + nw].copy_from_slice(&p.apanel);
        self.qinvu[1][n + 1..=n + nw].copy_from_slice(&p.q1);
        self.qinvu[2][n + 1..=n + nw].copy_from_slice(&p.q2);
        self.dij = p.dij;
        self.xte = 0.5 * (self.x[1] + self.x[n]);
        self.yte = 0.5 * (self.y[1] + self.y[n]);
        self.awake = p.alfa;
        self.lwake = true;
        self.lwdij = true;
        self.ladij = true;
        self.qw_valid = true;
    }
}
