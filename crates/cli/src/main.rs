//! `foiley` - headless command-line interface.
//!
//! Run `foiley --help` or `foiley <command> --help` for every option.
//! Without a command and with a script piped on stdin, Foiley behaves
//! exactly like XFoil (existing XFoil batch scripts run unchanged).

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod legacy;
mod progress_ui;

use clap::{Args, Parser, Subcommand, ValueEnum};
use foiley_app::{auto, JobKind, JobManager, JobOutput, Project};
use foiley_core::atmosphere::{isa, Fluid};
use foiley_core::Vec3;
use foiley_foil::analysis::PolarType;
use foiley_foil::batch::BatchCase;
use foiley_foil::foil::Paneling;
use foiley_foil::{AnalysisSettings, Foil, FoilPolar, Sweep};
use foiley_plane::analysis::{Method, PolarKind, Range, RefDims};
use foiley_plane::geom::Inertia;
use foiley_plane::{Plane, PlanePolar, PlanePolarDef};
use progress_ui::Ui;
use std::io::IsTerminal;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(name = "foiley", version, about = "Foiley: airfoil, wing and plane analysis (XFoil 2D, LLT/VLM/panel 3D, stability)", long_about = None)]
struct Cli {
    /// compute threads (0 = all cores but one)
    #[arg(long, global = true, default_value_t = 0)]
    threads: usize,
    /// only print warnings and errors on stderr
    #[arg(short, long, global = true)]
    quiet: bool,
    /// print per-point solver details
    #[arg(short, long, global = true)]
    verbose: bool,
    /// machine-readable progress events (JSON lines on stderr)
    #[arg(long, global = true)]
    json_progress: bool,
    #[command(subcommand)]
    cmd: Option<Cmd>,
}

#[derive(Subcommand)]
enum Cmd {
    /// XFoil-compatible menu mode (reads XFoil commands from stdin)
    Xfoil {
        /// airfoil file to load first
        file: Option<String>,
    },
    /// 2D airfoil tools and analyses
    #[command(subcommand)]
    Foil(FoilCmd),
    /// 3D wing / plane analyses
    #[command(subcommand)]
    Plane(PlaneCmd),
    /// Run a JSON job script (see `foiley run --help`)
    Run {
        /// job script: {"project": .., "steps": [..], "save": ..}
        script: PathBuf,
    },
    /// Show (and optionally convert) a project: .foiley, xflr5 .xfl, flow5 .fl5 or .plr polars
    Project {
        file: PathBuf,
        /// save as a Foiley project
        #[arg(long)]
        save: Option<PathBuf>,
    },
}

#[derive(Clone, Copy, ValueEnum, PartialEq)]
enum Format {
    Table,
    Csv,
    Json,
    /// XFoil polar text format (2D polars only)
    Xfoil,
}

#[derive(Args, Clone)]
struct Out {
    /// output file (default: stdout)
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// output format
    #[arg(short, long, value_enum, default_value_t = Format::Table)]
    format: Format,
}

#[derive(Args, Clone)]
struct FoilSrc {
    /// airfoil: a .dat file, "NACA 2412" / "naca4412", or a foil name in --project
    foil: String,
    /// project file to look the foil up in
    #[arg(long)]
    project: Option<PathBuf>,
}

#[derive(Args, Clone)]
struct XfoilOpts {
    /// Reynolds number (type 2: Re*sqrt(CL); type 3: Re*CL)
    #[arg(long, default_value_t = 1e6)]
    re: f64,
    /// Mach number
    #[arg(long, default_value_t = 0.0)]
    mach: f64,
    /// e^N transition criterion (9 = wind tunnel, 11-12 = quiet free air, 4-6 = turbulent)
    #[arg(long, default_value_t = 9.0)]
    ncrit: f64,
    /// forced transition x/c "top,bottom" (1 = free transition)
    #[arg(long, default_value = "1,1")]
    xtr: String,
    /// inviscid analysis (no boundary layer)
    #[arg(long)]
    inviscid: bool,
    /// Newton iteration limit per point
    #[arg(long, default_value_t = 100)]
    iter: usize,
    /// number of panel nodes (XFoil PANGEN, max 354)
    #[arg(long, default_value_t = 160)]
    panels: usize,
    /// do not try to recover unconverged points
    #[arg(long)]
    no_recover: bool,
    /// march in the given order instead of outwards from zero lift
    #[arg(long)]
    no_from_zero: bool,
    /// single-threaded sweep (bit-identical to XFoil's own sequence)
    #[arg(long)]
    sequential: bool,
}

impl XfoilOpts {
    fn settings(&self) -> Result<AnalysisSettings, String> {
        let x = parse_list(&self.xtr)?;
        if x.len() != 2 {
            return Err("--xtr needs two values: top,bottom".into());
        }
        Ok(AnalysisSettings {
            re: self.re,
            mach: self.mach,
            ncrit: self.ncrit,
            xtr_top: x[0],
            xtr_bot: x[1],
            viscous: !self.inviscid,
            max_iter: self.iter,
            paneling: Paneling { panels: self.panels, ..Default::default() },
            recover: !self.no_recover,
            from_zero: !self.no_from_zero,
            parallel: !self.sequential,
            ..Default::default()
        })
    }
}

#[derive(Subcommand)]
enum FoilCmd {
    /// Geometric properties of one or more foils
    Info {
        foils: Vec<String>,
        #[arg(long)]
        project: Option<PathBuf>,
    },
    /// Polar: sweep alpha, CL, Re (type 4) or flap angle (type 6)
    Polar {
        #[command(flatten)]
        src: FoilSrc,
        #[command(flatten)]
        x: XfoilOpts,
        /// polar type: 1 fixed speed, 2 fixed lift, 3 rubber chord, 4 fixed alpha (sweep Re), 6 flap sweep
        #[arg(long = "type", default_value_t = 1)]
        ptype: u8,
        /// alpha sweep "from:to:step" (degrees) or list "0,2,4"
        #[arg(long, allow_hyphen_values = true)]
        alpha: Option<String>,
        /// CL sweep "from:to:step"
        #[arg(long, allow_hyphen_values = true)]
        cl: Option<String>,
        /// Re sweep "from:to:step" (type 4)
        #[arg(long)]
        re_sweep: Option<String>,
        /// flap-angle sweep "from:to:step" (type 6)
        #[arg(long, allow_hyphen_values = true)]
        flap_sweep: Option<String>,
        /// fixed alpha for type 4 / 6
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        at_alpha: f64,
        /// flap hinge "x/c,y/t" for type 6
        #[arg(long, default_value = "0.75,0.5")]
        hinge: String,
        #[command(flatten)]
        out: Out,
    },
    /// Single operating point with Cp and boundary-layer distributions
    Point {
        #[command(flatten)]
        src: FoilSrc,
        #[command(flatten)]
        x: XfoilOpts,
        #[arg(long, allow_hyphen_values = true)]
        alpha: Option<f64>,
        #[arg(long, allow_hyphen_values = true)]
        cl: Option<f64>,
        #[command(flatten)]
        out: Out,
    },
    /// Polar mesh: foils x Re x Mach x Ncrit, in parallel
    Batch {
        /// foils (files, NACA codes or project names)
        foils: Vec<String>,
        #[arg(long)]
        project: Option<PathBuf>,
        /// Reynolds numbers: list "1e5,2e5" or geometric "from:to:count"
        #[arg(long, default_value = "1e5:1e6:6")]
        re: String,
        #[arg(long, default_value = "0")]
        mach: String,
        #[arg(long, default_value = "9")]
        ncrit: String,
        #[arg(long, default_value = "-6:16:0.5", allow_hyphen_values = true)]
        alpha: String,
        #[command(flatten)]
        x: XfoilOpts,
        /// directory for one CSV per polar
        #[arg(long)]
        out_dir: Option<PathBuf>,
        /// add the polars to this project file (created if missing)
        #[arg(long)]
        save: Option<PathBuf>,
    },
    /// Shape optimisation (MOPSO particle swarm or genetic algorithm, Hicks-Henne bumps)
    Optimize {
        #[command(flatten)]
        src: FoilSrc,
        #[command(flatten)]
        x: XfoilOpts,
        /// design angle of attack, degrees
        #[arg(long, default_value_t = 2.0, allow_hyphen_values = true)]
        alpha: f64,
        /// objectives "cl/cd:max,cd:min,cl:0.6[:tol],cm:0,cpmin:max"
        #[arg(long, default_value = "cl/cd:max")]
        objectives: String,
        /// mopso or ga
        #[arg(long, default_value = "mopso")]
        algorithm: String,
        #[arg(long, default_value_t = 12)]
        population: usize,
        #[arg(long, default_value_t = 25)]
        iterations: usize,
        /// Hicks-Henne bumps per surface
        #[arg(long, default_value_t = 3)]
        bumps: usize,
        /// max bump amplitude, fraction of chord
        #[arg(long, default_value_t = 0.01)]
        amplitude: f64,
        /// also optimise a TE flap angle "min,max" (deg)
        #[arg(long, allow_hyphen_values = true)]
        flap: Option<String>,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        /// output .dat of the best design
        #[arg(short, long)]
        out: Option<PathBuf>,
        /// CSV of the Pareto set
        #[arg(long)]
        pareto: Option<PathBuf>,
    },
    /// Inverse design: modify the surface-speed distribution, solve for the foil
    Inverse {
        #[command(flatten)]
        src: FoilSrc,
        /// design angle of attack, degrees
        #[arg(long, default_value_t = 2.0, allow_hyphen_values = true)]
        alpha: f64,
        /// mixed (QDES, one segment) or full (MDES, whole surface)
        #[arg(long, default_value = "mixed")]
        mode: String,
        /// segment "upper|lower,x_from,x_to"
        #[arg(long)]
        segment: String,
        /// multiply the speed in the segment by this factor
        #[arg(long)]
        scale: Option<f64>,
        /// add a bump of this amplitude (Q/Vinf) in the segment
        #[arg(long, allow_hyphen_values = true)]
        bump: Option<f64>,
        /// smooth the speed in the segment
        #[arg(long)]
        smooth: bool,
        #[arg(long, default_value_t = 10)]
        iterations: usize,
        /// write the speed distributions (before/after) as CSV
        #[arg(long)]
        speed_csv: Option<PathBuf>,
        #[arg(short, long)]
        out: Option<PathBuf>,
    },
    /// Modify a foil (XFoil GDES operations) and write it as .dat
    Design {
        #[command(flatten)]
        src: FoilSrc,
        /// trailing-edge flap "x/c,y/t,deg" (+ = down)
        #[arg(long, allow_hyphen_values = true)]
        flap: Option<String>,
        /// trailing-edge gap "gap[,blend]" (chord fractions)
        #[arg(long)]
        te_gap: Option<String>,
        /// leading-edge radius scale "factor[,blend]"
        #[arg(long)]
        le_radius: Option<String>,
        /// set max thickness (chord fraction)
        #[arg(long)]
        thickness: Option<f64>,
        /// set max camber (chord fraction)
        #[arg(long, allow_hyphen_values = true)]
        camber: Option<f64>,
        /// move max thickness / camber to "xt,xc"
        #[arg(long)]
        high_points: Option<String>,
        /// repanel with N nodes
        #[arg(long)]
        repanel: Option<usize>,
        /// de-rotate and scale to unit chord
        #[arg(long)]
        normalize: bool,
        /// interpolate with another foil "FOIL,fraction"
        #[arg(long)]
        mix: Option<String>,
        /// new name
        #[arg(long)]
        name: Option<String>,
        /// output .dat (default: stdout)
        #[arg(short, long)]
        out: Option<PathBuf>,
    },
}

#[derive(Args, Clone)]
struct PlaneSrc {
    /// plane: xflr5/flow5 .xml file, or a .foiley project
    plane: PathBuf,
    /// plane name inside a project (default: first plane)
    #[arg(long)]
    name: Option<String>,
    /// extra directory to look for foil .dat files
    #[arg(long)]
    foils: Option<PathBuf>,
}

#[derive(Subcommand)]
enum PlaneCmd {
    /// Geometry, reference dimensions and mass properties
    Info {
        #[command(flatten)]
        src: PlaneSrc,
    },
    /// Export: AVL geometry and mass files, STL mesh, xflr5 XML
    Export {
        #[command(flatten)]
        src: PlaneSrc,
        #[arg(long)]
        avl: Option<PathBuf>,
        #[arg(long)]
        mass: Option<PathBuf>,
        /// panel mesh as binary STL
        #[arg(long)]
        stl: Option<PathBuf>,
        /// thin (camber-surface) instead of thick STL mesh
        #[arg(long)]
        thin: bool,
        #[arg(long)]
        xml: Option<PathBuf>,
    },
    /// Replace the plane's fuselage with an STL / STEP surface and save the project
    ImportFuselage {
        #[command(flatten)]
        src: PlaneSrc,
        /// .stl, .step or .stp file
        #[arg(long)]
        file: PathBuf,
        /// file units: auto, m, mm, cm, in, ft
        #[arg(long, default_value = "auto")]
        units: String,
        /// target triangle count for STEP tessellation
        #[arg(long, default_value_t = 3000)]
        triangles: usize,
        /// structural mass of the fuselage, kg
        #[arg(long, default_value_t = 0.0)]
        mass: f64,
        /// output project
        #[arg(long)]
        save: PathBuf,
    },
    /// Run a 3D polar (types 1-8)
    Analyze {
        #[command(flatten)]
        src: PlaneSrc,
        /// t1 fixed speed, t2 fixed lift, t3 glide, t4 fixed alpha, t5 sideslip, t7 stability, t8 list
        #[arg(long = "type", default_value = "t1")]
        ptype: String,
        /// llt, vlm1, vlm2, panel (quads), tripanel (uniform triangles), trilinear (Galerkin)
        #[arg(long, default_value = "vlm2")]
        method: String,
        /// VLM free-wake relaxation iterations (0 = rigid wake)
        #[arg(long, default_value_t = 0)]
        free_wake: usize,
        /// free wake: segments per trailing line
        #[arg(long, default_value_t = 24)]
        wake_nseg: usize,
        /// free wake: under-relaxation of the node update, in (0, 1]
        #[arg(long, default_value_t = 0.5)]
        wake_relax: f64,
        /// free wake: convergence tolerance (node displacement / c_ref and relative circulation change)
        #[arg(long, default_value_t = 1e-3)]
        wake_tol: f64,
        /// T7: forced response "step", "pulse:SECONDS", "doublet:SECONDS" or "sine:HZ" (CSV to --response-out)
        #[arg(long)]
        response: Option<String>,
        /// amplitude of the control input, units of the control parameter
        #[arg(long, default_value_t = 0.2)]
        response_amp: f64,
        #[arg(long, default_value_t = 10.0)]
        response_time: f64,
        #[arg(long)]
        response_out: Option<PathBuf>,
        /// sweep "from:to:step" (alpha for t1-t3, speed for t4, beta for t5, control for t7)
        #[arg(long, default_value = "-4:10:1", allow_hyphen_values = true)]
        range: String,
        /// speed for t1/t5, m/s
        #[arg(long, default_value_t = 15.0)]
        speed: f64,
        /// alpha for t4/t5, degrees
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        alpha: f64,
        /// sideslip for t1-t4, degrees
        #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
        beta: f64,
        /// t8 points "a,b,V;a,b,V;..."
        #[arg(long)]
        list: Option<String>,
        /// skip viscous drag
        #[arg(long)]
        inviscid: bool,
        /// altitude for ISA air properties, m
        #[arg(long, default_value_t = 0.0)]
        altitude: f64,
        /// run in water (hydrofoils)
        #[arg(long)]
        water: bool,
        /// mass override, kg
        #[arg(long)]
        mass: Option<f64>,
        /// CoG override "x,y,z", m
        #[arg(long)]
        cog: Option<String>,
        /// ground height below the origin, m
        #[arg(long)]
        ground: Option<f64>,
        /// stability derivatives and modes for every point (needs mass)
        #[arg(long)]
        derivatives: bool,
        /// fail points whose viscous drag is clamped to the 2D polar envelope
        #[arg(long)]
        strict_envelope: bool,
        /// flap preset "wing,flap,deg" (repeatable)
        #[arg(long = "flap", allow_hyphen_values = true)]
        flaps: Vec<String>,
        /// T7 gain "wing,flap,deg_per_unit" (repeatable)
        #[arg(long = "gain", allow_hyphen_values = true)]
        gains: Vec<String>,
        /// reference dimensions: projected, planform or "area,span,chord"
        #[arg(long, default_value = "projected")]
        reference: String,
        /// number of Re values for automatic 2D polar meshes (0 = never compute missing polars)
        #[arg(long, default_value_t = 6)]
        auto_polars: usize,
        /// LLT stations
        #[arg(long, default_value_t = 40)]
        llt_stations: usize,
        /// span distribution output (CSV, one block per point)
        #[arg(long)]
        strips: Option<PathBuf>,
        /// save plane, foils, 2D and 3D polars into this project
        #[arg(long)]
        save: Option<PathBuf>,
        #[command(flatten)]
        out: Out,
    },
}

fn parse_list(s: &str) -> Result<Vec<f64>, String> {
    s.split(|c| c == ',' || c == ';')
        .map(|t| t.trim())
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<f64>().map_err(|_| format!("'{t}' is not a number")))
        .collect()
}

/// "a:b:s" -> (a, b, s); "x" -> (x, x, 1)
fn parse_range(s: &str) -> Result<(f64, f64, f64), String> {
    let v: Vec<&str> = s.split(':').collect();
    let num = |t: &str| t.trim().parse::<f64>().map_err(|_| format!("'{t}' is not a number in range '{s}'"));
    match v.len() {
        1 => {
            let x = num(v[0])?;
            Ok((x, x, 1.0))
        }
        3 => Ok((num(v[0])?, num(v[1])?, num(v[2])?)),
        _ => Err(format!("range '{s}' must be 'from:to:step' or a single value")),
    }
}

fn load_project(p: &Option<PathBuf>) -> Result<Project, String> {
    match p {
        Some(path) => Project::load(path),
        None => Ok(Project::new("cli")),
    }
}

/// Resolve a foil spec: file, NACA code or project name.
fn resolve_foil(spec: &str, proj: &Project) -> Result<Foil, String> {
    if let Some(f) = proj.foil(spec) {
        return Ok(f.clone());
    }
    let path = Path::new(spec);
    if path.exists() {
        return Foil::load(path).map_err(|e| e.to_string());
    }
    let l = spec.to_lowercase().replace([' ', '-', '_'], "");
    if let Some(code) = l.strip_prefix("naca") {
        return Foil::naca(code).map_err(|e| e.to_string());
    }
    if l.chars().all(|c| c.is_ascii_digit()) && (l.len() == 4 || l.len() == 5) {
        return Foil::naca(&l).map_err(|e| e.to_string());
    }
    Err(format!("'{spec}' is not a file, a NACA designation, or a foil in the project"))
}

fn write_out(out: &Out, text: &str) -> Result<(), String> {
    match &out.out {
        Some(p) => std::fs::write(p, text).map_err(|e| format!("{}: {e}", p.display())),
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

fn foil_polar_table(p: &FoilPolar) -> String {
    let mut s = format!("# {} | {} | Re {:.0} M {:.2} Ncrit {:.1} xtr {:.2}/{:.2}\n", p.foil, p.name, p.re, p.mach, p.ncrit, p.xtr_top, p.xtr_bot);
    s.push_str("  alpha      CL        CD       CDp       Cm    Xtr_top  Xtr_bot    L/D       Re   flag\n");
    for q in &p.points {
        s.push_str(&format!(
            "{:7.2} {:8.4} {:9.5} {:9.5} {:8.4} {:8.4} {:8.4} {:7.2} {:9.0} {}\n",
            if p.polar_type == 6 { q.flap } else { q.alpha },
            q.cl,
            q.cd,
            q.cdp,
            q.cm,
            q.xtr_top,
            q.xtr_bot,
            q.ld(),
            q.re,
            if q.recovered { "recovered" } else { "" }
        ));
    }
    if !p.failed.is_empty() {
        s.push_str(&format!("# not converged: {:?}\n", p.failed));
    }
    if let (Some(a0), Some(cla)) = (p.alpha0(), p.cl_alpha()) {
        s.push_str(&format!("# alpha0 = {a0:.3} deg, CLa = {cla:.4} /deg"));
        if let Some(m) = p.cl_max() {
            s.push_str(&format!(", CLmax = {:.4} at {:.2} deg", m.cl, m.alpha));
        }
        if let Some(m) = p.ld_max() {
            s.push_str(&format!(", L/Dmax = {:.1} at {:.2} deg", m.ld(), m.alpha));
        }
        s.push('\n');
    }
    s
}

fn plane_polar_table(p: &PlanePolar) -> String {
    let (s, b, c) = p.ref_dims;
    let mut t = format!("# {} | {} | Sref {s:.4} m2, bref {b:.3} m, cref {c:.4} m, {} unknowns\n", p.plane, p.name, p.unknowns);
    t.push_str("  alpha   beta      V    ctrl       CL        CD      CDi      CDv        Cm        Cn    L/D     Vz\n");
    for o in &p.points {
        t.push_str(&format!(
            "{:7.2} {:6.2} {:6.2} {:7.3} {:8.4} {:9.5} {:8.5} {:8.5} {:9.5} {:>9} {:6.2} {:6.3}\n",
            o.alpha, o.beta, o.velocity, o.ctrl, o.cl, o.cd, o.cdi, o.cdv, o.cm, o.cn.map(|v| format!("{v:.5}")).unwrap_or_else(|| "-".into()), o.ld(), o.vz()
        ));
    }
    for o in &p.points {
        if let Some(s) = &o.stability {
            t.push_str(&format!(
                "# ctrl {:.3}: trim alpha {:.3} deg, V {:.3} m/s | static margin {:.1}% | x_NP {:.4} m | Cma {:.3} Cmq {:.2} Clp {:.3} Cnb {:.4}\n",
                o.ctrl,
                s.alpha,
                s.u0,
                s.static_margin * 100.0,
                s.x_np,
                s.nondim.cma,
                s.nondim.cmq,
                s.nondim.clp,
                s.nondim.cnb
            ));
            if let Some(c) = &s.control {
                t.push_str(&format!("#   control derivatives per unit: X {:+.4} Y {:+.4} Z {:+.4} N | L {:+.4} M {:+.4} N {:+.4} N.m\n", c[0], c[1], c[2], c[3], c[4], c[5]));
            }
            for m in &s.modes {
                t.push_str(&format!(
                    "#   {:16} lambda = {:+9.4} {:+9.4}i  f = {:7.4} Hz  zeta = {:+.3}  {} {:.2} s\n",
                    m.name,
                    m.eigenvalue.re,
                    m.eigenvalue.im,
                    m.freq_hz,
                    m.damping,
                    if m.stable { "t_half" } else { "t_double" },
                    m.t_half_or_double
                ));
            }
        }
    }
    let mut warns: Vec<&String> = p.points.iter().flat_map(|o| o.warnings.iter()).collect();
    warns.sort();
    warns.dedup();
    for w in warns {
        t.push_str(&format!("# warning: {w}\n"));
    }
    for f in &p.failed {
        t.push_str(&format!("# failed: {f}\n"));
    }
    t
}

fn strips_csv(p: &PlanePolar) -> String {
    let mut s = String::from("alpha,V,ctrl,wing,span_pos,chord,Cl,Cdi,Cdv,ai_deg,Re,Xtr_top,Xtr_bot,gamma,lift_per_span,bending,out_of_polar,envelope\n");
    for o in &p.points {
        for r in &o.strips {
            s.push_str(&format!(
                "{:.3},{:.3},{:.3},{},{:.5},{:.5},{:.5},{:.6},{:.6},{:.4},{:.0},{:.4},{:.4},{:.5},{:.4},{:.4},{},{}\n",
                o.alpha, o.velocity, o.ctrl, r.wing, r.span_pos, r.chord, r.cl, r.cdi, r.cdv, r.ai, r.re, r.xtr_top, r.xtr_bot, r.gamma, r.lift_per_span, r.bending, r.out_of_polar, foiley_plane::foildb::Envelope::from_u8(r.envelope).label()
            ));
        }
    }
    s
}

/// Load a plane (xml or project) and the foils it needs.
fn load_plane(src: &PlaneSrc) -> Result<(Project, Plane), String> {
    let ext = src.plane.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
    let mut proj = if ext == "foiley" || ext == "json" { Project::load(&src.plane)? } else { Project::new("cli") };
    if ext == "xml" {
        let msg = proj.import(&src.plane)?;
        eprintln!("  imported {msg}");
    }
    let plane = match &src.name {
        Some(n) => proj.planes.iter().find(|p| &p.name == n).cloned().ok_or_else(|| format!("no plane named '{n}'"))?,
        None => proj.planes.first().cloned().ok_or("the file contains no plane")?,
    };
    // missing foils: look in --foils dir, then NACA designations
    for name in plane.foil_names() {
        if proj.foil(&name).is_some() {
            continue;
        }
        let mut found = None;
        if let Some(d) = &src.foils {
            let p = d.join(format!("{name}.dat"));
            if let Ok(f) = Foil::load(&p) {
                found = Some(f);
            }
        }
        if found.is_none() {
            let l = name.to_lowercase().replace([' ', '-', '_'], "");
            if let Some(code) = l.strip_prefix("naca") {
                found = Foil::naca(code).ok();
            }
        }
        match found {
            Some(mut f) => {
                f.name = name.clone();
                proj.foils.push(f);
            }
            None => return Err(format!("foil '{name}' used by the plane was not found (give --foils DIR with {name}.dat)")),
        }
    }
    Ok((proj, plane))
}

fn run_job(mgr: &mut JobManager, ui: &mut Ui, title: &str, kind: JobKind) -> Result<JobOutput, String> {
    let job = mgr.submit(title, kind);
    ui.wait(mgr, &job)?;
    job.take_output().ok_or_else(|| "job produced no output".into())
}

fn ensure_polars(mgr: &mut JobManager, ui: &mut Ui, proj: &mut Project, plane: &Plane, def: &PlanePolarDef, count: usize) -> Result<(), String> {
    if !(def.viscous || def.method == Method::Llt) || count == 0 {
        return Ok(());
    }
    let db = proj.foil_db();
    let plan = auto::plan(plane, &proj.foils, &db, def, count)?;
    if plan.is_empty() {
        return Ok(());
    }
    eprintln!("  computing {}", plan.describe());
    if let JobOutput::Batch(pols) = run_job(mgr, ui, "2D polar mesh", JobKind::AutoPolars { batches: plan.batches })? {
        for p in pols {
            proj.add_foil_polar(p);
        }
    }
    Ok(())
}

fn real_main() -> Result<(), String> {
    let cli = Cli::parse();
    let Some(cmd) = cli.cmd else {
        if !std::io::stdin().is_terminal() {
            legacy::run(None, cli.verbose);
            return Ok(());
        }
        return Err("no command given. Try `foiley --help`, or pipe an XFoil script on stdin (or use `foiley xfoil`).".into());
    };
    let mut mgr = JobManager::new(cli.threads);
    let mut ui = Ui::new(cli.quiet, cli.verbose, cli.json_progress);
    match cmd {
        Cmd::Xfoil { file } => legacy::run(file, cli.verbose),
        Cmd::Project { file, save } => {
            let ext = file.extension().map(|e| e.to_string_lossy().to_lowercase()).unwrap_or_default();
            let p = if ext == "foiley" || ext == "json" {
                Project::load(&file)?
            } else {
                let mut p = Project::new(&file.file_stem().map(|s| s.to_string_lossy().to_string()).unwrap_or_default());
                let msg = p.import(&file)?;
                eprintln!("  imported {msg}");
                p
            };
            if let Some(out) = &save {
                p.save(out)?;
                eprintln!("  saved {}", out.display());
            }
            println!("project '{}'", p.name);
            println!("  foils ({}): {}", p.foils.len(), p.foils.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(", "));
            println!("  2D polars: {}", p.foil_polars.len());
            for f in &p.foils {
                let n = p.foil_polars.iter().filter(|q| q.foil == f.name).count();
                if n > 0 {
                    println!("    {}: {n} polars", f.name);
                }
            }
            println!("  planes ({}): {}", p.planes.len(), p.planes.iter().map(|f| f.name.as_str()).collect::<Vec<_>>().join(", "));
            for q in &p.plane_polars {
                println!("    {} / {}: {} points", q.plane, q.name, q.points.len());
            }
        }
        Cmd::Foil(fc) => match fc {
            FoilCmd::Info { foils, project } => {
                let proj = load_project(&project)?;
                println!("{:24} {:>6} {:>8} {:>6} {:>8} {:>6} {:>8} {:>8} {:>8}", "foil", "points", "t/c", "x_t", "camber", "x_c", "LE rad", "TE gap", "TE ang");
                for s in &foils {
                    let f = resolve_foil(s, &proj)?;
                    let p = f.props();
                    println!(
                        "{:24} {:>6} {:>7.2}% {:>6.3} {:>7.2}% {:>6.3} {:>8.5} {:>8.5} {:>8.2}",
                        f.name,
                        p.points,
                        p.thickness * 100.0,
                        p.x_thickness,
                        p.camber * 100.0,
                        p.x_camber,
                        p.le_radius,
                        p.te_gap,
                        p.te_angle
                    );
                }
            }
            FoilCmd::Polar { src, x, ptype, alpha, cl, re_sweep, flap_sweep, at_alpha, hinge, out } => {
                let proj = load_project(&src.project)?;
                let foil = resolve_foil(&src.foil, &proj)?;
                let mut st = x.settings()?;
                st.polar_type = match ptype {
                    1 => PolarType::FixedSpeed,
                    2 => PolarType::FixedLift,
                    3 => PolarType::RubberChord,
                    4 => PolarType::FixedAlpha,
                    6 => PolarType::FlapSweep,
                    _ => return Err("--type must be 1, 2, 3, 4 or 6".into()),
                };
                st.alpha_fixed = at_alpha;
                let h = parse_list(&hinge)?;
                if h.len() == 2 {
                    st.flap_hinge = (h[0], h[1]);
                }
                let sweep = match (alpha, cl, re_sweep, flap_sweep) {
                    (_, _, Some(r), _) => {
                        let (a, b, s) = parse_range(&r)?;
                        st.polar_type = PolarType::FixedAlpha;
                        Sweep::Re { from: a, to: b, step: s }
                    }
                    (_, _, _, Some(r)) => {
                        let (a, b, s) = parse_range(&r)?;
                        st.polar_type = PolarType::FlapSweep;
                        Sweep::Flap { from: a, to: b, step: s }
                    }
                    (_, Some(r), _, _) => {
                        let (a, b, s) = parse_range(&r)?;
                        Sweep::Cl { from: a, to: b, step: s }
                    }
                    (Some(r), _, _, _) if r.contains(',') => Sweep::AlphaList(parse_list(&r)?),
                    (Some(r), _, _, _) => {
                        let (a, b, s) = parse_range(&r)?;
                        Sweep::Alpha { from: a, to: b, step: s }
                    }
                    _ => Sweep::Alpha { from: -5.0, to: 15.0, step: 0.5 },
                };
                if let JobOutput::FoilPolar { polar, .. } = run_job(&mut mgr, &mut ui, &format!("{} polar", foil.name), JobKind::FoilPolar { foil, settings: st, sweep })? {
                    let text = match out.format {
                        Format::Table => foil_polar_table(&polar),
                        Format::Csv => polar.to_csv(),
                        Format::Json => serde_json::to_string_pretty(&polar).unwrap(),
                        Format::Xfoil => polar.to_xfoil_text(),
                    };
                    write_out(&out, &text)?;
                }
            }
            FoilCmd::Point { src, x, alpha, cl, out } => {
                let proj = load_project(&src.project)?;
                let foil = resolve_foil(&src.foil, &proj)?;
                let st = x.settings()?;
                let (value, is_alpha) = match (alpha, cl) {
                    (Some(a), _) => (a, true),
                    (None, Some(c)) => (c, false),
                    _ => return Err("give --alpha or --cl".into()),
                };
                if let JobOutput::FoilPoint(d) = run_job(&mut mgr, &mut ui, "point", JobKind::FoilPoint { foil: foil.clone(), settings: st, value, is_alpha })? {
                    let d = d.ok_or("the point did not converge (try another alpha, Re or Ncrit)")?;
                    let text = match out.format {
                        Format::Json => serde_json::to_string_pretty(&d).unwrap(),
                        _ => {
                            let q = &d.point;
                            let mut s = format!("# {} alpha {:.3} CL {:.5} CD {:.6} CDp {:.6} Cm {:.5} Xtr {:.4}/{:.4}\n", foil.name, q.alpha, q.cl, q.cd, q.cdp, q.cm, q.xtr_top, q.xtr_bot);
                            s.push_str("x,y,Cp_inviscid,Cp_viscous\n");
                            for i in 0..d.x.len() {
                                s.push_str(&format!("{:.6},{:.6},{:.5},{:.5}\n", d.x[i], d.y[i], d.cp_inv[i], d.cp_visc[i]));
                            }
                            for (name, b) in [("top", &d.top), ("bottom", &d.bot), ("wake", &d.wake)] {
                                s.push_str(&format!("\n# boundary layer: {name}\nx,s,Ue/Vinf,dstar,theta,H,Cf,N_or_sqrtCtau\n"));
                                for i in 0..b.x.len() {
                                    s.push_str(&format!("{:.6},{:.6},{:.5},{:.4e},{:.4e},{:.4},{:.5e},{:.4}\n", b.x[i], b.s[i], b.ue[i], b.dstar[i], b.theta[i], b.h[i], b.cf[i], b.n_or_ctau[i]));
                                }
                            }
                            s
                        }
                    };
                    write_out(&out, &text)?;
                }
            }
            FoilCmd::Batch { foils, project, re, mach, ncrit, alpha, x, out_dir, save } => {
                let mut proj = match &save {
                    Some(p) if p.exists() => Project::load(p)?,
                    _ => load_project(&project)?,
                };
                let fs: Vec<Foil> = foils.iter().map(|s| resolve_foil(s, &proj)).collect::<Result<_, _>>()?;
                if fs.is_empty() {
                    return Err("give at least one foil".into());
                }
                let res = if re.contains(':') {
                    let (a, b, n) = parse_range(&re)?;
                    foiley_foil::batch::re_list_geometric(a, b, n as usize)
                } else {
                    parse_list(&re)?
                };
                let mut cases = Vec::new();
                for &r in &res {
                    for &m in &parse_list(&mach)? {
                        for &n in &parse_list(&ncrit)? {
                            cases.push(BatchCase { re: r, mach: m, ncrit: n });
                        }
                    }
                }
                let (a, b, s) = parse_range(&alpha)?;
                let st = x.settings()?;
                if let JobOutput::Batch(pols) = run_job(&mut mgr, &mut ui, "batch", JobKind::Batch { foils: fs.clone(), cases, settings: st, sweep: Sweep::Alpha { from: a, to: b, step: s } })? {
                    if let Some(d) = &out_dir {
                        std::fs::create_dir_all(d).map_err(|e| e.to_string())?;
                    }
                    for p in &pols {
                        if let Some(d) = &out_dir {
                            let fname = format!("{}_{}.csv", p.foil.replace([' ', '/', '\\'], "_"), p.name);
                            std::fs::write(d.join(&fname), p.to_csv()).map_err(|e| e.to_string())?;
                        } else {
                            print!("{}", foil_polar_table(p));
                        }
                    }
                    if let Some(path) = &save {
                        for f in fs {
                            if proj.foil(&f.name).is_none() {
                                proj.foils.push(f);
                            }
                        }
                        for p in pols {
                            proj.add_foil_polar(p);
                        }
                        proj.save(path)?;
                        eprintln!("  saved {}", path.display());
                    }
                }
            }
            FoilCmd::Optimize { src, x, alpha, objectives, algorithm, population, iterations, bumps, amplitude, flap, seed, out, pareto } => {
                use foiley_foil::optim::{Algorithm, Goal, Objective, ObjectiveSpec, OptimSettings};
                let proj = load_project(&src.project)?;
                let foil = resolve_foil(&src.foil, &proj)?;
                let mut st = OptimSettings::default();
                st.algorithm = match algorithm.to_lowercase().as_str() {
                    "ga" | "genetic" => Algorithm::Genetic,
                    _ => Algorithm::Mopso,
                };
                st.analysis = x.settings()?;
                st.alpha = alpha;
                st.population = population;
                st.iterations = iterations;
                st.bumps_per_side = bumps;
                st.max_amplitude = amplitude;
                st.seed = seed;
                if let Some(f) = flap {
                    let v = parse_list(&f)?;
                    if v.len() != 2 {
                        return Err("--flap needs min,max".into());
                    }
                    st.flap = Some((v[0], v[1]));
                }
                st.objectives = objectives
                    .split(',')
                    .filter(|t| !t.trim().is_empty())
                    .map(|t| {
                        let parts: Vec<&str> = t.trim().split(':').collect();
                        let obj = match parts[0].to_lowercase().as_str() {
                            "cl" => Objective::Cl,
                            "cd" => Objective::Cd,
                            "cl/cd" | "ld" => Objective::ClCd,
                            "cpmin" => Objective::Cpmin,
                            "cm" => Objective::Cm,
                            o => return Err(format!("unknown objective '{o}'")),
                        };
                        let goal = match parts.get(1).map(|s| s.to_lowercase()) {
                            None => Goal::Maximise,
                            Some(g) if g == "max" => Goal::Maximise,
                            Some(g) if g == "min" => Goal::Minimise,
                            Some(g) => Goal::Target { value: g.parse().map_err(|_| format!("bad target '{g}'"))?, tolerance: parts.get(2).and_then(|t| t.parse().ok()).unwrap_or(1e-3) },
                        };
                        Ok(ObjectiveSpec { objective: obj, enabled: true, goal })
                    })
                    .collect::<Result<_, String>>()?;
                if let JobOutput::Optimise(r) = run_job(&mut mgr, &mut ui, &format!("optimise {}", foil.name), JobKind::Optimise { foil: foil.clone(), settings: st })? {
                    let best = r.best.as_ref().ok_or("no feasible design found")?;
                    println!("# {} evaluations, {} Pareto designs", r.evaluations, r.pareto.len());
                    println!("# base:  see `foiley foil point`; best: CL {:.4} CD {:.5} CL/CD {:.2} Cm {:.4} Cpmin {:.3}", best.cl, best.cd, best.cl / best.cd.max(1e-12), best.cm, best.cpmin);
                    if let Some(p) = pareto {
                        let mut s = String::from("index,CL,CD,CL/CD,Cm,Cpmin\n");
                        for (i, c) in r.pareto.iter().enumerate() {
                            s.push_str(&format!("{i},{:.5},{:.6},{:.3},{:.5},{:.4}\n", c.cl, c.cd, c.cl / c.cd.max(1e-12), c.cm, c.cpmin));
                        }
                        std::fs::write(p, s).map_err(|e| e.to_string())?;
                    }
                    let mut f = best.foil.clone();
                    f.name = format!("{} opt", foil.name);
                    match out {
                        Some(o) => f.save(&o).map_err(|e| e.to_string())?,
                        None => print!("{}", f.to_dat()),
                    }
                }
            }
            FoilCmd::Inverse { src, alpha, mode, segment, scale, bump, smooth, iterations, speed_csv, out } => {
                use foiley_foil::inverse::{modify, InverseMode, ModKind, Modification};
                let proj = load_project(&src.project)?;
                let foil = resolve_foil(&src.foil, &proj)?;
                let mode = if mode.to_lowercase().starts_with('f') { InverseMode::Full } else { InverseMode::Mixed };
                let parts: Vec<&str> = segment.split(',').collect();
                if parts.len() != 3 {
                    return Err("--segment must be upper|lower,x_from,x_to".into());
                }
                let upper = parts[0].trim().to_lowercase().starts_with('u');
                let (x0, x1): (f64, f64) = (parts[1].trim().parse().map_err(|_| "bad x_from")?, parts[2].trim().parse().map_err(|_| "bad x_to")?);
                let kind = match (scale, bump, smooth) {
                    (Some(f), _, _) => ModKind::Scale(f),
                    (_, Some(a), _) => ModKind::Bump(a),
                    (_, _, true) => ModKind::Smooth,
                    _ => return Err("give --scale, --bump or --smooth".into()),
                };
                let pan = foiley_foil::foil::Paneling::default();
                let d0 = match run_job(&mut mgr, &mut ui, "speed distribution", JobKind::Inverse { foil: foil.clone(), alpha, mode, paneling: pan, target: None, iterations })? {
                    JobOutput::Inverse(_, d) => d,
                    _ => unreachable!(),
                };
                let (t, seg) = modify(&d0, &Modification { x_from: x0, x_to: x1, upper, kind });
                if seg.is_none() {
                    return Err("the segment contains fewer than 2 points".into());
                }
                if let JobOutput::Inverse(Some(f), d1) = run_job(&mut mgr, &mut ui, "inverse design", JobKind::Inverse { foil: foil.clone(), alpha, mode, paneling: pan, target: Some((t.clone(), seg)), iterations })? {
                    let p = f.props();
                    eprintln!("  new foil: t/c {:.2}% at {:.3}, camber {:.2}% at {:.3}, CL {:.4} at {alpha} deg", p.thickness * 100.0, p.x_thickness, p.camber * 100.0, p.x_camber, d1.cl);
                    if let Some(path) = speed_csv {
                        let mut s = String::from("x,s,q_initial,q_target\n");
                        for i in 0..d0.q.len() {
                            s.push_str(&format!("{:.6},{:.6},{:.6},{:.6}\n", d0.x[i], d0.s[i], d0.q[i], t[i]));
                        }
                        std::fs::write(path, s).map_err(|e| e.to_string())?;
                    }
                    match out {
                        Some(o) => f.save(&o).map_err(|e| e.to_string())?,
                        None => print!("{}", f.to_dat()),
                    }
                }
            }
            FoilCmd::Design { src, flap, te_gap, le_radius, thickness, camber, high_points, repanel, normalize, mix, name, out } => {
                let proj = load_project(&src.project)?;
                let mut f = resolve_foil(&src.foil, &proj)?;
                let e = |r: Result<Foil, foiley_foil::FoilError>| r.map_err(|e| e.to_string());
                if normalize {
                    f = f.normalized();
                }
                if let Some(m) = mix {
                    let (other, t) = m.rsplit_once(',').ok_or("--mix needs FOIL,fraction")?;
                    let g = resolve_foil(other, &proj)?;
                    let t: f64 = t.trim().parse().map_err(|_| "bad mix fraction")?;
                    f = e(Foil::interpolate(&f, &g, t, &f.name.clone()))?;
                }
                if thickness.is_some() || camber.is_some() {
                    f = e(f.with_tc(thickness, camber))?;
                }
                if let Some(h) = high_points {
                    let v = parse_list(&h)?;
                    if v.len() != 2 {
                        return Err("--high-points needs xt,xc".into());
                    }
                    f = e(f.with_high_points(v[0], v[1]))?;
                }
                if let Some(g) = le_radius {
                    let v = parse_list(&g)?;
                    f = e(f.with_le_radius(v[0], *v.get(1).unwrap_or(&0.1)))?;
                }
                if let Some(g) = te_gap {
                    let v = parse_list(&g)?;
                    f = e(f.with_te_gap(v[0], *v.get(1).unwrap_or(&0.8)))?;
                }
                if let Some(fl) = flap {
                    let v = parse_list(&fl)?;
                    if v.len() != 3 {
                        return Err("--flap needs x/c,y/t,deg".into());
                    }
                    f = e(f.with_flap(v[0], v[1], v[2]))?;
                }
                if let Some(n) = repanel {
                    f = e(f.repaneled(&Paneling { panels: n, ..Default::default() }))?;
                }
                if let Some(n) = name {
                    f.name = n;
                }
                let p = f.props();
                eprintln!("  {}: {} points, t/c {:.2}% at {:.3}, camber {:.2}% at {:.3}", f.name, p.points, p.thickness * 100.0, p.x_thickness, p.camber * 100.0, p.x_camber);
                match out {
                    Some(o) => f.save(&o).map_err(|e| e.to_string())?,
                    None => print!("{}", f.to_dat()),
                }
            }
        },
        Cmd::Plane(pc) => match pc {
            PlaneCmd::Info { src } => {
                let (proj, plane) = load_plane(&src)?;
                let db = proj.foil_db();
                println!("plane '{}'", plane.name);
                for w in &plane.wings {
                    println!(
                        "  {:14} {:10} span {:.4} m  area {:.4} m2  MAC {:.4} m  AR {:.2}  taper {:.3}  sweep {:.2} deg  sections {}",
                        w.name,
                        w.kind.label(),
                        w.planform_span(),
                        w.planform_area(),
                        w.mac(),
                        w.aspect_ratio(),
                        w.taper(),
                        w.sweep(),
                        w.sections.len()
                    );
                }
                for f in &plane.fuses {
                    println!("  fuselage {:10} length {:.3} m  wetted area {:.4} m2", f.name, f.length(), f.wetted_area());
                }
                let i = plane.inertia(&db);
                println!("  mass {:.4} kg  CoG ({:.4}, {:.4}, {:.4}) m", i.mass, i.cog.x, i.cog.y, i.cog.z);
                println!("  inertia about CoG: Ixx {:.5}  Iyy {:.5}  Izz {:.5}  Ixz {:.5} kg.m2", i.ixx, i.iyy, i.izz, i.ixz);
                let (h, v) = plane.tail_volumes();
                println!("  tail volumes: horizontal {h:.3}, vertical {v:.4}");
                println!("  foils: {}", plane.foil_names().join(", "));
            }
            PlaneCmd::Export { src, avl, mass, stl, thin, xml } => {
                let (proj, plane) = load_plane(&src)?;
                let db = proj.foil_db();
                if let Some(p) = avl {
                    std::fs::write(&p, foiley_plane::export::avl(&plane, &db, 0.0)).map_err(|e| e.to_string())?;
                    let dir = p.parent().unwrap_or(Path::new("."));
                    for (name, text) in foiley_plane::export::avl_foil_files(&plane, &proj.foils) {
                        std::fs::write(dir.join(name), text).map_err(|e| e.to_string())?;
                    }
                    eprintln!("  wrote {}", p.display());
                }
                if let Some(p) = mass {
                    std::fs::write(&p, foiley_plane::export::avl_mass(&plane, &db)).map_err(|e| e.to_string())?;
                    eprintln!("  wrote {}", p.display());
                }
                if let Some(p) = stl {
                    std::fs::write(&p, foiley_plane::export::stl_mesh(&plane, &db, !thin).map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
                    eprintln!("  wrote {}", p.display());
                }
                if let Some(p) = xml {
                    std::fs::write(&p, foiley_plane::xml::write_plane(&plane)).map_err(|e| e.to_string())?;
                    eprintln!("  wrote {}", p.display());
                }
            }
            PlaneCmd::ImportFuselage { src, file, units, triangles, mass, save } => {
                let (mut proj, plane) = load_plane(&src)?;
                let idx = match proj.planes.iter().position(|p| p.name == plane.name) {
                    Some(i) => i,
                    None => {
                        proj.planes.push(plane.clone());
                        proj.planes.len() - 1
                    }
                };
                let scale = match units.to_lowercase().as_str() {
                    "auto" => None,
                    "m" => Some(1.0),
                    "mm" => Some(1e-3),
                    "cm" => Some(1e-2),
                    "in" | "inch" => Some(0.0254),
                    "ft" | "foot" => Some(0.3048),
                    u => return Err(format!("unknown unit '{u}'")),
                };
                let msg = proj.import_fuselage(&file, idx, scale, triangles)?;
                if let Some(f) = proj.planes[idx].fuses.first_mut() {
                    f.mass = mass;
                }
                eprintln!("  {msg}");
                proj.save(&save)?;
                eprintln!("  saved {}", save.display());
            }
            PlaneCmd::Analyze { src, ptype, method, free_wake, wake_nseg, wake_relax, wake_tol, response, response_amp, response_time, response_out, range, speed, alpha, beta, list, inviscid, altitude, water, mass, cog, ground, derivatives, strict_envelope, flaps, gains, reference, auto_polars, llt_stations, strips, save, out } => {
                let (mut proj, plane) = load_plane(&src)?;
                let kind = match ptype.to_lowercase().trim_start_matches('t') {
                    "1" => PolarKind::FixedSpeed,
                    "2" => PolarKind::FixedLift,
                    "3" => PolarKind::Glide,
                    "4" => PolarKind::FixedAlpha,
                    "5" => PolarKind::Sideslip,
                    "6" => return Err("T6 control polars are defined in JSON job scripts (`foiley run`), since they need several ranges".into()),
                    "7" => PolarKind::Stability,
                    "8" => PolarKind::List,
                    _ => return Err("--type must be t1..t8".into()),
                };
                let method = match method.to_lowercase().as_str() {
                    "llt" => Method::Llt,
                    "vlm1" => Method::Vlm1,
                    "vlm2" | "vlm" => Method::Vlm2,
                    "panel" | "panels" | "3dpanel" => Method::Panel,
                    "tripanel" | "triangles" | "tri" => Method::TriPanel,
                    "trilinear" => Method::TriLinear,
                    m => return Err(format!("unknown method '{m}' (llt, vlm1, vlm2, panel, tripanel, trilinear)")),
                };
                let (a, b, s) = parse_range(&range)?;
                let fluid = if water { Fluid::water() } else { let air = isa(altitude, 0.0); Fluid { rho: air.rho, mu: air.mu, a: air.a } };
                let mut def = PlanePolarDef { kind, method, velocity: speed, alpha, beta, range: Range { from: a, to: b, step: s }, viscous: !inviscid, fluid, ground_height: ground, derivatives, ..Default::default() };
                def.llt.stations = llt_stations;
                def.free_wake = free_wake;
                def.wake_nseg = wake_nseg;
                def.wake_relax = wake_relax;
                def.wake_tol = wake_tol;
                def.strict_envelope = strict_envelope;
                if kind == PolarKind::Stability && a == -4.0 && b == 10.0 {
                    def.range = Range { from: 0.0, to: 0.0, step: 1.0 };
                }
                if let Some(l) = list {
                    def.list = l
                        .split(';')
                        .filter(|t| !t.trim().is_empty())
                        .map(|t| {
                            let v = parse_list(t)?;
                            if v.len() != 3 {
                                return Err(format!("t8 point '{t}' needs alpha,beta,V"));
                            }
                            Ok((v[0], v[1], v[2]))
                        })
                        .collect::<Result<_, String>>()?;
                }
                def.reference = match reference.to_lowercase().as_str() {
                    "projected" => RefDims::Projected,
                    "planform" => RefDims::Planform,
                    r => {
                        let v = parse_list(r)?;
                        if v.len() != 3 {
                            return Err("--reference must be projected, planform or area,span,chord".into());
                        }
                        RefDims::Custom { area: v[0], span: v[1], chord: v[2] }
                    }
                };
                if mass.is_some() || cog.is_some() {
                    let mut i: Inertia = plane.inertia(&proj.foil_db());
                    if let Some(m) = mass {
                        let k = if i.mass > 0.0 { m / i.mass } else { 1.0 };
                        i.ixx *= k;
                        i.iyy *= k;
                        i.izz *= k;
                        i.ixz *= k;
                        i.mass = m;
                    }
                    if let Some(c) = cog {
                        let v = parse_list(&c)?;
                        if v.len() != 3 {
                            return Err("--cog needs x,y,z".into());
                        }
                        i.cog = Vec3::new(v[0], v[1], v[2]);
                    }
                    def.inertia = Some(i);
                }
                use foiley_plane::analysis::ControlTarget;
                for f in &flaps {
                    let v = parse_list(f)?;
                    if v.len() != 3 {
                        return Err("--flap needs wing,flap,deg".into());
                    }
                    def.flaps.push(((v[0] as usize, v[1] as usize), v[2]));
                }
                for g in &gains {
                    let v = parse_list(g)?;
                    if v.len() != 3 {
                        return Err("--gain needs wing,flap,deg_per_unit".into());
                    }
                    def.gains.push((ControlTarget::Flap { wing: v[0] as usize, index: v[1] as usize }, v[2]));
                }
                def.name = def.auto_name();
                // the shared Method::supports table, before any 2D polar work
                def.validate()?;
                ensure_polars(&mut mgr, &mut ui, &mut proj, &plane, &def, auto_polars)?;
                let db = proj.foil_db();
                let title = format!("{} {}", plane.name, def.name);
                if let JobOutput::PlanePolar(pol) = run_job(&mut mgr, &mut ui, &title, JobKind::PlanePolar { plane: plane.clone(), db, def })? {
                    let text = match out.format {
                        Format::Table | Format::Xfoil => plane_polar_table(&pol),
                        Format::Csv => pol.to_csv(),
                        Format::Json => serde_json::to_string_pretty(&pol).unwrap(),
                    };
                    write_out(&out, &text)?;
                    if let Some(spec) = &response {
                        use foiley_plane::stability::{forced_response, ControlInput};
                        let (kind, arg) = spec.split_once(':').unwrap_or((spec.as_str(), "0.5"));
                        let a: f64 = arg.parse().map_err(|_| "bad --response argument")?;
                        let input = match kind {
                            "step" => ControlInput::Step,
                            "pulse" => ControlInput::Pulse { duration: a },
                            "doublet" => ControlInput::Doublet { duration: a },
                            "sine" => ControlInput::Sine { hz: a },
                            _ => return Err("--response must be step, pulse:T, doublet:T or sine:HZ".into()),
                        };
                        let s = pol.points.iter().find_map(|o| o.stability.as_ref()).ok_or("no stability results (use --type t7 with --gain)")?;
                        if s.control.is_none() {
                            return Err("no control derivatives: give --gain wing,flap,deg_per_unit".into());
                        }
                        let dt = (response_time / 2000.0).max(1e-4);
                        let rl = forced_response(&s.a_long, &s.b_long, input, response_amp, response_time, dt);
                        let rt = forced_response(&s.a_lat, &s.b_lat, input, response_amp, response_time, dt);
                        let mut csv = String::from("t,u_mps,w_mps,q_degps,theta_deg,v_mps,p_degps,r_degps,phi_deg,input\n");
                        let r2d = 180.0 / std::f64::consts::PI;
                        for ((t, a), (_, b)) in rl.iter().zip(&rt) {
                            csv.push_str(&format!("{t:.4},{:.5},{:.5},{:.5},{:.5},{:.5},{:.5},{:.5},{:.5},{:.4}\n", a[0], a[1], a[2] * r2d, a[3] * r2d, b[0], b[1] * r2d, b[2] * r2d, b[3] * r2d, response_amp * input.value(*t)));
                        }
                        match &response_out {
                            Some(p) => std::fs::write(p, csv).map_err(|e| e.to_string())?,
                            None => print!("{csv}"),
                        }
                    }
                    if let Some(sp) = strips {
                        std::fs::write(&sp, strips_csv(&pol)).map_err(|e| e.to_string())?;
                    }
                    if let Some(path) = save {
                        if proj.planes.iter().all(|p| p.name != plane.name) {
                            proj.planes.push(plane);
                        }
                        proj.add_plane_polar(pol);
                        proj.save(&path)?;
                        eprintln!("  saved {}", path.display());
                    }
                }
            }
        },
        Cmd::Run { script } => run_script(&mut mgr, &mut ui, &script)?,
    }
    Ok(())
}

/// JSON job script.
#[derive(serde::Deserialize)]
struct Script {
    #[serde(default)]
    project: Option<PathBuf>,
    #[serde(default)]
    save: Option<PathBuf>,
    #[serde(default)]
    output_dir: Option<PathBuf>,
    steps: Vec<Step>,
}

#[derive(serde::Deserialize)]
#[serde(rename_all = "snake_case")]
enum Step {
    /// import a file (.dat, .xml, polar .txt, project)
    Import(PathBuf),
    /// add a NACA foil
    Naca(String),
    FoilPolar {
        foil: String,
        #[serde(default)]
        settings: AnalysisSettings,
        sweep: Sweep,
    },
    Batch {
        foils: Vec<String>,
        re: Vec<f64>,
        #[serde(default)]
        mach: Vec<f64>,
        #[serde(default)]
        ncrit: Vec<f64>,
        #[serde(default)]
        settings: AnalysisSettings,
        sweep: Sweep,
    },
    PlanePolar {
        #[serde(default)]
        plane: Option<String>,
        def: PlanePolarDef,
        #[serde(default = "six")]
        auto_polars: usize,
    },
    /// optimise a foil; the best design is added to the project as "<foil> opt"
    Optimize {
        foil: String,
        #[serde(default)]
        settings: foiley_foil::optim::OptimSettings,
    },
    /// replace a plane's fuselage with an STL / STEP surface
    ImportFuselage {
        file: PathBuf,
        #[serde(default)]
        plane: Option<String>,
        /// file-unit to metre factor (default: automatic)
        #[serde(default)]
        scale: Option<f64>,
    },
}

fn six() -> usize {
    6
}

fn run_script(mgr: &mut JobManager, ui: &mut Ui, path: &Path) -> Result<(), String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let sc: Script = serde_json::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
    let base = path.parent().unwrap_or(Path::new("."));
    let rel = |p: &Path| if p.is_absolute() { p.to_path_buf() } else { base.join(p) };
    let mut proj = match &sc.project {
        Some(p) => Project::load(&rel(p))?,
        None => Project::new("script"),
    };
    let outdir = sc.output_dir.as_ref().map(|d| rel(d));
    if let Some(d) = &outdir {
        std::fs::create_dir_all(d).map_err(|e| e.to_string())?;
    }
    let n = sc.steps.len();
    for (k, step) in sc.steps.into_iter().enumerate() {
        eprintln!("step {}/{n}", k + 1);
        match step {
            Step::Import(p) => eprintln!("  imported {}", proj.import(&rel(&p))?),
            Step::Naca(c) => {
                let f = Foil::naca(&c).map_err(|e| e.to_string())?;
                eprintln!("  added {}", f.name);
                proj.upsert_foil(f);
            }
            Step::FoilPolar { foil, settings, sweep } => {
                let f = resolve_foil(&foil, &proj)?;
                if let JobOutput::FoilPolar { polar, .. } = run_job(mgr, ui, &format!("{} polar", f.name), JobKind::FoilPolar { foil: f, settings, sweep })? {
                    if let Some(d) = &outdir {
                        std::fs::write(d.join(format!("{}_{}.csv", polar.foil.replace(' ', "_"), polar.name)), polar.to_csv()).map_err(|e| e.to_string())?;
                    }
                    proj.add_foil_polar(polar);
                }
            }
            Step::Batch { foils, re, mach, ncrit, settings, sweep } => {
                let fs: Vec<Foil> = foils.iter().map(|s| resolve_foil(s, &proj)).collect::<Result<_, _>>()?;
                let mach = if mach.is_empty() { vec![0.0] } else { mach };
                let ncrit = if ncrit.is_empty() { vec![9.0] } else { ncrit };
                let mut cases = Vec::new();
                for &r in &re {
                    for &m in &mach {
                        for &nc in &ncrit {
                            cases.push(BatchCase { re: r, mach: m, ncrit: nc });
                        }
                    }
                }
                if let JobOutput::Batch(pols) = run_job(mgr, ui, "batch", JobKind::Batch { foils: fs, cases, settings, sweep })? {
                    for p in pols {
                        if let Some(d) = &outdir {
                            std::fs::write(d.join(format!("{}_{}.csv", p.foil.replace(' ', "_"), p.name)), p.to_csv()).map_err(|e| e.to_string())?;
                        }
                        proj.add_foil_polar(p);
                    }
                }
            }
            Step::Optimize { foil, settings } => {
                let f = resolve_foil(&foil, &proj)?;
                if let JobOutput::Optimise(r) = run_job(mgr, ui, &format!("optimise {}", f.name), JobKind::Optimise { foil: f.clone(), settings })? {
                    if let Some(b) = r.best {
                        let mut g = b.foil.clone();
                        g.name = proj.unique_foil_name(&format!("{} opt", f.name));
                        eprintln!("  best: CL {:.4} CD {:.5} -> foil '{}'", b.cl, b.cd, g.name);
                        proj.foils.push(g);
                    }
                }
            }
            Step::ImportFuselage { file, plane, scale } => {
                let idx = match &plane {
                    Some(n) => proj.planes.iter().position(|p| &p.name == n).ok_or_else(|| format!("no plane '{n}'"))?,
                    None => 0,
                };
                eprintln!("  {}", proj.import_fuselage(&rel(&file), idx, scale, 3000)?);
            }
            Step::PlanePolar { plane, def, auto_polars } => {
                let pl = match &plane {
                    Some(n) => proj.planes.iter().find(|p| &p.name == n).cloned().ok_or_else(|| format!("no plane '{n}'"))?,
                    None => proj.planes.first().cloned().ok_or("no plane loaded")?,
                };
                for name in pl.foil_names() {
                    if proj.foil(&name).is_none() {
                        let l = name.to_lowercase().replace([' ', '-', '_'], "");
                        if let Some(code) = l.strip_prefix("naca") {
                            let mut f = Foil::naca(code).map_err(|e| e.to_string())?;
                            f.name = name.clone();
                            proj.foils.push(f);
                        }
                    }
                }
                let mut def = def;
                if def.name.is_empty() {
                    def.name = def.auto_name();
                }
                ensure_polars(mgr, ui, &mut proj, &pl, &def, auto_polars)?;
                let db = proj.foil_db();
                if let JobOutput::PlanePolar(pol) = run_job(mgr, ui, &format!("{} {}", pl.name, def.name), JobKind::PlanePolar { plane: pl, db, def })? {
                    print!("{}", plane_polar_table(&pol));
                    if let Some(d) = &outdir {
                        std::fs::write(d.join(format!("{}_{}.csv", pol.plane.replace(' ', "_"), pol.name)), pol.to_csv()).map_err(|e| e.to_string())?;
                        if pol.points.iter().any(|q| !q.strips.is_empty()) {
                            std::fs::write(d.join(format!("{}_{}_strips.csv", pol.plane.replace(' ', "_"), pol.name)), strips_csv(&pol)).map_err(|e| e.to_string())?;
                        }
                    }
                    proj.add_plane_polar(pol);
                }
            }
        }
    }
    if let Some(s) = &sc.save {
        proj.save(&rel(s))?;
        eprintln!("saved {}", rel(s).display());
    }
    Ok(())
}

fn main() {
    match real_main() {
        Ok(()) => {}
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(if e == "cancelled" { 130 } else { 1 });
        }
    }
}
