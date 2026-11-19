use super::StudioModule;
use clap::{Args, ValueEnum};
use foiley_gui::snap::{self, Crop, Options, Preset};
use std::path::PathBuf;

#[derive(Args)]
#[command(
    after_help = "Steps: settle [seconds], wait FRAMES, module NAME, snap NAME [full|viewport], key [MOD+]KEY, \
text \"TEXT\", click POS [BUTTON] or click \"Label\", hover POS, drag [MOD+]BUTTON POS DX DY [FRAMES], \
scroll DY [POS], pinch FACTOR [POS], trace NAME, end.\n\
Positions are window pixels x,y or viewport fractions vp:fx,fy. In a file, # starts a comment.\n\
Exit status: 0 success, 1 argument or script error, 2 settle timeout or job failure, 3 no render adapter."
)]
pub struct SnapArgs {
    /// Project, geometry or polar files, as for Studio
    files: Vec<PathBuf>,
    /// Module to start in
    #[arg(long, value_enum)]
    module: Option<StudioModule>,
    /// Run a panel analysis of the first plane
    #[arg(long)]
    run: bool,
    /// Image size in window pixels
    #[arg(long, default_value = "1400x900", value_parser = size_parser)]
    size: [u32; 2],
    /// Pixels per egui point
    #[arg(long, default_value_t = 1.0)]
    scale: f32,
    /// Studio theme
    #[arg(long, default_value = "dark", value_parser = ["dark", "light"])]
    theme: String,
    /// Viewport navigation preset
    #[arg(long, value_parser = ["cad", "standard"])]
    nav: Option<String>,
    /// Use perspective. Use --persp=false for orthographic
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "true")]
    persp: Option<bool>,
    /// Show the ground grid. Use --grid=false to hide it
    #[arg(long, num_args = 0..=1, require_equals = true, default_missing_value = "true")]
    grid: Option<bool>,
    /// Output directory (defaults to target/snap/<timestamp>)
    #[arg(long)]
    out: Option<PathBuf>,
    /// Script file with one step per line
    #[arg(long, conflicts_with = "steps")]
    script: Option<PathBuf>,
    /// Step or semicolon-separated steps. Repeat to append steps
    #[arg(long = "do")]
    steps: Vec<String>,
    /// Built-in viewport check
    #[arg(long, value_enum)]
    preset: Option<SnapPreset>,
    /// Also write an animated GIF for each trace
    #[arg(long)]
    gif: bool,
    /// Default crop for captures
    #[arg(long, default_value = "full", value_parser = Crop::parse)]
    crop: Crop,
}

#[derive(Clone, Copy, ValueEnum)]
enum SnapPreset {
    Views,
    Orbit,
    Nav,
    Analysis,
    Narrow,
}

impl From<SnapPreset> for Preset {
    fn from(preset: SnapPreset) -> Self {
        match preset {
            SnapPreset::Views => Self::Views,
            SnapPreset::Orbit => Self::Orbit,
            SnapPreset::Nav => Self::Nav,
            SnapPreset::Analysis => Self::Analysis,
            SnapPreset::Narrow => Self::Narrow,
        }
    }
}

fn size_parser(value: &str) -> Result<[u32; 2], String> {
    let error = || format!("expected a size in WxH pixels, got {value:?}");
    let (width, height) = value.split_once(['x', 'X']).ok_or_else(error)?;
    Ok([
        width.parse().map_err(|_| error())?,
        height.parse().map_err(|_| error())?,
    ])
}

impl SnapArgs {
    pub fn run(self, threads: usize) {
        let launch = foiley_gui::Launch {
            files: self.files,
            module: self.module.map(Into::into),
            run: self.run,
            threads,
        };
        match snap::run(Options {
            launch,
            size: self.size,
            scale: self.scale,
            light: self.theme == "light",
            nav: self.nav,
            perspective: self.persp,
            grid: self.grid,
            out: self.out,
            script: self.script,
            steps: self.steps,
            preset: self.preset.map(Into::into),
            gif: self.gif,
            crop: self.crop,
        }) {
            Ok(files) => {
                for file in files {
                    println!("{}", file.display());
                }
            }
            Err(error) => {
                eprintln!("foiley snap: {error}");
                std::process::exit(error.code());
            }
        }
    }
}
