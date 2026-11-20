//! Windows (MSVC): embed `foiley.manifest` in `foiley.exe`, so that started
//! from Explorer it gets no console window (see `src/launch.rs`), and give its
//! main thread an 8 MB stack.

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed=foiley.manifest");
    let var = |k: &str| std::env::var(k).unwrap_or_default();
    if var("CARGO_CFG_TARGET_OS") == "windows" && var("CARGO_CFG_TARGET_ENV") == "msvc" {
        let manifest = std::path::Path::new(&var("CARGO_MANIFEST_DIR")).join("foiley.manifest");
        println!("cargo:rustc-link-arg-bin=foiley=/MANIFEST:EMBED");
        println!("cargo:rustc-link-arg-bin=foiley=/MANIFESTINPUT:{}", manifest.display());
        // `real_main` dispatches every subcommand in one frame; unoptimised
        // builds overflow MSVC's default 1 MB main-thread stack. 8 MB matches Linux.
        println!("cargo:rustc-link-arg-bin=foiley=/STACK:8388608");
    }
}
