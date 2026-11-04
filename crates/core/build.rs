//! Embeds the git commit, build profile and compile-time target features
//! for run manifests (`foiley_core::provenance`).

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    let out = Command::new("git").args(args).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let s = String::from_utf8(out.stdout).ok()?.trim().to_string();
    (!s.is_empty()).then_some(s)
}

fn main() {
    let commit = std::env::var("FOILEY_GIT_COMMIT").ok().filter(|s| !s.is_empty()).or_else(|| git(&["rev-parse", "HEAD"])).unwrap_or_else(|| "unknown".into());
    println!("cargo:rustc-env=FOILEY_GIT_COMMIT={commit}");
    println!("cargo:rerun-if-env-changed=FOILEY_GIT_COMMIT");
    // rebuild when HEAD moves (worktrees keep HEAD in their own git dir and
    // the refs in the common one)
    for (flag, files) in [("--git-dir", &["HEAD"][..]), ("--git-common-dir", &["refs/heads", "packed-refs"][..])] {
        if let Some(dir) = git(&["rev-parse", flag]) {
            for f in files {
                let p = std::path::Path::new(&dir).join(f);
                if p.exists() {
                    println!("cargo:rerun-if-changed={}", p.display());
                }
            }
        }
    }
    let profile = std::env::var("PROFILE").unwrap_or_else(|_| "unknown".into());
    println!("cargo:rustc-env=FOILEY_BUILD_PROFILE={profile}");
    let features = std::env::var("CARGO_CFG_TARGET_FEATURE").unwrap_or_default();
    println!("cargo:rustc-env=FOILEY_TARGET_FEATURES={features}");
    let target = std::env::var("TARGET").unwrap_or_default();
    println!("cargo:rustc-env=FOILEY_TARGET={target}");
    println!("cargo:rerun-if-changed=build.rs");
}
