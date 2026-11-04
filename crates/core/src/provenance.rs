//! Run manifests: what produced a 2D or 3D polar (build, threads, settings,
//! geometry, mesh and polar families), hashed once per run.

use crate::hash::{Hash32, Hasher};
use serde::{Deserialize, Serialize};

/// What this binary is.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct BuildInfo {
    pub foiley_version: String,
    /// `git rev-parse HEAD` at build time ("unknown" outside a checkout)
    pub git_commit: String,
    /// cargo profile ("release", "debug")
    pub build_profile: String,
    /// target triple
    pub target: String,
    /// compile-time target features (`CARGO_CFG_TARGET_FEATURE`)
    pub target_features: Vec<String>,
}

impl BuildInfo {
    pub fn current() -> BuildInfo {
        BuildInfo {
            foiley_version: env!("CARGO_PKG_VERSION").into(),
            git_commit: env!("FOILEY_GIT_COMMIT").into(),
            build_profile: env!("FOILEY_BUILD_PROFILE").into(),
            target: env!("FOILEY_TARGET").into(),
            target_features: env!("FOILEY_TARGET_FEATURES").split(',').filter(|s| !s.is_empty()).map(String::from).collect(),
        }
    }
}

/// Threads available to the run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct ThreadConfig {
    /// rayon worker threads
    pub rayon_threads: usize,
    /// logical CPUs of the machine
    pub available: usize,
}

impl ThreadConfig {
    pub fn current() -> ThreadConfig {
        ThreadConfig { rayon_threads: rayon::current_num_threads(), available: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) }
    }
}

/// A 2D polar family a 3D run used.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct FamilyRef {
    pub foil: String,
    /// index in the run's foil database (only meaningful within the run)
    pub family_id: u32,
    /// hash of the family key (`FamilyKey::hash`)
    pub key_hash: Hash32,
    pub foil_geom_hash: Hash32,
    pub flap_deg: f64,
}

/// Provenance of one polar.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct RunManifest {
    pub foiley_version: String,
    pub git_commit: String,
    pub build_profile: String,
    pub target: String,
    pub target_features: Vec<String>,
    pub thread_config: ThreadConfig,
    /// hash of the settings snapshot stored on the polar
    /// (`FoilPolar::settings`, `PlanePolar::def`)
    pub settings: Hash32,
    /// foil (2D) or plane (3D) geometry
    pub geometry_hash: Hash32,
    /// 3D panel methods: the mesh of the base configuration
    #[serde(default)]
    pub mesh_hash: Option<Hash32>,
    /// 3D viscous runs: the 2D polar families used
    #[serde(default)]
    pub families: Vec<FamilyRef>,
    /// UTC, RFC 3339
    pub started_at: String,
}

impl RunManifest {
    /// A manifest for a run starting now.
    pub fn new(settings: Hash32, geometry_hash: Hash32) -> RunManifest {
        let b = BuildInfo::current();
        RunManifest {
            foiley_version: b.foiley_version,
            git_commit: b.git_commit,
            build_profile: b.build_profile,
            target: b.target,
            target_features: b.target_features,
            thread_config: ThreadConfig::current(),
            settings,
            geometry_hash,
            mesh_hash: None,
            families: vec![],
            started_at: now_rfc3339(),
        }
    }
}

/// Hash of a value's JSON serialization (settings snapshots, plane geometry).
pub fn hash_json<T: Serialize>(tag: &str, v: &T) -> Hash32 {
    let bytes = serde_json::to_vec(v).unwrap_or_default();
    Hasher::new().tag(tag).bytes(&bytes).finish()
}

/// The current UTC time as RFC 3339 (seconds).
pub fn now_rfc3339() -> String {
    let secs = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs() as i64).unwrap_or(0);
    rfc3339(secs)
}

/// Seconds since the Unix epoch as RFC 3339 UTC (proleptic Gregorian).
pub fn rfc3339(secs: i64) -> String {
    let days = secs.div_euclid(86_400);
    let rem = secs.rem_euclid(86_400);
    // civil-from-days (H. Hinnant)
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };
    format!("{y:04}-{m:02}-{d:02}T{:02}:{:02}:{:02}Z", rem / 3600, rem % 3600 / 60, rem % 60)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rfc3339_known_dates() {
        assert_eq!(rfc3339(0), "1970-01-01T00:00:00Z");
        assert_eq!(rfc3339(951_782_400), "2000-02-29T00:00:00Z");
        assert_eq!(rfc3339(1_793_836_800 + 3_723), "2026-11-05T01:02:03Z");
    }
    #[test]
    fn build_info_embedded() {
        let b = BuildInfo::current();
        assert!(!b.git_commit.is_empty() && !b.build_profile.is_empty());
    }
}
