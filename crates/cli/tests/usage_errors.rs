//! Usage errors exit with clap's status 2; `foiley snap` documents status 1.
#![cfg(feature = "snap")]

use std::process::Command;

fn status(args: &[&str]) -> Option<i32> {
    let output = Command::new(env!("CARGO_BIN_EXE_foiley")).args(args).output().unwrap();
    assert!(!output.stderr.is_empty(), "{args:?}: no usage message");
    output.status.code()
}

#[test]
fn usage_errors_exit_with_status_2_and_snap_usage_errors_with_1() {
    assert_eq!(status(&["foil", "polar", "--bogus"]), Some(2));
    // a script file named "snap" does not make this a snap usage error
    assert_eq!(status(&["run", "snap", "--bogus"]), Some(2));
    assert_eq!(status(&["snap", "--bogus"]), Some(1));
    assert_eq!(status(&["--threads", "1", "snap", "--size", "big"]), Some(1));
}
