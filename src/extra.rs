//! Secondary menus and commands (GDES, QDES, MDES, file utilities, dumps).

use crate::cli::Session;

/// Top-level commands not handled in cli.rs. Returns true if handled.
pub fn top_level(s: &mut Session, com: &str, arg: &str) -> bool {
    match com {
        "GDES" => crate::gdes::gdes(s),
        "QDES" => crate::qdes::qdes(s),
        "MDES" => crate::mdes::mdes(s),
        _ => {
            if !crate::opextra::top_level(s, com, arg) {
                crate::outln!(" {:<4} command not recognized.  Type a \"?\" for list", com);
            }
        }
    }
    true
}

/// OPER-level commands not handled in cli.rs. Returns true if handled.
pub fn oper_cmd(s: &mut Session, com: &str, arg: &str) -> bool {
    crate::opextra::oper_cmd(s, com, arg)
}
