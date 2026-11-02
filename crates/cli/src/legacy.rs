//! XFoil-compatible menu mode: reads XFoil commands from stdin, exactly like
//! the original program (existing batch scripts run unchanged).

use foiley::cli::{Input, Session};

pub fn run(file: Option<String>, verbose: bool) {
    if verbose {
        foiley::set_verbose(true);
    }
    let mut s = Session::new(Input::stdin());
    s.run(file);
    foiley::flush_out();
}
