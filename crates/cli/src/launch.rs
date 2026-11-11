//! One `foiley` binary is both the command-line tool and the desktop app.
//!
//! * `foiley studio [files]` (alias `gui`) opens Foiley Studio.
//! * `foiley <file>`, an existing file that is not a command name, is short
//!   for `foiley studio <file>`.
//! * Started from the desktop (Windows: double-click, "Open with", a file
//!   dropped on the exe, a shortcut), any start that is not a command opens
//!   Studio.
//! * Everything else is the command line, unchanged.  A start whose stdin is
//!   a pipe or a file never opens a window: with no command it is XFoil.
//!
//! Windows: `foiley.exe` is a console program, so shells wait for it and
//! print its output.  `foiley.manifest` (embedded by `build.rs`) sets
//! `consoleAllocationPolicy = detached`: started from Explorer on Windows 11
//! 24H2 or later it gets no console window at all.  Older Windows creates one
//! for it, and Studio frees it at start-up, so it only flashes.

use std::ffi::OsString;
use std::path::Path;

/// Insert `studio` into `args` (argv, program name first) when this start
/// means the GUI.  `commands` are the subcommand names and aliases.
pub fn route(args: &mut Vec<OsString>, commands: &[String]) {
    if stdin_has_input() {
        return;
    }
    let studio = match args.get(1).map(|a| a.to_string_lossy()) {
        None => from_desktop(),
        Some(a) => !a.starts_with('-') && !commands.iter().any(|c| *c == a) && (from_desktop() || Path::new(&*a).is_file()),
    };
    if studio {
        args.insert(1, "studio".into());
    }
}

/// Run Foiley Studio until its window is closed.
pub fn studio(launch: foiley_gui::Launch) -> Result<(), String> {
    let desktop = from_desktop();
    #[cfg(windows)]
    if desktop {
        win::free_console();
    }
    let r = foiley_gui::run(launch);
    if let Err(e) = &r {
        if desktop {
            foiley_gui::error_dialog(e);
        }
    }
    r
}

/// Started from the desktop, not from a terminal.
#[cfg(windows)]
fn from_desktop() -> bool {
    // 0: no console (the detached policy at work), 1: a console Windows
    // created just for this process.  From a terminal the shell shares it.
    win::console_processes() <= 1
}

#[cfg(not(windows))]
fn from_desktop() -> bool {
    // Linux / macOS desktops start programs through a launcher entry,
    // which can say `foiley studio %F` itself.
    false
}

#[cfg(windows)]
fn stdin_has_input() -> bool {
    win::stdin_is_pipe_or_file()
}

#[cfg(not(windows))]
fn stdin_has_input() -> bool {
    use std::io::IsTerminal;
    !std::io::stdin().is_terminal()
}

#[cfg(windows)]
mod win {
    use std::ffi::c_void;
    use std::os::windows::io::AsRawHandle;

    #[link(name = "kernel32")]
    extern "system" {
        fn GetConsoleProcessList(list: *mut u32, count: u32) -> u32;
        fn FreeConsole() -> i32;
        fn GetFileType(file: *mut c_void) -> u32;
        fn GetStdHandle(which: u32) -> *mut c_void;
        fn SetStdHandle(which: u32, handle: *mut c_void) -> i32;
    }

    const FILE_TYPE_DISK: u32 = 1;
    const FILE_TYPE_CHAR: u32 = 2;
    const FILE_TYPE_PIPE: u32 = 3;
    /// STD_INPUT_HANDLE, STD_OUTPUT_HANDLE, STD_ERROR_HANDLE
    const STD_HANDLES: [u32; 3] = [-10i32 as u32, -11i32 as u32, -12i32 as u32];

    /// Processes attached to this process's console (0: no console).
    pub fn console_processes() -> u32 {
        let mut ids = [0u32; 4];
        // SAFETY: the buffer holds `ids.len()` entries; a larger count is
        // returned, not written.
        unsafe { GetConsoleProcessList(ids.as_mut_ptr(), ids.len() as u32) }
    }

    pub fn stdin_is_pipe_or_file() -> bool {
        let h = std::io::stdin().as_raw_handle();
        // SAFETY: GetFileType accepts any handle value, invalid ones included
        !h.is_null() && matches!(unsafe { GetFileType(h) }, FILE_TYPE_DISK | FILE_TYPE_PIPE)
    }

    /// Leave the console Windows created for this process (it closes), and
    /// clear the standard handles that pointed at it, which is what a
    /// windowed program starts with; redirected handles are kept.
    pub fn free_console() {
        // SAFETY: plain Win32 calls on this process's own standard handles
        unsafe {
            let console: Vec<u32> = STD_HANDLES.into_iter().filter(|&w| GetFileType(GetStdHandle(w)) == FILE_TYPE_CHAR).collect();
            FreeConsole();
            for w in console {
                SetStdHandle(w, std::ptr::null_mut());
            }
        }
    }
}
