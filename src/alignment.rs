//! The embedded alignment worker owned by this executable.
use crate::model::Result;
use std::{path::PathBuf, process::Command};

pub struct Adapter {
    pub executable: PathBuf,
}
impl Adapter {
    pub fn command(&self) -> Command {
        let mut command = Command::new(&self.executable);
        command.arg("internal-align");
        command
    }
}

#[cfg(feature = "native-aligner")]
pub fn worker(input: &std::path::Path) -> Result<()> {
    use std::{ffi::CString, os::unix::ffi::OsStrExt};
    extern "C" {
        fn pathallele_align_file(path: *const std::os::raw::c_char) -> std::os::raw::c_int;
    }
    let input = CString::new(input.as_os_str().as_bytes()).map_err(|e| e.to_string())?;
    // C++ owns all alignment objects and catches exceptions. Run this only in
    // an isolated child: upstream vcflib may exit on an alignment error.
    let status = unsafe { pathallele_align_file(input.as_ptr()) };
    if status == 0 {
        Ok(())
    } else {
        Err(format!("Native alignment failed ({status})"))
    }
}

#[cfg(not(feature = "native-aligner"))]
pub fn worker(_input: &std::path::Path) -> Result<()> {
    Err("This build does not include native-aligner".into())
}
