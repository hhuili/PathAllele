use std::{env, path::PathBuf, process::Command};

fn main() {
    println!("cargo:rerun-if-changed=native/align_alleles.cpp");
    for key in [
        "PATHALLELE_NATIVE_PREFIX",
        "CONDA_PREFIX",
        "CXX",
        "AR",
        "PATHALLELE_CXX11_ABI",
    ] {
        println!("cargo:rerun-if-env-changed={key}");
    }
    if env::var_os("CARGO_FEATURE_NATIVE_ALIGNER").is_none() {
        return;
    }
    let prefix = PathBuf::from(
        env::var_os("PATHALLELE_NATIVE_PREFIX")
            .or_else(|| env::var_os("CONDA_PREFIX"))
            .expect("native-aligner requires PATHALLELE_NATIVE_PREFIX or CONDA_PREFIX"),
    );
    let prefix = prefix.canonicalize().expect("native dependency prefix");
    let headers = prefix.join("include/wfa2lib");
    assert!(
        prefix.join("include/vcflib/Variant.h").is_file(),
        "Missing vcflib headers"
    );
    assert!(
        headers.join("wavefront/wfa.hpp").is_file(),
        "Missing WFA2 headers"
    );
    let out = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    let object = out.join("align_alleles.o");
    let mut cxx = Command::new(env::var_os("CXX").unwrap_or_else(|| "c++".into()));
    let abi = env::var("PATHALLELE_CXX11_ABI").unwrap_or_else(|_| "1".into());
    assert!(
        abi == "0" || abi == "1",
        "PATHALLELE_CXX11_ABI must be 0 or 1"
    );
    cxx.arg(format!("-D_GLIBCXX_USE_CXX11_ABI={abi}"));
    cxx.args(["-std=c++17", "-O2", "-fPIC", "-c"])
        .arg("native/align_alleles.cpp")
        .arg("-o")
        .arg(&object);
    for include in [
        &prefix.join("include/vcflib"),
        &prefix.join("include"),
        &headers,
    ] {
        cxx.arg("-I").arg(include);
    }
    assert!(
        cxx.status().expect("C++ compiler").success(),
        "Native adapter compilation failed"
    );
    assert!(
        Command::new(env::var_os("AR").unwrap_or_else(|| "ar".into()))
            .arg("crs")
            .arg(out.join("libpathallele_native.a"))
            .arg(object)
            .status()
            .expect("archiver")
            .success(),
        "Native archive failed"
    );
    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=pathallele_native");
    println!(
        "cargo:rustc-link-search=native={}",
        prefix.join("lib").display()
    );
    println!("cargo:rustc-link-lib=dylib=vcflib");
    println!("cargo:rustc-link-lib=dylib=wfa2");
    println!(
        "cargo:rustc-link-lib=dylib={}",
        if env::var("CARGO_CFG_TARGET_OS").unwrap() == "macos" {
            "c++"
        } else {
            "stdc++"
        }
    );
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,{}",
        prefix.join("lib").display()
    );
}
