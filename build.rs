use std::env;
use std::path::PathBuf;

fn main() {
    // read environment variables
    let cint_dir = env::var("REST_CINT_DIR");
    let cint_src = env::var("REST_CINT_SRC");
    let cint_ver = env::var("REST_CINT_VER").unwrap_or("6.1.2".to_string());
    let build_temporary = env::var("REST_BUILD_TEMPORARY").unwrap_or("debug".to_string());

    // process sequence:
    // 1. Use `REST_CINT_DIR` as already compiled directory
    // 2. Use `REST_CINT_SRC` as source code directory to be compiled
    // 3. Use `REST_CINT_VER` to download and compile libcint from github
    //    Use `REST_BUILD_TEMPORARY` to save scratch of downloaded or compiled files

    if let Ok(cint_dir) = cint_dir {
        println!("cargo:rustc-link-search=native={cint_dir}");
    } else {
        if let Ok(cint_src) = cint_src {
            build_cint(cint_src);
        } else {
            // download libcint source code if cint_src is not provided
            // also see https://github.com/blas-lapack-rs/openblas-src/blob/master/openblas-build/src/download.rs
            let url = format!("https://github.com/sunqm/libcint/archive/refs/tags/v{cint_ver}.tar.gz");
            let dest = PathBuf::from(build_temporary.clone()).join(format!("libcint-{cint_ver}"));
            let buf = ureq::AgentBuilder::new()
                .tls_connector(std::sync::Arc::new(
                    native_tls::TlsConnector::new().expect("failed to create TLS connector"),
                )).build().get(&url).call().unwrap().into_reader();
            let gz_stream = flate2::read::GzDecoder::new(buf);
            let mut ar = tar::Archive::new(gz_stream);
            ar.unpack(build_temporary.clone()).unwrap();
            assert!(dest.exists());
            let cint_src = dest.to_str().unwrap().to_string();
            build_cint(cint_src);
        }
    }
    println!("cargo:rustc-link-lib=cint");
    println!("cargo:rustc-link-lib=quadmath");

    build_ecp();
}

/// Builds the `libcint` library from the provided path.
fn build_cint(p: String) {
    let dst = cmake::Config::new(p)
        .define("WITH_FORTRAN", "OFF")
        .define("ENABLE_STATIC", "ON")
        .build();
    println!("cargo:rustc-link-search=native={}/lib", dst.display());
}

/// Builds the `cecp` code from PySCF.
fn build_ecp() {
    cc::Build::new()
        .file("src/cecp/f2c_dgemm.c")
        .file("src/cecp/nr_ecp.c")
        .file("src/cecp/nr_ecp_deriv.c")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-implicit-function-declaration")
        .flag_if_supported("-Wno-parentheses")
        .compile("cecp");
    println!("cargo::rerun-if-changed=src/cecp");
}