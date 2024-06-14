use std::env;
use std::path::PathBuf;

fn main() {
    // read environment variables
    let cint_dir = env::var("REST_CINT_DIR").unwrap_or("".to_string());
    let build_temporary = env::var("REST_BUILD_TEMPORARY").unwrap_or("debug".to_string());
    let mut cint_src = env::var("REST_CINT_SRC").unwrap_or("".to_string());
    let is_build_cint = feature_enabled("build_cint");

    if is_build_cint {
        const LIBCINT_VERSION: &str = "6.1.2";
        if cint_src.is_empty() {
            // download libcint source code if cint_src is not provided
            // also see https://github.com/blas-lapack-rs/openblas-src/blob/master/openblas-build/src/download.rs
            let url = format!("https://github.com/sunqm/libcint/archive/refs/tags/v{LIBCINT_VERSION}.tar.gz");
            let dest = PathBuf::from(build_temporary.clone()).join(format!("libcint-{LIBCINT_VERSION}"));
            let buf = ureq::AgentBuilder::new()
                .tls_connector(std::sync::Arc::new(
                    native_tls::TlsConnector::new().expect("failed to create TLS connector"),
                )).build().get(&url).call().unwrap().into_reader();
            let gz_stream = flate2::read::GzDecoder::new(buf);
            let mut ar = tar::Archive::new(gz_stream);
            ar.unpack(build_temporary.clone()).unwrap();
            assert!(dest.exists());
            cint_src = dest.to_str().unwrap().to_string();
        }
        // build cint
        build_cint(cint_src);
    } else if !cint_dir.is_empty() {
        // link cint by specifing directory containing library
        println!("cargo:rustc-link-search=native={cint_dir}");
    }
    println!("cargo:rustc-link-lib=cint");
    println!("cargo:rustc-link-lib=quadmath");
    println!("cargo:rustc-link-lib=gomp");

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

/// Checks if a specific feature is enabled in the current build.
/// 
/// Features are the ones defined in `Cargo.toml`.
fn feature_enabled(feature: &str) -> bool {
    env::var(format!("CARGO_FEATURE_{}", feature.to_uppercase())).is_ok()
}