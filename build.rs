extern crate dunce;
use std::env;

fn main() {
    let blas_dir = env::var("REST_BLAS_DIR").unwrap_or("".to_string());
    let cint_dir = env::var("REST_CINT_DIR").unwrap_or("".to_string());
    
    let library_path = [
        dunce::canonicalize(blas_dir).unwrap(),
        dunce::canonicalize(cint_dir).unwrap(),
    ];
    library_path.iter().for_each(|path| {
        println!("cargo:rustc-link-search=native={}",env::join_paths(&[path]).unwrap().to_str().unwrap())
    });
    println!("cargo:rustc-link-lib=cint");
    println!("cargo:rustc-link-lib=openblas");

    // ecp build
    cc::Build::new()
        .file("src/cecp/nr_ecp.c")
        .file("src/cecp/nr_ecp_deriv.c")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-O3")
        .flag_if_supported("-Wno-implicit-function-declaration")
        .compile("cecp");
    println!("cargo::rerun-if-changed=src/cecp");
}