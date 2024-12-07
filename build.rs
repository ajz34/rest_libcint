use std::{error::Error, path::PathBuf};

/// Generate link search paths from a list of paths.
///
/// This allows paths like `/path/to/lib1:/path/to/lib2` to be split into individual paths.
fn generate_link_search_paths(paths: &[Result<String, impl Error + Clone>]) -> Vec<String> {
    paths
        .iter()
        .map(|path| {
            path.clone()
                .unwrap_or_default()
                .split(":")
                .map(|path| path.to_string())
                .collect::<Vec<_>>()
        })
        .into_iter()
        .flatten()
        .filter(|path| !path.is_empty())
        .collect::<Vec<_>>()
}

/// Check if the library is found in the given paths.
fn check_library_found(
    lib_name: &str,
    lib_paths: &[String],
    lib_extension: &[String],
) -> Option<String> {
    for path in lib_paths {
        for ext in lib_extension {
            let lib_path = PathBuf::from(&path).join(format!("lib{}.{}", lib_name, ext));
            if lib_path.exists() {
                return Some(lib_path.to_string_lossy().to_string());
            }
        }
    }
    return None;
}


fn main() {
    build_cint();
    build_ecp();
}

/// Builds the `libcint` library from the provided path.
fn build_cint() {
    // search dirs
    for key in ["CINT_DIR", "REST_EXT_DIR", "LD_LIBRARY_PATH"].iter() {
        println!("cargo:rerun-if-env-changed={}", key);
    }
    let lib_paths = generate_link_search_paths(&[
        std::env::var("CINT_DIR"),
        std::env::var("REST_EXT_DIR"),
        std::env::var("LD_LIBRARY_PATH"),
    ]);
    if let Some(path) = check_library_found("cint", &lib_paths, &["a".to_string(), "so".to_string()]) {
        let path = std::fs::canonicalize(path).unwrap();
        let path = path.parent().unwrap().display();
        println!("cargo:rustc-link-search=native={}", path);
    } else {
        let dst = cmake::Config::new("external_deps").build();
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    }
    println!("cargo:rustc-link-lib=cint");
    println!("cargo:rustc-link-lib=quadmath");
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
