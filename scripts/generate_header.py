# # Auto-generate some parts of `libcint` rust interface

# This file is organized to be compatible to jupyter/jupytext, open as notebook.

# - This file must be manually modified (specify `dir_libcint_src`) to be executed.
# - Please let bindgen (`sudo apt install bindgen` or `cargo install bindgen-cli`) in bash `$PATH`.

# ## User Specifications

# -- libcint source code location --
# the same to `git clone https://github.com/sunqm/libcint`
dir_libcint_src = "/home/a/rest_pack/deps/libcint-6.1.2-src/"

# ## Preparasion

# ### python import

import os
import glob
import subprocess
import re

# +
with open(f"{dir_libcint_src}/include/cint.h.in", "r") as f:
    raw_cint_h_in = f.read()

with open(f"{dir_libcint_src}/include/cint_funcs.h", "r") as f:
    raw_cint_funcs_h = f.read()

with open(f"{dir_libcint_src}/CMakeLists.txt", "r") as f:
    raw_cmakelists = f.read()
# -

# ### configure `cint.h`

# -- configure version of cint.h --
ver = {
    "cint_VERSION_MAJOR": None,
    "cint_VERSION_MINOR": None,
    "cint_VERSION_PATCH": None,
}
for line in raw_cmakelists.split("\n"):
    for key in ver:
        if line.startswith("set(" + key):
            ver[key] = int(line.strip().split()[1].replace(")", "").replace("\"", ""))
            break
cint_h = raw_cint_h_in \
    .replace("@cint_VERSION@", f"{ver['cint_VERSION_MAJOR']}.{ver['cint_VERSION_MINOR']}.{ver['cint_VERSION_PATCH']}") \
    .replace("@cint_SOVERSION@", f"{ver['cint_VERSION_MAJOR']}") \
    .replace("#cmakedefine I8", "/* #undef I8 */") \
    .replace("#cmakedefine CACHE_SIZE_I8", "/* #undef CACHE_SIZE_I8 */")

# ### obtain known intors

## -- configure list of integrals known to `cint_funcs.h` --
# known intor of `cint_funcs.h`
known_intor = []
for line in raw_cint_funcs_h.split("\n"):
    if line.startswith("extern CINTOptimizerFunction"):
        known_intor.append(line.strip().replace(";", "").split()[-1].replace("_optimizer", ""))
known_intor = sorted(set(known_intor))

## -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
actual_intor = []
for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
    with open(src_path, "r") as f:
        for line in f.readlines():
            if line.startswith("ALL_CINT(") or line.startswith("ALL_CINT1E("):
                actual_intor.append(line.split("(")[-1].replace(")", "").replace(";", "").strip())
actual_intor = sorted(set(actual_intor))

# assert known intors are all included in all intors
assert len(set(known_intor) - set(actual_intor)) == 0

# ## `bindgen` auto-generation

# ### configure `cint_funcs.h`

# +
# for unknown intors (mostly not in autocode, but in main source code), generate signature for header
token = """
extern CINTOptimizerFunction {0:}_optimizer;
extern CINTIntegralFunction {0:}_cart;
extern CINTIntegralFunction {0:}_sph;
extern CINTIntegralFunction {0:}_spinor;
"""
cint_funcs = raw_cint_funcs_h + "".join([token.format(intor) for intor in sorted(set(actual_intor) - set(known_intor))])    
cint_funcs = cint_funcs.replace("#include <cint.h>", "#include \"cint.h\"")

# libcint2 interface compability
cint_funcs += """
typedef CACHE_SIZE_T CCINTIntegralFunction(double *out, FINT *shls,
                                  FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                                  CINTOpt *opt);
"""
token = """
extern CINTOptimizerFunction c{0:}_optimizer;
extern CCINTIntegralFunction c{0:}_cart;
extern CCINTIntegralFunction c{0:}_sph;
extern CCINTIntegralFunction c{0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([token.format(intor) for intor in sorted(actual_intor)])

# +
with open("cint.h", "w") as f:
    f.write(cint_h)

with open("cint_funcs.h", "w") as f:
    f.write(cint_funcs)
# -

# ### `bindgen` transpile header to rust

subprocess.run([
    "bindgen", "cint_funcs.h",
    "-o", "cint.rs",
    "--no-layout-tests",  # disable redundant tests
    "--wasm-import-module-name=cint",  # to be modified to #[link(name = ...)]
    "--blocklist-function", "[_].*",  # exclude some function of complex.c
    "--blocklist-var", "_.*",  # exclude variables from complex.c and stdint.c
    "--blocklist-type", "_.*",  # exclude types from complex.c
])

# ### modify bindgen auto-generated file for our purpose

with open("cint.rs", "r") as f:
    token = f.read()

# exclude some function of complex.c
res = re.findall(r"(?:\#.*?\})", token.replace("\n", "NEWLINE"))
for r in res:
    if "pub fn c" in r and "pub fn cint" not in r:
        token = token.replace(r.replace("NEWLINE", "\n") + "\n", "")

# change to #[link(name = "cint")]
token = token.replace("wasm_import_module", "name")

with open("cint.rs", "w") as f:
    f.write(token)

# ### cleanup

subprocess.run(["mv", "cint.rs", "../src"])
subprocess.run(["rm", "cint.h", "cint_funcs.h"])

# ## Auto-generate intor_optimizer

token_optimizer = """
use crate::cint;
use crate::CINTR2CDATA;

/// optimizer macro rules
macro_rules! impl_intor_optimizer {
    ($name:expr, $optim_rust:ident, $optim_cint:ident) => {
        pub fn $optim_rust(&mut self){
            self.cint_del_optimizer_rust();
            unsafe {
                cint::$optim_cint(
                    &mut self.c_opt.0, 
                    self.c_atm.0, self.c_natm, 
                    self.c_bas.0, self.c_nbas, 
                    self.c_env.0);
            }
        }
    };
}

impl CINTR2CDATA {
SUBSTITUTE_OPTIMIZER
}
""".strip()

token_sub = "\n".join([
    "    impl_intor_optimizer!(\"INTOR\", INTOR_optimizer_rust, INTOR_optimizer);" \
    .replace("INTOR", intor)
    for intor in actual_intor
])

with open("../src/cint_optimizer_rust.rs", "w") as f:
    f.write(token_optimizer.replace("SUBSTITUTE_OPTIMIZER", token_sub))
