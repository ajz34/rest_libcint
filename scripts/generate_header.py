# # Generate `libcint` Header of Necessary Functions

# - This file must be manually modified (specify `dir_libcint_src`) to be executed.
# - Please let bindgen (`sudo apt install bindgen` or `cargo install bindgen-cli`) in bash `$PATH`.

# ## User Specifications

# -- libcint source code location --
# the same to `git clone https://github.com/sunqm/libcint`
dir_libcint_src = "/home/a/rest_pack/deps/libcint-6.1.2-src/"

# ## Preparasion

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

# ## Configure `cint.h`

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

# ## Configure `cint_funcs.h`

# +
## -- configure list of integrals known to `cint_funcs.h` --
# known intor of `cint_funcs.h`
known_intor = []
for line in raw_cint_funcs_h.split("\n"):
    if line.startswith("extern CINTOptimizerFunction"):
        known_intor.append(line.strip().replace(";", "").split()[-1].replace("_optimizer", ""))
# all intor generated in source code (autocode included)
actual_intor = []
for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
    with open(src_path, "r") as f:
        for line in f.readlines():
            if line.startswith("ALL_CINT(") or line.startswith("ALL_CINT1E("):
                actual_intor.append(line.split("(")[-1].replace(")", "").replace(";", "").strip())

# assert known intors are all included in all intors
assert len(set(known_intor) - set(actual_intor)) == 0
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
# -

# ## Write out header file

# +
with open("cint.h", "w") as f:
    f.write(cint_h)

with open("cint_funcs.h", "w") as f:
    f.write(cint_funcs)
# -

# ## `bindgen` transpile header to rust

subprocess.run([
    "bindgen", "cint_funcs.h",
    "-o", "cint.rs",
    "--no-layout-tests",  # disable redundant tests
    "--wasm-import-module-name=cint",  # to be modified to #[link(name = ...)]
    "--blocklist-function", "[_].*",  # exclude some function of complex.c
    "--blocklist-var", "_.*",  # exclude variables from complex.c and stdint.c
    "--blocklist-type", "_.*",  # exclude types from complex.c
])

# ## modify generated file

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

# ## Cleanup

subprocess.run(["mv", "cint.rs", "../src"])
subprocess.run(["rm", "cint.h", "cint_funcs.h"])
