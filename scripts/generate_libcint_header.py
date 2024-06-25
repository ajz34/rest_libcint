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

# ### obtain known intors

# -- configure list of integrals known to `cint_funcs.h` --
# known intor of `cint_funcs.h`
known_intor = []
for line in raw_cint_funcs_h.split("\n"):
    if line.startswith("extern CINTOptimizerFunction"):
        known_intor.append(line.strip().replace(";", "").split()[-1].replace("_optimizer", ""))
known_intor = sorted(set(known_intor))

# +
# -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
# pattern: ALL_CINT(<intor>)
# pattern: ALL_CINT1E(<intor>)
# actual_intor = []
# for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
#     with open(src_path, "r") as f:
#         for line in f.readlines():
#             if line.startswith("ALL_CINT(") or line.startswith("ALL_CINT1E("):
#                 actual_intor.append(line.split("(")[-1].replace(")", "").replace(";", "").strip())
# actual_intor = sorted(set(actual_intor))

# +
# -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
# pattern: 
#   CACHE_SIZE_T <intor>(...
#   FINT ng[] = {...};
# pattern rule:
#   <intor> ends with _sph, _cart, _spinor; no # in <intor>

actual_intor = {}
valid_suffix = ("_sph", "_cart", "_spinor")
for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
    with open(src_path, "r") as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if line.startswith("CACHE_SIZE_T int"):
                # match pattern for intor
                intor = line.split("(")[0].split()[-1].strip()
                is_intor_valid = False
                for suffix in valid_suffix:
                    if intor.endswith(suffix):
                        intor = intor.replace(suffix, "")
                        is_intor_valid = True
                if "#" in intor or intor in actual_intor or not is_intor_valid:
                    continue
                # match pattern for ng
                ng = None
                for l in lines[n:n+5]:
                    if l.strip().startswith("FINT ng[]"):
                        ng = [int(i) for i in l.strip().replace(";", "").replace("}", "").split("{")[-1].split(",")]
                if ng is None:
                    raise ValueError(f"{intor}")
                actual_intor[intor] = ng
# -

# assert known intors are all included in all intors
assert len(set(known_intor) - set(actual_intor)) == 0

# ## `bindgen` auto-generation

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

# ### configure `cint_funcs.h`

# necessary header (modified from libcint 6.1.1)
cint_funcs = """
#include "cint.h"

#if !defined HAVE_DEFINED_CINTINTEGRALFUNCTION
#define HAVE_DEFINED_CINTINTEGRALFUNCTION
typedef void CINTOptimizerFunction(
            CINTOpt **opt,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env);
typedef CACHE_SIZE_T CINTIntegralFunctionReal(
            double *out, const FINT *dims, const FINT *shls,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env,
            const CINTOpt *opt, double *cache);
typedef CACHE_SIZE_T CINTIntegralFunctionComplex(
            double complex *out, const FINT *dims, const FINT *shls,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env,
            const CINTOpt *opt, double *cache);
#endif
"""

# for unknown intors (mostly not in autocode, but in main source code), generate signature for header
token = """
extern CINTOptimizerFunction       {0:}_optimizer;
extern CINTIntegralFunctionReal    {0:}_cart;
extern CINTIntegralFunctionReal    {0:}_sph;
extern CINTIntegralFunctionComplex {0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([token.format(intor) for intor in sorted(set(actual_intor))])
cint_funcs = cint_funcs.replace("#include <cint.h>", "#include \"cint.h\"")

# libcint2 interface compability
token = """
extern CINTOptimizerFunction       c{0:}_optimizer;
extern CINTIntegralFunctionReal    c{0:}_cart;
extern CINTIntegralFunctionReal    c{0:}_sph;
extern CINTIntegralFunctionComplex c{0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([token.format(intor) for intor in sorted(actual_intor) if intor not in ["int2e"]])

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

# ## Auto-generate `cint_wrapper.rs`

# ### IntorBase

# +
# initial information

token_wrapper = """
/* Generated by python scripts for libcint low-level wrapper. */
/* Should not modify manually. */

use crate::cint;
use crate::cint::CINTOpt;
use std::os::raw::c_int;
"""

# +
# IntorBase: basic definition

token_wrapper += """
pub trait IntorBase {
    unsafe fn optimizer(
        opt: *mut *mut CINTOpt,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64);
    unsafe fn integral_sph(
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const CINTOpt,
        cache: *mut f64) -> c_int;
    unsafe fn integral_cart(
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const CINTOpt,
        cache: *mut f64) -> c_int;
    unsafe fn integral_spinor(
        out: *mut cint::__BindgenComplex<f64>,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const CINTOpt,
        cache: *mut f64) -> c_int;
    fn n_comp() -> usize;
    fn n_spinor_comp() -> usize;
    fn n_center() -> usize;
    fn ng() -> Vec<i32>;
    fn intor_type() -> &'static str;
    fn name() -> &'static str;
}

macro_rules! impl_intorbase {
    (
        $intor: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $integral_cart: ident,
        $integral_spinor: ident,
        $n_comp: expr,
        $n_spinor_comp: expr,
        $n_center: expr,
        $ng: expr,
        $intor_type: literal,
        $name: literal
    ) => {
#[allow(non_camel_case_types)]
pub struct $intor;
impl IntorBase for $intor {
    unsafe fn optimizer(
            opt: *mut *mut CINTOpt,
            atm: *const c_int,
            natm: c_int,
            bas: *const c_int,
            nbas: c_int,
            env: *const f64) {
        cint::$optimizer(opt, atm, natm, bas, nbas, env)
    }
    unsafe fn integral_sph(
            out: *mut f64,
            dims: *const c_int,
            shls: *const c_int,
            atm: *const c_int,
            natm: c_int,
            bas: *const c_int,
            nbas: c_int,
            env: *const f64,
            opt: *const CINTOpt,
            cache: *mut f64) -> c_int {
        cint::$integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
    }
    unsafe fn integral_cart(
            out: *mut f64,
            dims: *const c_int,
            shls: *const c_int,
            atm: *const c_int,
            natm: c_int,
            bas: *const c_int,
            nbas: c_int,
            env: *const f64,
            opt: *const CINTOpt,
            cache: *mut f64) -> c_int {
        cint::$integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
    }
    unsafe fn integral_spinor(
        out: *mut cint::__BindgenComplex<f64>,
            dims: *const c_int,
            shls: *const c_int,
            atm: *const c_int,
            natm: c_int,
            bas: *const c_int,
            nbas: c_int,
            env: *const f64,
            opt: *const CINTOpt,
            cache: *mut f64) -> c_int {
        cint::$integral_spinor(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
    }
    fn n_comp() -> usize { $n_comp as usize }
    fn n_spinor_comp() -> usize { $n_spinor_comp as usize }
    fn n_center() -> usize { $n_center as usize }
    fn ng() -> Vec<i32> { $ng }
    fn intor_type() -> &'static str { $intor_type }
    fn name() -> &'static str { $name }
}
    };
}
"""


# +
# IntorBase: macro expansion

def gen_impl_intorbase(intor):
    token = """impl_intorbase!(
    {0:},
    {0:}_optimizer,
    {0:}_sph,
    {0:}_cart,
    {0:}_spinor,
    {2:}, {3:}, {4:}, vec!{5:},
    "{1:}", "{0:}");\n"""
    intor_type = intor.split("_")[0]
    intor_center = intor[:5]
    if intor_center == "int2e":
        intor_center = "int4c"
    if intor_center == "int1e":
        intor_center = "int2c"
    n_center = int(intor_center[3])
    ng = actual_intor[intor]
    assert len(ng) == 8
    comp_1e, comp_2e, comp_tensor = ng[-3:]
    comp_all = max(comp_1e, 1) * max(comp_2e, 1) * comp_tensor
    token = token.format(intor, intor_type, comp_tensor, comp_all, n_center, ng)
    return token

for intor in actual_intor:
    token_wrapper += gen_impl_intorbase(intor)
# -

# ### (fallback support) Optimizer

token_optimizer = """
/* Generated by python scripts for optimizer. */
/* Should not modify manually. */
/* The following code will possibly to be deprecated. */

use crate::CINTR2CDATA;

/// optimizer macro rules
macro_rules! impl_intor_optimizer {
    ($name:expr, $optim_rust:ident, $coptim_rust:ident, $optim_cint:ident) => {
        pub fn $optim_rust(&mut self){
            self.cint_del_optimizer_rust();
            unsafe {
                cint::$optim_cint(
                    &mut self.c_opt,
                    self.c_atm.as_mut_ptr(), self.c_natm,
                    self.c_bas.as_mut_ptr(), self.c_nbas,
                    self.c_env.as_mut_ptr());
            }
        }
        pub fn $coptim_rust(&mut self){
            self.cint_del_optimizer_rust();
            unsafe {
                cint::$optim_cint(
                    &mut self.c_opt,
                    self.c_atm.as_mut_ptr(), self.c_natm,
                    self.c_bas.as_mut_ptr(), self.c_nbas,
                    self.c_env.as_mut_ptr());
            }
        }
    };
}

impl CINTR2CDATA {
SUBSTITUTE_OPTIMIZER_FN
}
"""

token_sub = "\n".join([
    "    impl_intor_optimizer!(\"INTOR\", INTOR_optimizer_rust, cINTOR_optimizer_rust, INTOR_optimizer);".replace("INTOR", intor)
    for intor in actual_intor
])
token_optimizer = token_optimizer.replace("SUBSTITUTE_OPTIMIZER_FN", token_sub)

with open("../src/cint_wrapper.rs", "w") as f:
    f.write(token_wrapper + token_optimizer)
