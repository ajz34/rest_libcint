# Introduction

`rest_libcint` crate is a library for providing gaussian type orbital (GTO) electronic integrals, for audiences in computational chemistry and physics.

## Example

Example of code usage (two-fold symmetry, with shell slices specified, $(\mu \nu | \nabla_r \kappa \lambda)$ with upper-triangular symmetric $\mu, \nu$ in f-contiguous):
```rust
// cint_data:             CINTR2CDATA instance
// out, ecp:              output integral as Vec<f64>
// out_shape, ecp_shape:  shape of integral, f-contiguous, as Vec<usize>

// include all useful functions or traits
use rest_libcint::prelude::*;

let shl_slices = vec![[10, 50], [10, 50], [127, 168], [215, 272]];
let (out, out_shape) = cint_data.integral_s2ij::<int2e_ip2>(Some(&shl_slices));
let (ecp, ecp_shape) = cint_data.ecp_integral_s1::<ECPscalar_ipipnuc>(Some(&shl_slices));
```

This user guide will further clarify code this example in
- [Data Initialization](data-init.md):
  What is `cint_data`? How is it defined?
- [Dimensionality Convention](dimensionality.md):
  How is the integral tensor stored?
  What is the relationship between integral tensor and formula?
- [Integral Options](integral.md):
  What does `integral_s2ij::<int2e_ip2>` exactly do?
  How to perform other kinds of integrals?

Through the entire document, variable naming conventions will also be the same to the above example.

## Functionality and Limitations

This library provides several functionalities:
- Wrappers for C library libcint, which is written and maintained on <https://github.com/sunqm/libcint>.
- Wrappers for ECP code from PySCF ([nr_ecp.c](https://github.com/pyscf/pyscf/blob/v2.5.0/pyscf/lib/gto/nr_ecp.c) and [nr_ecp_deriv.c](https://github.com/pyscf/pyscf/blob/v2.5.0/pyscf/lib/gto/nr_ecp_deriv.c)). 
- Generating electronic integrals in parallel on single node.

With some preliminary tests, parallel efficiency of `rest_libcint` should be equally or slightly faster than PySCF wrapper (`pyscf.gto`) to libcint.

Currently some limitations still exists:
- Only gives F-contiguous output tensor. C-contiguous or arbitary strided tensors is not supported.
  API user may handle tensor transpositions properly by themselves.
- Only supports no symmetry (`s1`) or two-fold symmetry (`s2ij`). Higher symmetry will be implemented in future.
  Two-fold symmetry should be useful in resolution-of-identity algorithms, but not enough for conventional integrals.
- Spin-orbital integrators (complex integrals) are not handled properly, and will be implemented in future.

