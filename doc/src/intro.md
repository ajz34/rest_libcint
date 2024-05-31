# Introduction

`rest_libcint` crate is a library for providing gaussian type orbital (GTO) electronic integrals, for audiences in computational chemistry and physics.

## Example

Example of code usage (two-fold symmetry, with shell slices specified, $(\mu \nu | \nabla_r \kappa \lambda)$ with upper-triangular symmetric $\mu, \nu$ in f-contiguous):
```rust
// cint_data: CINTR2CDATA instance
// out: output integral as Vec<f64>
// out_shape: shape of integral, f-contiguous, as Vec<usize>

let shl_slices = vec![[10, 50], [10, 50], [127, 168], [215, 272]];
let (out, out_shape) = cint_data.integral_s2ij::<int2e_ip2>(Some(&shl_slices));
```

## Functionality and Limitations

This library provides several functionalities:
- Wrappers for C library libcint, which is written and maintained on <https://github.com/sunqm/libcint>.
- Generating electronic integrals in parallel on single node.

With some preliminary tests, parallel efficiency of `rest_libcint` should be equally or slightly faster than PySCF wrapper (`pyscf.gto`) to libcint.

Currently some limitations still exists:
- Only gives f-contiguous output tensor. Support of c-contiguous or arbitary strided tensors is not supported.
  API user may handle tensor transpositions properly by themselves.
- Only supports no symmetry (`s1`) or two-fold symmetry (`s2ij`). Higher symmetry will be implemented in future.
  Two-fold symmetry should be useful in resolution-of-identity algorithms, but not enough for conventional integrals.
- Spin-orbital integrators (complex integrals) are not handled properly, and will be implemented in future.

