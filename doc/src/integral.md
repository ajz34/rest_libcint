# Integral Options

## Optional parameter: `shl_slices`

`shl_slices` is optional. By default, all shells are chosen to be integrated.

For example of `int3c2e_ip1`, which is $(\partial_t \mu \nu | P)$,
```rust
let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(None);
```
To specify certain shells to be integrated, one may call integrator by
```rust
let shl_slices = vec![[0, 5], [0, 7], [3, 8]];
let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(Some(&shl_slices));
```

## Symmetry: `s1` and `s2ij`

Currently, we have implemented no symmetry (`s1`), and symmetric in first two atomic orbitals (`s2ij`).
For the current water/cc-pVDZ case of `int3c2e_ip2` $(\mu \nu | \partial_t P)$,
- shape of `s1`: $(\mu, \nu, P, t)$, which is (24, 24, 24, 3).
    ```rust
    let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip2>(None);
    ```
- shape of `s2ij`: $(\mathrm{tp}(\mu \nu), P, t)$, which is (300, 24, 3),
    where $\mathrm{tp}$ means triangular-packed.
    ```rust
    let (out, out_shape) = cint_data.integral_s2ij::<int3c2e_ip2>(None);
    ```

Note that triangular-packed is upper-triangular in F-contiguous, or lower-triangular in C-contiguous.

```admonish warning

For symmetric integrals, we only check whether shell slices are the same (for example of `s2ij`,
the first and second shell slices must be the same).
However, we do not check whether the integrator gives the symmetric integral.

For example, `int3c2e_ip2` is symmetric on the first two atomic orbitals
$(\mu \nu | \partial_t P) = (\nu \mu | \partial_t P)$, thus applicable to symmetry `s2ij`.

However, `int3c2e_ip1` is not symmetric, since $(\partial_t \mu \nu | P) \neq (\partial_t \nu \mu | P)$.
Using `int3c2e_ip1` with symmetry `s2ij` is not correct, but the program may not panic in this case.

```

## Integrators

We use the same naming convention in libcint for integrators. To name a few important integrators,

- `int2e`: $(\mu \nu | \kappa \lambda)$,
- `int3c2e`: $(\mu \nu | P)$,
- `int1e_kin`: $(\mu | - \frac{1}{2} \nabla^2 | \nu)$,
- `int2e_ip2`: $(\mu \nu | \partial_t \kappa \lambda)$

Some of the integrals are also listed in [PySCF document](https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.moleintor), and we refer to that document for more details.
Also see the implemented structures (with trait bound `Integrator`) in `rest_libcint::cint_wrapper`.

## ECP integrals

ECP integral code is quite similar to that of general integral with some differences.

ECP integral is performed by `integral_ecp_s1` function. For example of `ECPscalar` type:

```rust
let (out, shape) = cint_data.integral_ecp_s1::<ECPscalar>(None);
```

Currently only `s1` symmetry is available.

ECP Integrators are implemented in `rest_libcint::cecp_wrapper` (with trait bound `ECPIntegrator`).
To name a few important integrators,
- `ECPscalar`
- `ECPscalar_iprinvip`

## Spinor Integration

Spinor integral code is also similar to that of general integral with some differences.

Spinor integral is performed by `integral_spinor_s1` function. For example of `int2e_ip1ip2` type:

```rust
let (out, shape) = cint_data.integral_spinor_s1::<int2e_ip1ip2>(Some(&shl_slices));
```

Please note that the output data is `Complex<f64>` instead of double float `f64`.
Definition of complex type comes from crate `num-complex`.

Currently, for high-level API, only `s1` symmetry is available (`integral_spinor_s1`).
We expose functions for other kind of symmetries for spinor integral, though it is not desired in general.
