# Dimensionality Convention

We will clarify tensor dimensionality and contiguous conventions.
We will demonstrate some code to perform integral computation here, which will be further elaborated in
[Integral Options](integral.md).

This is important to both audiences from REST and PySCF, if persuit to efficiently calling electronic integral
tensors is of great concern.

```admonish warning

`rest_libcint` use F-contiguous convention universally.

For PySCF, the shape of output tensor is generally not accordance to contiguous; 2-center or 3-center integrals with F-contiguous tensor, 4-center integrals with C-contiguous tensor.

If API callers need different convention or custom contiguous for electronic integral tensor,
they may wish to handle tensor transposition and data copy themselves after calling electronic integral engine.
This is the same to `rest_libcint` or PySCF.

```

## Definition

For a tensor $A_{ijk}$, if its shape is **defined** as $(i, j, k)$, then
- C-contiguous: $A_{ijk} = A[(i \times \mathrm{dim}_j + j) \times \mathrm{dim}_k + k]$, index from $k$;
- F-contiguous: $A_{ijk} = A[i + \mathrm{dim}_i \times (j + \mathrm{dim}_j \times k)]$, index from $i$.

As a result, if $A_{ijk} = B_{kji}$, and

- shape of $A_{ijk}$ is defined as $(i, j, k)$, stored in C-contiguous,
- shape of $B_{kji}$ is defined as $(k, j, i)$, stored in F-contiguous,

then $A_{ijk}$ and $B_{kji}$ share the exact same data layout in memory.

## Integral tensor

Taking derivative to 3c-2e (three-center, two-electron) ERI (electron-repulsion integral) as example. The
integral of concern is `int3c2e_ip1`
$$
    (\nabla \mu \nu | P)
$$
In our notation convention, if variable $t \in \{x, y, z\}$ indicates component of electron coordinate
$\bm{r}$, then this integral could be written as
$$
    (\partial_t \mu \nu | P) = \iint \frac{\partial \phi_\mu (\bm{r})}{\partial t} \phi_\mu (\bm{r}) \frac{1}{| \bm{r} - \bm{r}' |} \phi_P (\bm{r}') \, \mathrm{d} \bm{r} \, \mathrm{d} \bm{r}'
$$
By convention from PySCF, we name
- component: $t$,
- `i, j, k`: $\mu, \nu, P$

Note that up to now, we have not defined the shape of output tensor.
Shape of output tensor can be different for different programs with different conventions.

## Convention of `rest_libcint`

In `rest_libcint`, we use F-contiguous in general. For the specific case of `int3c2e_ip1` of water/cc-pVDZ ([initialization code](data-init.md#initialize-cint_data-in-rust_libcint))
```rust
use rest_libcint::prelude::*;
let shl_slices = vec![[0, 5], [0, 7], [3, 8]];
let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(Some(&shl_slices));
```
The `out_shape` is $(\mu, \nu, P, t)$, or specifically `[14, 16, 13, 3]`.

This is also true for other 2, 3, 4-center integrals.

Also note that by our convention, if number of components is 1 (such as integrals with no derivatives),
dimension of component will be squeezed out.
The component dimension will always be the last dimension in a F-contiguous tensor,
i.e., the most discontiguous dimension.

## Export as tensor with proper shape

Currently, this crate gives output electronic integral tensor by 1-dimensional vector (`Vec<f64>`) with
f-contiguous shape (`Vec<usize>`). We do not assume or recommend any tensor backends. User may handle
output tensors by feeding integral data `out` and shape `out_shape` to their favourate tensor backends.

For example of ndarray (note `.f()` means the input shape is f-contiguous):
```rust
use rest_libcint::prelude::*;
use ndarray::prelude::*;
let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(None);
let out = Array::from_shape_vec(out_shape.f(), out).unwrap();
```

## Convention of PySCF

In PySCF (v2.5), following code could invoke this tensor (for specific shell slices in water/cc-pVDZ):
```python
>>> mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="cc-pVDZ").build()
>>> out = mol.intor("int3c2e_ip1", shls_slice=(0, 5, 0, 7, 3, 8))
>>> out.shape
(3, 14, 16, 13)
```
Shape of this tensor corresponds indices $(t, \mu, \nu, P)$. But as an important note,
**this tensor is not contiguous by current layout**:
```python
>>> out.flags
C_CONTIGUOUS : False
F_CONTIGUOUS : False
```
In fact, to use this tensor C-contiguously (which is default layout in numpy), one should make
transposition to $(t, P, \nu, \mu)$ as C-contiguous; or $(\mu, \nu, P, t)$ as F-contiguous:
```python
>>> out.transpose(0, 3, 2, 1).flags
C_CONTIGUOUS : True
F_CONTIGUOUS : False
>>> out.transpose(1, 2, 3, 0).flags
C_CONTIGUOUS : False
F_CONTIGUOUS : True
```
This is also true to other kinds of 3-center and 2-center integrals. But also note that,
4-center integrals are by default C-contiguous; for example $(\partial_t \mu \nu | \kappa \lambda)$
is C-contiguous for shape $(t, \mu, \nu, \kappa, \lambda)$:
```python
>>> mol.intor("int2e_ip1", shls_slice=(0, 5, 0, 5, 3, 8, 4, 6), aosym="s2ij").flags
C_CONTIGUOUS : True
F_CONTIGUOUS : False
```
As a final note, PySCF always let the component dimension (in this specific case, the dimension which refers to $t$) to
be the most discontiguous dimension.

## Relation of `rest_libcint` and PySCF

For 3-center and 2-center integrals, `rest_libcint` and PySCF shares the same data layout in memory.

More specifically, for integral $(\partial_t \mu \nu | P)$,

| program | shape | contiguous |
|--|--|--|
| rest_libcint | $(\mu, \nu, P, t)$ | F-contiguous |
| PySCF | $(t, \mu, \nu, P)$ | - |
| PySCF (transposed) | $(t, P, \nu, \mu)$ | C-contiguous |

For 4-center integrals, `rest_libcint` and PySCF does not share the same data layout in memory.
