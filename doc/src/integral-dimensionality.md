# Integral Dimensionality

We will introduce the electronic integral computation utilities, and at the same time clarify tensor dimensionality.
This is important to both audiences from REST and PySCF, if persuit to efficiently calling electronic integral
tensors is of great concern.

<div class="warning">

`rest_libcint` use fortran convention with f-contiguous tensor.

PySCF use C convention; 2-center or 3-center integrals with f-contiguous tensor, 4-center integrals with c-contiguous tensor.

If API callers need different convention or custom contiguous for electronic integral tensor,
they may wish to handle tensor transposition and data copy themselves after calling electronic integral engine.
This is the same to `rest_libcint` or PySCF.

</div>

## Dimensionality clarification

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

In real program,

- In `rest_libcint`, we use f-contiguous in general. For the specific case of `int3c2e_ip1` of water/cc-pVDZ,
    ```rust
    let shl_slices = vec![[0, 5], [0, 7], [3, 8]];
    let (out, out_shape) = cint_data.integral_s1::<int3c2e>(Some(&shl_slices));
    ```
    The output shape is $(\mu, \nu, P, t)$, or specifically `[14, 16, 13, 3]`.

    This is also true for other 2, 3, 4-center integrals.

    Also note that by our convention, if number of components is 1 (such as integrals with no derivatives),
    dimension of component will be squeezed out.

- In PySCF (v2.5), following code could invoke this tensor (for specific shell slices in water/cc-pVDZ):
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
    In fact, to use this tensor c-contiguously (which is default layout in numpy), one should make
    transposition to $(t, P, \nu, \mu)$; or $(\mu, \nu, P, t)$ as f-contiguous:
    ```python
    >>> out.transpose(0, 3, 2, 1).flags
      C_CONTIGUOUS : True
      F_CONTIGUOUS : False
    >>> out.transpose(1, 2, 3, 0).flags
      C_CONTIGUOUS : False
      F_CONTIGUOUS : True
    ```
    This is also true to other kinds of 3-center and 2-center integrals. But also note that,
    4-center integrals are by default c-contiguous; for example $(\partial_t \mu \nu | \kappa \lambda)$
    is c-contiguous for shape $(t, \mu, \nu, \kappa, \lambda)$:
    ```python
    >>> mol.intor("int2e_ip1", shls_slice=(0, 5, 0, 5, 3, 8, 4, 6), aosym="s2ij").flags
      C_CONTIGUOUS : True
      F_CONTIGUOUS : False
    ```
    As a final note, PySCF always let the component dimension (in this specific case, the dimension which refers to $t$) to
    be the most discontiguous dimension.
