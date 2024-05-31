# Integral and Dimensionality

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
    use rest_libcint::prelude::*;
    let shl_slices = vec![[0, 5], [0, 7], [3, 8]];
    let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(Some(&shl_slices));
    ```
    The output shape is $(\mu, \nu, P, t)$, or specifically `[14, 16, 13, 3]`.

    This is also true for other 2, 3, 4-center integrals.

    Also note that by our convention, if number of components is 1 (such as integrals with no derivatives),
    dimension of component will be squeezed out.
    The component dimension will always be the last dimension in a f-contiguous tensor,
    i.e., the most discontiguous dimension.

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

## More on integrals in `rest_libcint`

1. **Optional parameter: `shl_slices`**

    `shl_slices` is optional. By default, all shells are chosen to be integrated:
    ```rust
    use rest_libcint::prelude::*;
    let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(None);
    ```

2. **Symmetry: `s1` and `s2ij`**

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
    
    Note that triangular-packed is upper-triangular in f-contiguous, lower-triangular in c-contiguous.

    <div class="warning">

    For symmetric integrals, we only check whether shell slices are the same (for example of `s2ij`,
    the first and second shell slices must be the same).
    However, we do not check whether the integrator gives the symmetric integral.

    For example, `int3c2e_ip2` is symmetric on the first two atomic orbitals
    $(\mu \nu | \partial_t P) = (\nu \mu | \partial_t P)$, thus applicable to symmetry `s2ij`.
    
    However, `int3c2e_ip1` is not symmetric, since $(\partial_t \mu \nu | P) \neq (\partial_t \nu \mu | P)$.
    Using `int3c2e_ip1` with symmetry `s2ij` is not correct, but the program may not panic in this case.

    </div>

3. **Export as tensor with proper shape**

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

4. **Integrators**

    We use the same naming convention in libcint for integrators. To name a few important integrators,

    - `int2e`: $(\mu \nu | \kappa \lambda)$,
    - `int3c2e`: $(\mu \nu | P)$,
    - `int1e_kin`: $(\mu | - \frac{1}{2} \nabla^2 | \nu)$,

    Some of the integrals are also listed in [PySCF document](https://pyscf.org/pyscf_api_docs/pyscf.gto.html#module-pyscf.gto.moleintor).
    Also see the implemented structures (with trait bound `Integrator`) in `rest_libcint::cint_wrapper`.
