# Installation

## Install libcint

Currently, this crate does not ship C library libcint by default.
User may need to manually install libcint.
We refer <https://github.com/sunqm/libcint> for more details.

This crate uses C/rust bindings.
Current created binding files are based on libcint [v6.1.2](https://github.com/sunqm/libcint/releases/tag/v6.1.2).
We suggest users to install libcint with version no less than v6.1.2.

## Using `rest_libcint`

Please add following directories into dynamic libaray search path and rust build flags:
- `$CINT_DIR`: libcint library directory, which contains `$CINT_DIR/libcint.so`
- `$BLAS_DIR`: OpenBLAS library directory, which contains `$BLAS_DIR/libopenblas.so`

```bash
export RUSTFLAGS="-L$CINT_DIR -L$BLAS_DIR $RUSTFLAGS"
export LD_LIBRARY_PATH=$CINT_DIR:$BLAS_DIR:$LD_LIBRARY_PATH
```

```admonish note
BLAS library is applied for ECP integrals. 

In our testing, for OpenBLAS, it is better to use user-compiled version > 0.3.20.
The OpenBLAS prebuilt and shipped by Ubuntu may cause some minor efficiency problem,
though not affecting correctness of ECP integral tensors.
```