# Installation

## Install libcint

Currently, this crate does not ship C library libcint by default.
User may need to manually install libcint.
We refer <https://github.com/sunqm/libcint> for more details.

This crate uses C/rust bindings.
Current created binding files are based on libcint [v6.1.2](https://github.com/sunqm/libcint/releases/tag/v6.1.2).
We suggest users to install libcint with version no less than v6.1.2.

## Using `rest_libcint`

After installation of libcint, one obtains a directory with shared-library file `$CINT_DIR/libcint.so`.
Please add `$CINT_DIR` into dynamic libaray search path, and rust build flags:
```bash
export RUSTFLAGS="-L$CINT_DIR" $RUSTFLAGS
export LD_LIBRARY_PATH=$CINT_DIR:$LD_LIBRARY_PATH
```
