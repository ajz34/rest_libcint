# Installation

By default, [Libcint](https://github.com/sunqm/libcint) is downloaded through internet and compiled as static library.

To install this crate, [quadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) (Quad-Precision Math Library) is required currently.
`quadmath` is generally shipped with gcc.

## Recommanded installation

If there's internet access to github, then simply change to `rest_libcint` project root directory, and then simply
```bash
$ cargo build
```

## Use locally stored `libcint` source code or library

If internet access to github is blocked, then there possibly two ways:

1. Use locally stored `libcint` source code. In this case, declare environment variable
    ```bash
    $ export REST_CINT_SRC=<path of libcint-6.x.x source code>
    ```
    `REST_CINT_SRC` should be a path containing `CMakeLists.txt`. Then simply compile by
    ```bash
    $ cargo build
    ```
    Through this way, a temporary directory `debug` will be generated,
    containing libcint source code downloaded from internet. You can declare
    `REST_BUILD_TEMPORARY` to change this behavior.

2. Use locally compiled `libcint` library. In this case, declare environment variable
    ```bash
    $ export REST_CINT_DIR=<path of libcint/lib>
    ```
    `REST_CINT_DIR` should a path containing either `libcint.so` or `libcint.a`. Then compile by
    ```bash
    $ cargo build --no-default-features
    ```
    If `rest_libcint` is used by other libraries, then dependency declaration in `Cargo.toml` should be something like
    ```toml
    [dependencies]
    rest_libcint = { version = "*", default-features = false }
    ```
    If dynamic library (`libcint.so`) is actually linked, please also append `REST_CINT_DIR` to `LD_LIBRARY_PATH` at runtime.

    For more details of installing `libcint`, we refer to <https://github.com/sunqm/libcint>.

    This crate uses C/rust bindings.
    Current created binding files are based on libcint [v6.1.2](https://github.com/sunqm/libcint/releases/tag/v6.1.2).
    We suggest users to install libcint with version no less than v6.1.2.

```admonish note

As a notice, if there's any problem in installation, one can always specify `RUSTFLAGS`
and `LD_LIBRARY_PATH` environment variables to manually link some libraries.

For more details of installation, see `build.rs`.

```