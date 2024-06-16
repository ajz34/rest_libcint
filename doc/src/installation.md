# Installation

By default, [Libcint](https://github.com/sunqm/libcint) is downloaded through internet and compiled as static library.

To install this crate, [quadmath](https://gcc.gnu.org/onlinedocs/libquadmath/) (Quad-Precision Math Library) is required currently.
`quadmath` is generally shipped with gcc.

## Recommanded installation

If there's internet access to github, then then simply
```bash
$ cargo build
```

## Advanced options

Advanced options are useful if internet access to github is blocked, or if you are developer of REST.

Building process can be controlled by bash environment variables.
The following building workflows are listed by priority;
i.e., if environment variable `REST_CINT_DIR` is defined (for 1st workflow),
then value of `REST_CINT_SRC` (for 2nd workflow) or `REST_CINT_VER` (for 3rd workflow) will be ignored.

1. Use locally compiled `libcint` library. In this case, declare environment variable
    ```bash
    $ export REST_CINT_DIR=<path of libcint/lib>
    ```
    `REST_CINT_DIR` should a path containing either `libcint.so` or `libcint.a`.
    Then simply compile by `cargo build`.

    If dynamic library (`libcint.so`) is actually linked, please also append `REST_CINT_DIR` to `LD_LIBRARY_PATH` at runtime.

    For more details of installing `libcint`, we refer to <https://github.com/sunqm/libcint>.

    This crate uses C/rust bindings.
    Current created binding files are based on libcint [v6.1.2](https://github.com/sunqm/libcint/releases/tag/v6.1.2).
    We suggest users to install libcint with version no less than v6.1.2.

2. Use locally stored `libcint` source code. In this case, declare environment variable
    ```bash
    $ export REST_CINT_SRC=<path of libcint-6.x.x source code>
    ```
    `REST_CINT_SRC` should be a path containing `CMakeLists.txt`.
    Then simply compile by `cargo build`.

    Through this way, a temporary directory `debug` will be generated,
    containing libcint source code downloaded from internet. You can declare
    `REST_BUILD_TEMPORARY` to change this behavior.

3. Download from internet for specified version of libcint.
    This can be controlled by environment variable `REST_CINT_VER`:
    ```bash
    $ export REST_CINT_VER="6.1.2"
    ```

```admonish note

As a notice, if there's any problem in installation, one can always specify `RUSTFLAGS`
and `LD_LIBRARY_PATH` environment variables to manually link some libraries.

For more details of installation, see `build.rs`.

```