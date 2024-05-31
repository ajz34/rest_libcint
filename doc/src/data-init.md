# Data Initialization

To give valid electronic integral, a data structure for describing atomic (molecular) system is required.
In `rest_libcint`, this data structure is `CINTR2CDATA`; variable naming convention of this instance is `cint_data`.
We will provide simple illustrations on this data structure in this document.

## Data Structure

In the following contents, we use cc-pVDZ water as example.

```python
from pyscf import gto
mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="cc-pVDZ").build()
```

What we need to build `CINTR2CDATA` instance is similar (or virtually the same) to PySCF.
The data is composed by three components: atomic data, shell data, and parameter data.

- Parameter data

    Data block stored in double floats (f64). Meaning of these data are given by atomic data
    and shell data.

    For the specific example of water/cc-pVDZ, parameter data is

    ```python
    array([    0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
               0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
               0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
               0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
               0.    ,     0.    ,     0.    ,     0.    ,     1.7763,
               0.    ,     0.    ,     0.    ,    -0.4448,     0.    ,
               1.7198,     0.    , 11720.    ,  1759.    ,   400.8   ,
             113.7   ,    37.03  ,    13.27  ,     5.025 ,     1.013 ,
               2.0195,     3.7518,     6.2968,     9.2147,    10.7299,
               7.8782,     2.2964,     0.0394,    -0.8949,    -1.7033,
              -2.7874,    -4.4459,    -5.2862,    -5.7102,    -1.949 ,
               2.7944,     0.3023,     1.03  ,    17.7   ,     3.854 ,
               1.046 ,     6.6386,     5.2543,     2.2875,     0.2753,
               0.5818,     1.185 ,     3.5119,    13.01  ,     1.962 ,
               0.4446,     0.5798,     0.9834,     1.1193,     0.122 ,
               0.5215,     0.727 ,     1.9584])
    ```
    
    In PySCF, parameter data is stored in `mol._env`;
    in `rest_libcint`, this is stored in `cint_data.c_env`.

- Atomic data

    For the specific example of water/cc-pVDZ, atomic data is

    ```python
    array([[ 8, 20,  1, 23,  0,  0],
           [ 1, 24,  1, 27,  0,  0],
           [ 1, 28,  1, 31,  0,  0]], dtype=int32)
    ```

    This is an matrix with shape 3x6, indicating 3 atoms, with 6 data slots for each atom.
    Location of data slots are defined in `cint.h` in libcint:
    
    - `CHARGE_OF = 0`: Atomic charge (8 for oxygen, 1 for hydrogen).
    - `PTR_COORD = 1`: Pointer of coordinate (\\( x, y, z \\)) in Bohr (a.u.).
        For example, starting pointer of coordinate of the second hydrogen is 28, so cartesian
        coordinate of this atom is `mol._env[28:31] = [-0.4448,  0.    ,  1.7198]`.
    - `NUC_MOD_OF = 2`: Nuclear type. `1` refers to point-charge nuclear; 2 refers to gaussian-charge
        nuclear. In general cases, this value should be `1`.
    - `PTR_ZETA = 3`: Manually modified atomic charge distribution by STO function.
        This value is generally set to zero if no manually modification required.
    - `PTR_FRAC_CHARGE = 4`: Fractional atomic charge. This is generally not used.
    - `RESERVE_ATMSLOT = 5`.
    
    In PySCF, parameter data is stored in `mol._atm`;
    in `rest_libcint`, this is stored in `cint_data.c_atm`.

- Shell data

    For the specific example of water/cc-pVDZ, shell data is

    ```python
    array([[ 0,  0,  8,  2,  0, 32, 40,  0],
           [ 0,  0,  1,  1,  0, 56, 57,  0],
           [ 0,  1,  3,  1,  0, 58, 61,  0],
           [ 0,  1,  1,  1,  0, 64, 65,  0],
           [ 0,  2,  1,  1,  0, 66, 67,  0],
           [ 1,  0,  3,  1,  0, 68, 71,  0],
           [ 1,  0,  1,  1,  0, 74, 75,  0],
           [ 1,  1,  1,  1,  0, 76, 77,  0],
           [ 2,  0,  3,  1,  0, 68, 71,  0],
           [ 2,  0,  1,  1,  0, 74, 75,  0],
           [ 2,  1,  1,  1,  0, 76, 77,  0]], dtype=int32)
    ```

    Shell data is the most relavent data to electronic integral tensor computation.

    This is an matrix with shape 11x8, indicating 11 shells, with 8 data slots for each atom.
    Location of data slots are defined in `cint.h` in libcint:
    
    - `ATOM_OF = 0`: This value corresponds to index of relavent atom in atomic data.
    - `ANG_OF = 1`: Angular momentum of shell \\( l \\). This value also indicates how many contracted GTOs
        (CGTOs), or number of atomic orbitals (AOs).
        For example, \\( l = 2 \\) meaning that this is \\( d \\) basis function; 6 CGTOs for cartesian or
        5 CGTOs for spheric.
    - `NCTR_OF = 3`: Number of contracted functions. This value could be larger than one, when
        multiple contracted GTOs utilizes same exponents but different coefficients. In this example
        of water/cc-pVDZ, only the basis function on oxygen atom have two contracted functions.
    - `KAPPA_OF = 4`: This value is related to two-components integrals.
    - `PTR_EXP = 5`: Pointer of exponents to parameter data.
    - `PTR_COEFF = 6`: Pointer of CGTO coefficients to parameter data.
    - `RESERVE_BASLOT = 7`.
    
    In PySCF, parameter data is stored in `mol._bas`;
    in `rest_libcint`, this is stored in `cint_data.c_bas`.

## Initialize `cint_data` in `rust_libcint`

Initialization functions are `CINTR2CDATA::new` and `CINTR2CDATA::initial_r2c`.

Example code:

```rust
use rest_libcint::CINTR2CDATA;
use itertools::Itertools;

fn initialize() -> CINTR2CDATA {
    // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="cc-pVDZ").build()
    let c_atm= vec![
        [ 8, 20,  1, 23,  0,  0],
        [ 1, 24,  1, 27,  0,  0],
        [ 1, 28,  1, 31,  0,  0],
    ];
    let c_bas = vec![
        [ 0,  0,  8,  2,  0, 32, 40,  0],
        [ 0,  0,  1,  1,  0, 56, 57,  0],
        [ 0,  1,  3,  1,  0, 58, 61,  0],
        [ 0,  1,  1,  1,  0, 64, 65,  0],
        [ 0,  2,  1,  1,  0, 66, 67,  0],
        [ 1,  0,  3,  1,  0, 68, 71,  0],
        [ 1,  0,  1,  1,  0, 74, 75,  0],
        [ 1,  1,  1,  1,  0, 76, 77,  0],
        [ 2,  0,  3,  1,  0, 68, 71,  0],
        [ 2,  0,  1,  1,  0, 74, 75,  0],
        [ 2,  1,  1,  1,  0, 76, 77,  0],
    ];
    let c_env = vec![
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     1.7763,
        0.    ,     0.    ,     0.    ,    -0.4448,     0.    ,
        1.7198,     0.    , 11720.    ,  1759.    ,   400.8   ,
        113.7   ,    37.03  ,    13.27  ,     5.025 ,     1.013 ,
        2.0195,     3.7518,     6.2968,     9.2147,    10.7299,
        7.8782,     2.2964,     0.0394,    -0.8949,    -1.7033,
        -2.7874,    -4.4459,    -5.2862,    -5.7102,    -1.949 ,
        2.7944,     0.3023,     1.03  ,    17.7   ,     3.854 ,
        1.046 ,     6.6386,     5.2543,     2.2875,     0.2753,
        0.5818,     1.185 ,     3.5119,    13.01  ,     1.962 ,
        0.4446,     0.5798,     0.9834,     1.1193,     0.122 ,
        0.5215,     0.727 ,     1.9584];

    let c_atm = c_atm.iter().map(|&v| v.to_vec()).collect_vec();
    let c_bas = c_bas.iter().map(|&v| v.to_vec()).collect_vec();
    let mut cint_data = CINTR2CDATA::new();
    cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
    return cint_data;
}
```

