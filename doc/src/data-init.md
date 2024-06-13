# Data Initialization

To give valid electronic integral, a data structure for describing atomic (molecular) system is required.
In `rest_libcint`, this data structure is `CINTR2CDATA`; variable naming convention of this instance is `cint_data`.
We will provide simple illustrations on this data structure in this document.

## Data Structure

In the following contents, we use cc-pVDZ water as example.

To initialize such kind of data structure, in PySCF, one can execuate following code:

```python
from pyscf import gto
mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="cc-pVDZ").build()
```

What we need to build `CINTR2CDATA` instance is similar (or virtually the same) to PySCF.
In the [last section](#initialize-cint_data-in-rust_libcint) of document, we will see how
`rust_libcint` initializes such data structure. It may seems that code in rust is much more
redundant compared to PySCF; but as a reminder, functionality of this crate is **making
electronic integral**, instead of managing information of molecule or unit cell.
So we provide similar data structure as used by libcint, instead of user-friendly REST or PySCF.

The data is composed by three components: atomic data, shell data, and parameter data.

- Parameter data

    Data block stored in double floats (f64). Meaning of these data are given by atomic data
    and shell data.

    For the specific example of water/cc-pVDZ, parameter data is

    ```python
    array([
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     1.7763,     0.    ,     0.    ,     0.    ,    -0.4448,     0.    ,     1.7198,     0.    , 11720.    ,  1759.    ,   400.8   ,   113.7   ,    37.03  ,    13.27  ,     5.025 ,     1.013 ,
        2.0195,     3.7518,     6.2968,     9.2147,    10.7299,     7.8782,     2.2964,     0.0394,    -0.8949,    -1.7033,    -2.7874,    -4.4459,    -5.2862,    -5.7102,    -1.949 ,     2.7944,     0.3023,     1.03  ,    17.7   ,     3.854 ,
        1.046 ,     6.6386,     5.2543,     2.2875,     0.2753,     0.5818,     1.185 ,     3.5119,    13.01  ,     1.962 ,     0.4446,     0.5798,     0.9834,     1.1193,     0.122 ,     0.5215,     0.727 ,     1.9584])
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
    - `PTR_COORD = 1`: Pointer of coordinate ($x, y, z$) in Bohr (a.u.).
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
    - `ANG_OF = 1`: Angular momentum of shell $l$. This value also indicates how many contracted GTOs
        (CGTOs), or number of atomic orbitals (AOs).
        For example, $l = 2$ meaning that this is $d$ basis function; 6 CGTOs for cartesian or
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
        0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,     0.    ,
        0.    ,     0.    ,     0.    ,     0.    ,     1.7763,     0.    ,     0.    ,     0.    ,    -0.4448,     0.    ,     1.7198,     0.    , 11720.    ,  1759.    ,   400.8   ,   113.7   ,    37.03  ,    13.27  ,     5.025 ,     1.013 ,
        2.0195,     3.7518,     6.2968,     9.2147,    10.7299,     7.8782,     2.2964,     0.0394,    -0.8949,    -1.7033,    -2.7874,    -4.4459,    -5.2862,    -5.7102,    -1.949 ,     2.7944,     0.3023,     1.03  ,    17.7   ,     3.854 ,
        1.046 ,     6.6386,     5.2543,     2.2875,     0.2753,     0.5818,     1.185 ,     3.5119,    13.01  ,     1.962 ,     0.4446,     0.5798,     0.9834,     1.1193,     0.122 ,     0.5215,     0.727 ,     1.9584];

    let c_atm = c_atm.iter().map(|&v| v.to_vec()).collect_vec();
    let c_bas = c_bas.iter().map(|&v| v.to_vec()).collect_vec();
    let mut cint_data = CINTR2CDATA::new();
    cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
    return cint_data;
}
```

## ECP Data Structure

For basis set with ECP, in structure, this requires additional ECP shell data `c_ecp`
to be feed into `CINTR2CDATA`.

In PySCF, ECP shell data is stored in `mol._ecpbas`;
in `rest_libcint`, this is stored in `cint_data.c_ecp`.

Note that adding ECP basis information will also change `cint_data.c_env` in `rest_libcint`
or `mol._env` in PySCF. This crate does not provide functionality to generate molecule
data, so we will not dive into detail on this topic, and refer to code in `pyscf.gto`.

For example of TeH<sub>2</sub>:

```rust
use rest_libcint::CINTR2CDATA;
use itertools::Itertools;

fn initialize() -> CINTR2CDATA {
    // mol = gto.Mole(
    //     atom="Te; H 1 0.94; H 1 0.94 2 104.5",
    //     basis="def2-TZVP", ecp={"Te": "def2-TZVP"}).build()
    let c_atm= vec![
        [24, 20,  4, 23,  0,  0],
        [ 1, 24,  1, 27,  0,  0],
        [ 1, 28,  1, 31,  0,  0],
    ];
    let c_bas = vec![
        [  0,   0,   5,   1,   0,  32,  37,   0],
        [  0,   0,   2,   1,   0,  42,  44,   0],
        [  0,   0,   1,   1,   0,  46,  47,   0],
        [  0,   0,   1,   1,   0,  48,  49,   0],
        [  0,   0,   1,   1,   0,  50,  51,   0],
        [  0,   0,   1,   1,   0,  52,  53,   0],
        [  0,   1,   3,   1,   0,  54,  57,   0],
        [  0,   1,   3,   1,   0,  60,  63,   0],
        [  0,   1,   1,   1,   0,  66,  67,   0],
        [  0,   1,   1,   1,   0,  68,  69,   0],
        [  0,   1,   1,   1,   0,  70,  71,   0],
        [  0,   2,   6,   1,   0,  72,  78,   0],
        [  0,   2,   1,   1,   0,  84,  85,   0],
        [  0,   2,   1,   1,   0,  86,  87,   0],
        [  0,   3,   1,   1,   0,  88,  89,   0],
        [  0,   3,   1,   1,   0,  90,  91,   0],
        [  1,   0,   3,   1,   0,  92,  95,   0],
        [  1,   0,   1,   1,   0,  98,  99,   0],
        [  1,   0,   1,   1,   0, 100, 101,   0],
        [  1,   1,   1,   1,   0, 102, 103,   0],
        [  2,   0,   3,   1,   0,  92,  95,   0],
        [  2,   0,   1,   1,   0,  98,  99,   0],
        [  2,   0,   1,   1,   0, 100, 101,   0],
        [  2,   1,   1,   1,   0, 102, 103,   0],
    ];
    let c_env = vec![
        0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,    0.    ,
        0.    ,    0.    ,    1.7763,    0.    ,    0.    ,    0.    ,   -0.4448,    0.    ,    1.7198,    0.    , 6213.2002,  920.8964,  199.2804,   24.7742,   14.8382,    0.3401,    0.5575,    0.5374,   -1.8554,   20.2697,   12.2788,    6.3808,
       11.501 ,    3.275 ,    2.2228,    4.5993,    1.0776,    2.6721,    0.2814,    0.976 ,    0.1078,    0.4754,  204.294 ,   18.2088,    9.9211,    4.1073,   29.6827,  -63.3954,    3.1442,    1.7221,    0.891 ,    4.5487,    2.7516,    0.5249,
        0.398 ,    0.9224,    0.1654,    0.3077,    0.0651,    0.0959,  121.5106,   32.9688,   19.2499,    4.7198,    2.3428,    1.1135,    8.0596,    8.0036,   -4.4908,    8.6974,    5.4425,    1.3178,    0.492 ,    0.7542,    0.16  ,    0.1056,
        0.3826,    0.227 ,    1.85  ,    7.8731,   34.0613,    5.1236,    1.1647,    0.9062,    1.6355,    2.4145,    0.3272,    1.0931,    0.1031,    0.4596,    0.8   ,    2.2072,   15.2062,   15.2017,  -15.7454,  -20.7424,   16.8145,   15.2062,
       15.2017,    8.7935,  281.0458,   15.7454,   20.7424,   61.6207,   15.2062,   15.2017,   14.8778,   14.2697,    8.7244,    8.2915,   15.7454,   20.7424,   67.4495,  134.9043,   14.6895,   29.4151,   15.2258,   15.2062,   15.205 ,   15.2017,
        6.0718,    5.8048,   53.1357,   15.7454,   35.4321,   20.7424,    9.0698,   13.1223];
    let c_ecp = vec![
        [  0,  -1,   2,   2,   0, 104, 106,   0],
        [  0,   0,   4,   2,   0, 108, 112,   0],
        [  0,   1,   6,   2,   0, 116, 122,   0],
        [  0,   2,   6,   2,   0, 128, 134,   0],
    ];

    let c_atm = c_atm.iter().map(|&v| v.to_vec()).collect_vec();
    let c_bas = c_bas.iter().map(|&v| v.to_vec()).collect_vec();
    let c_ecp = c_ecp.iter().map(|&v| v.to_vec()).collect_vec();
    let mut cint_data = CINTR2CDATA::new();
    cint_data.initial_r2c_with_ecp(
        &c_atm, c_atm.len() as i32,
        &c_bas, c_bas.len() as i32,
        &c_ecp, c_ecp.len() as i32,
        &c_env);
    return cint_data;
}
```

```admonish note title="Note to developer"

Actual ECP computation are performed by struct `rest_libcint::cecp_crafter::ECPData`,
and actual integrator type is `ECPIntegrator` instead of `Integrator`.

Perhaps in the future, we will simply these traits to handle ECP and general integrator
properly. But up to now, these code works.

This should be not of concern for API user.

```
