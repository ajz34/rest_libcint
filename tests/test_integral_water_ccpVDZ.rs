use rest_libcint::CINTR2CDATA;
use itertools::Itertools;
use rest_libcint::cint_wrapper::*;

#[cfg(test)]
mod test_c15_tzvp {
    use super::*;

    #[test]
    fn test_contiguous() {

        let mut cint_data = initialize();
        let shl_slices = vec![[0, 5], [0, 7], [3, 8]];
        let (out, out_shape) = cint_data.integral_s1::<int3c2e_ip1>(Some(&shl_slices));
        
        println!("{:?}", out_shape);
    }

    fn initialize() -> CINTR2CDATA {
        // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="6-31G").build()
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
}