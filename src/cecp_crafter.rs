use std::ptr::null_mut;
use itertools::Itertools;
use crate::{cecp::ECPOpt, CintType, CINTR2CDATA};
use crate::cecp::*;
use crate::cecp_wrapper::*;

#[derive(Clone)]
pub struct ECPData {
    c_opt: *mut ECPOpt,
    c_natm: i32,
    c_nbas: i32,
    cint_type: CintType,
    c_atm: Vec<i32>,
    c_bas: Vec<i32>,
    c_env: Vec<f64>,
}

impl ECPData {

    pub fn new() -> ECPData {
        ECPData {
            c_opt: null_mut(),
            c_natm: 0,
            c_nbas: 0,
            cint_type: CintType::Spheric,
            c_atm: Vec::new(),
            c_bas: Vec::new(),
            c_env: Vec::new(),
        }
    }

    pub fn from_cint_data(cint_data: &CINTR2CDATA) -> ECPData {
        // initialize
        let mut ecp_data = ECPData {
            c_opt: null_mut(),
            c_natm: cint_data.c_natm,
            c_nbas: cint_data.c_nbas,
            cint_type: cint_data.cint_type,
            c_atm: cint_data.c_atm.clone(),
            c_bas: cint_data.c_bas.clone(),
            c_env: cint_data.c_env.clone(),
        };
        // concate ECP data
        ecp_data.c_bas.append(&mut cint_data.c_ecp.clone());
        ecp_data.c_nbas += cint_data.c_necp;
        ecp_data.c_env[AS_ECPBAS_OFFSET as usize] = cint_data.c_nbas as f64;
        ecp_data.c_env[AS_NECPBAS as usize] = cint_data.c_necp as f64;
        return ecp_data;
    }

    pub fn optimizer_destruct(&mut self) {
        unsafe { ECPdel_optimizer(&mut self.c_opt); }
    }

    pub fn optimizer<T> (&mut self)
    where
        T: ECPIntegrator
    {
        self.optimizer_destruct();
        unsafe {
            T::optimizer(
                &mut self.c_opt,
                self.c_atm.as_ptr(), self.c_natm,
                self.c_bas.as_ptr(), self.c_nbas,
                self.c_env.as_ptr());
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_optimizer() {
        let cint_data = initialize();
        ECPData::from_cint_data(&cint_data).optimizer::<ECPscalar_ignuc>();
    }

    fn initialize() -> CINTR2CDATA {
        let c_atm = vec![
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
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            1.77634256e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
           -4.44760657e-01,  0.00000000e+00,  1.71976186e+00,  0.00000000e+00,
            6.21320017e+03,  9.20896400e+02,  1.99280427e+02,  2.47742331e+01,
            1.48381992e+01,  3.40121652e-01,  5.57475265e-01,  5.37379946e-01,
           -1.85537899e+00,  2.02696629e+01,  1.22787620e+01,  6.38078455e+00,
            1.15009935e+01,  3.27495172e+00,  2.22284052e+00,  4.59934781e+00,
            1.07760434e+00,  2.67214195e+00,  2.81366490e-01,  9.76043793e-01,
            1.07815733e-01,  4.75364374e-01,  2.04294009e+02,  1.82087594e+01,
            9.92110243e+00,  4.10733734e+00,  2.96827246e+01, -6.33954414e+01,
            3.14415287e+00,  1.72208840e+00,  8.90989457e-01,  4.54869822e+00,
            2.75157372e+00,  5.24936750e-01,  3.98047196e-01,  9.22364551e-01,
            1.65387852e-01,  3.07690369e-01,  6.50826956e-02,  9.58994923e-02,
            1.21510552e+02,  3.29687944e+01,  1.92498625e+01,  4.71984073e+00,
            2.34280614e+00,  1.11353794e+00,  8.05961882e+00,  8.00355797e+00,
           -4.49082543e+00,  8.69743670e+00,  5.44246301e+00,  1.31784089e+00,
            4.92000615e-01,  7.54170023e-01,  1.60000000e-01,  1.05618344e-01,
            3.82557570e-01,  2.27027341e-01,  1.85000000e+00,  7.87311782e+00,
            3.40613410e+01,  5.12357460e+00,  1.16466260e+00,  9.06184461e-01,
            1.63547849e+00,  2.41451283e+00,  3.27230410e-01,  1.09308835e+00,
            1.03072410e-01,  4.59591351e-01,  8.00000000e-01,  2.20722637e+00,
            1.52061680e+01,  1.52017020e+01, -1.57454500e+01, -2.07424480e+01,
            1.68144730e+01,  1.52061680e+01,  1.52017020e+01,  8.79352600e+00,
            2.81045843e+02,  1.57454500e+01,  2.07424480e+01,  6.16206560e+01,
            1.52061680e+01,  1.52017020e+01,  1.48778010e+01,  1.42697310e+01,
            8.72443500e+00,  8.29151500e+00,  1.57454500e+01,  2.07424480e+01,
            6.74494640e+01,  1.34904304e+02,  1.46895470e+01,  2.94150630e+01,
            1.52258480e+01,  1.52061680e+01,  1.52050080e+01,  1.52017020e+01,
            6.07176900e+00,  5.80476000e+00,  5.31356870e+01,  1.57454500e+01,
            3.54320570e+01,  2.07424480e+01,  9.06980200e+00,  1.31223040e+01,
        ];
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
}
