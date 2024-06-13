use std::ptr::{null, null_mut};
use std::sync::Mutex;
use itertools::Itertools;
use rayon::prelude::*;
use rayon::{max_num_threads, current_thread_index};
use crate::{cecp::ECPOpt, CintType, CINTR2CDATA};
use crate::cecp::*;
use crate::cecp_wrapper::*;
use crate::cint;
use crate::utilities::*;

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

unsafe impl Send for ECPData {}
unsafe impl Sync for ECPData {}

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

    /* #region cgto size */

    /// Size (number of atomic orbitals) of spherical CGTO at certain basis index.
    /// 
    /// Given `None` input, this function will return the maximum size of CGTO,
    /// which may be helpful for allocating output buffer size of [`Self::integral_block`].
    pub fn cgto_size_sph(&self, id_bas: i32) -> usize {
        use cint::{ANG_OF, NCTR_OF, BAS_SLOTS};
        let loc_ang = ANG_OF as usize;
        let loc_nctr = NCTR_OF as usize;
        let loc_bas = BAS_SLOTS as usize;
        let id_bas = id_bas as usize;
        let nctr = self.c_bas[id_bas * loc_bas + loc_nctr];
        let nang = self.c_bas[id_bas * loc_bas + loc_ang] * 2 + 1;
        return (nctr * nang) as usize;
    }

    /// Size (number of atomic orbitals) of cartesian CGTO at certain basis index.
    /// 
    /// Given `None` input, this function will return the maximum size of CGTO,
    /// which may be helpful for allocating output buffer size of [`Self::integral_block`].
    pub fn cgto_size_cart(&self, id_bas: i32) -> usize {
        use cint::{ANG_OF, NCTR_OF, BAS_SLOTS};
        let loc_ang = ANG_OF as usize;
        let loc_nctr = NCTR_OF as usize;
        let loc_bas = BAS_SLOTS as usize;
        let id_bas = id_bas as usize;
        let nctr = self.c_bas[id_bas * loc_bas + loc_nctr];
        let l = self.c_bas[id_bas * loc_bas + loc_ang];
        let nang = (l + 1) * (l + 2) / 2;
        return (nctr * nang) as usize;
    }

    pub fn cgto_size(&self, id_bas: i32) -> usize {
        return match self.cint_type {
            CintType::Spheric => self.cgto_size_sph(id_bas),
            CintType::Cartesian => self.cgto_size_cart(id_bas),
        }
    }

    /* #endregion */

    /* #region cgto loc */

    /// Location of atomic orbitals, for basis configuration of the current molecule.
    pub fn cgto_loc(&self) -> Vec<usize> {
        let size_vec = (0..self.c_nbas).map(|idx| self.cgto_size(idx)).collect::<Vec<usize>>();
        let mut loc = vec![0; size_vec.len() + 1];
        for idx in (0..size_vec.len() as usize) {
            loc[idx + 1] = loc[idx] + size_vec[idx]; 
        }
        return loc;
    }

    /// Location of atomic orbitals, for specified slice of shell.
    pub fn cgto_loc_slice(&self, shl_slice: &[i32; 2]) -> Vec<usize> {
        let loc = self.cgto_loc();
        return loc[shl_slice[0] as usize .. shl_slice[1] as usize + 1].to_vec();
    }

    /// Location of atomic orbitals, for specified slices of shell.
    pub fn cgto_loc_slices(&self, shl_slices: &Vec<[i32; 2]>) -> Vec<Vec<usize>> {
        return shl_slices.iter().map(|shl_slice| self.cgto_loc_slice(shl_slice)).collect();
    }

    /// Location of atomic orbitals, relative to the first AO index (start index to be 0),
    /// for specified slice of shell.
    pub fn cgto_loc_slice_relative(&self, shl_slice: &[i32; 2]) -> Vec<usize> {
        let loc_slice = self.cgto_loc_slice(shl_slice);
        return loc_slice.iter().map(|x| x - loc_slice[0]).collect();
    }

    /// Location of atomic orbitals, relative to the first AO index (start index to be 0),
    /// for specified slice of shell.
    pub fn cgto_loc_slices_relative(&self, shl_slices: &Vec<[i32; 2]>) -> Vec<Vec<usize>> {
        return shl_slices.iter().map(|shl_slice| self.cgto_loc_slice_relative(shl_slice)).collect();
    }

    /* #endregion */

    /* #region shl_slices sanity check */

    pub fn check_shl_slices<T> (&self, shl_slices: &Vec<[i32; 2]>) -> Result<(), String>
    where
        T: ECPIntegrator
    {
        for shl_slice in shl_slices {
            if shl_slice[1] < shl_slice[0] {
                return Err(format!("in {shl_slice:?}, second number is smaller than first number"));
            }
            if shl_slice[1] > self.c_nbas {
                return Err(format!("in {:?}, number of shells exceeds {:}", shl_slice, self.c_nbas));
            }
            if shl_slice[0] < 0 {
                return Err(format!("in {shl_slice:?}, shell number should not be negative"));
            }
        }
        return Ok(());
    }

    /* #endregion */

    /* #region cgto shape and buffer size */

    /// Shape of integral (in atomic orbital basis), for specified slices of shell.
    pub fn cgto_shape<T> (&self, shl_slices: &Vec<[i32; 2]>) -> Vec<usize>
    where
        T: ECPIntegrator
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        let shape = shl_slices.iter().map(|shl_slice| {
            let loc = self.cgto_loc_slice(shl_slice);
            loc.last().unwrap().clone() - loc.first().unwrap().clone() as usize
        }).collect::<Vec<usize>>();
        return shape;
    }

    /// Obtain cache size for integral.
    /// 
    /// This function should be used with the shell slice one desired.
    /// 
    /// If the shell slice is not known to you currently (molecule information has passed into
    /// `ECPData`), just pass empty `shls_slice = vec![]`, then it should give the maximum
    /// cache size for this molecule/intor.
    pub fn size_of_cache<T> (&mut self, shls_slice: &Vec<[i32; 2]>) -> usize
    where
        T: ECPIntegrator
    {
        let shls_min = shls_slice.iter().map(|x| x[0]).min().unwrap_or(0);
        let shls_max = shls_slice.iter().map(|x| x[1]).max().unwrap_or(self.c_nbas);
        let cache_size = (shls_min..shls_max).into_iter().map(|shl| {
            let mut shls = [shl; 4];
            match self.cint_type {
                CintType::Spheric => unsafe {
                    T::integral_sph(
                        null_mut(), null(), shls.as_mut_ptr(),
                        self.c_atm.as_ptr(), self.c_natm,
                        self.c_bas.as_ptr(), self.c_nbas,
                        self.c_env.as_ptr(), null(), null_mut()) as usize
                    },
                CintType::Cartesian => unsafe {
                    T::integral_cart(
                        null_mut(), null(), shls.as_mut_ptr(),
                        self.c_atm.as_ptr(), self.c_natm,
                        self.c_bas.as_ptr(), self.c_nbas,
                        self.c_env.as_ptr(), null(), null_mut()) as usize
                    },
            }
        }).max().unwrap();
        return cache_size;
    }

    // Obtain buffer size for integral.
    pub fn size_of_buffer<T> (&self, shl_slices: &Vec<[i32; 2]>) -> usize
    where
        T: ECPIntegrator
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        T::n_comp() * shl_slices.iter().map(|&[shl_0, shl_1]| {
            (shl_0..shl_1).map(|shl| self.cgto_size(shl)).max().unwrap() as usize
        }).product::<usize>()
    }

    /* #endregion */

    /// Smallest unit of electron-integral function from libcint.
    /// 
    /// This is not a safe wrapper function, though it is pure rust. Use with caution.
    /// 
    /// * `out` -
    ///     Output integral buffer, need to be allocated enough space before calling this function.
    /// * `shls_slice` -
    ///     shell indices, which size should be n_center; dimension check should be performed
    ///     before calling this function.
    /// * `cache` -
    ///     Cache buffer, need to be allocated enough space before calling this function; simply
    ///     using `vec![]` should also works, which lets libcint manages cache and efficiency decreases.
    ///     See Also [`Self::max_cache_size`] for guide of properly allocate cache.
    pub unsafe fn integral_block<T> (&self, out: &mut [f64], shls: &[i32], shape: &[i32], cache: &mut [f64])
    where
        T: ECPIntegrator
    {
        let cache_ptr = match cache.len() {
            0 => null_mut(),
            _ => cache.as_mut_ptr(),
        };
        let shape_ptr = match shape.len() {
            0 => null_mut(),
            _ => shape.as_ptr(),
        };
        match self.cint_type {
            CintType::Spheric => unsafe {
                T::integral_sph(
                    out.as_mut_ptr(), shape_ptr, shls.as_ptr(),
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(), self.c_opt, cache_ptr)
                },
            CintType::Cartesian => unsafe {
                T::integral_cart(
                    out.as_mut_ptr(), shape_ptr, shls.as_ptr(),
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(), self.c_opt, cache_ptr)
                },
        };
    }

    /// Main integral engine for s1 symmetry.
    /// 
    /// This function a low-level API, which is not intended to be called by user.
    /// This function only works for f-contiguous integral (PySCF convention).
    pub fn integral_s1_inplace<T> (&mut self, out: &mut Vec<f64>, shl_slices: &Vec<[i32; 2]>)
    where
        T: ECPIntegrator
    {
        /* #region 1. dimension definition and sanity check */

        self.check_shl_slices::<T>(shl_slices).unwrap();

        let n_comp = T::n_comp();
        let cgto_shape = self.cgto_shape::<T>(shl_slices);
        let cgto_shape_i32 = cgto_shape.iter().map(|&v| v as i32).collect::<Vec<i32>>();
        let index_shape = shl_slices.iter().map(|[shl_start, shl_stop]| (shl_stop - shl_start) as usize).collect_vec();
        let cgto_locs_rel = self.cgto_loc_slices_relative(shl_slices);

        /* #endregion */

        /* #region 2. preparation for integral engine */

        // optimizer (make integral faster)
        self.optimizer::<T>();

        // cache: thread-local
        let cache_size = self.size_of_cache::<T>(shl_slices);
        let thread_cache = (0..rayon::current_num_threads()).map(|n| {Mutex::new(vec![0.; cache_size])}).collect_vec();

        // out: enable mut vector by passing immut slice
        let out_const_slice = out.as_slice();

        // reverse iteration for f-contiguous
        let index_shape_rev = index_shape.into_iter().rev().collect_vec();
        let shl_slices_rev = shl_slices.iter().rev().collect_vec();
        let cgto_locs_rel_rev = cgto_locs_rel.iter().rev().collect_vec();
        let cgto_shape_rev = cgto_shape.iter().rev().collect_vec();

        /* #endregion */

        /* #region 3. parallel integral generation */

        // Following code of parallel is not fearless.
        // Variable `out` will be written in parallel, which should be considered racing,
        // when calling `self.integral_block::<T>`,
        // and racing would not actually happen if I am careful.

        (0..(index_shape_rev[0] * index_shape_rev[1])).into_par_iter().for_each(|idx_01| {
            let idx_0 = idx_01 / index_shape_rev[1];
            let idx_1 = idx_01 % index_shape_rev[1];
            let shl_0 = idx_0 as i32 + shl_slices_rev[0][0];
            let shl_1 = idx_1 as i32 + shl_slices_rev[1][0];
            let cgto_0 = cgto_locs_rel_rev[0][idx_0];
            let cgto_1 = cgto_locs_rel_rev[1][idx_1];
            
            let thread_index = current_thread_index().unwrap_or(0);
            let mut cache = thread_cache[thread_index].lock().unwrap();
            
            let shls = [shl_1, shl_0];
            let offset = cgto_1 + cgto_shape_rev[1] * cgto_0;

            unsafe {
                let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                self.integral_block::<T>(out_with_offset, &shls, &cgto_shape_i32, &mut cache);
            }
        });
        /* #endregion */

        /* #region 4. cleanup */
        self.optimizer_destruct();
        /* #endregion */
    }

    pub fn integral_s1<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> (Vec<f64>, Vec<usize>)
    where
        T: ECPIntegrator
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices.clone(),
            None => {
                let nbas = self.c_env[AS_ECPBAS_OFFSET as usize] as i32;
                vec![[0, nbas]; 2]
            },
        };
        let mut out_shape = self.cgto_shape::<T>(&shl_slices);
        if T::n_comp() > 1 { out_shape.push(T::n_comp()); }
        let out_size = out_shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(out_size);
        unsafe { out.set_len(out_size) };
        self.integral_s1_inplace::<T>(&mut out, &shl_slices);
        return (out, out_shape);
    }

}

impl CINTR2CDATA {
    pub fn ecp_integral_s1<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> (Vec<f64>, Vec<usize>)
    where
        T: ECPIntegrator
    {
        if self.c_necp == 0 {
            panic!("No ECP data found in the molecule.");
        }
        let mut ecp_data = ECPData::from_cint_data(&self);
        ecp_data.integral_s1::<T>(shl_slices)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_s1() {
        let cint_data = initialize();
        let mut ecp_data = ECPData::from_cint_data(&cint_data);
        let (out, out_shape) = ecp_data.integral_s1::<ECPscalar>(None);
        println!("{:?}", out_shape);
        println!("{:?}", out);
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