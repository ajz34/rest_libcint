use core::panic;
use std::slice::from_raw_parts_mut;
use std::{ptr::null, ptr::null_mut};
use itertools::{Format, Itertools};
use rayon::prelude::*;
use rayon::{max_num_threads, current_thread_index};
use crate::cint_wrapper::*;
use crate::cint;
use crate::{CintType, CINTR2CDATA};

unsafe impl Send for CINTR2CDATA {}
unsafe impl Sync for CINTR2CDATA {}

unsafe fn cast_mut_slice<T> (slc: &[T]) -> &mut [T] {
    let len = slc.len();
    let ptr = slc.as_ptr() as *mut T;
    let mut mut_vector = from_raw_parts_mut(ptr, len);
    return mut_vector;
}

/* #region indices computation and copy */

#[inline(always)]
fn get_f_index_4d(indices: &[usize; 4], shape: &[usize; 4]) -> usize {
    return indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2] + shape[2] * indices[3]));
}

#[inline(always)]
fn get_f_index_4d_s2ij(indices: &[usize; 4], shape: &[usize; 3]) -> usize {
    return indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2] + shape[1] * indices[3]);
}

#[inline(always)]
fn copy_4d_s2ij_offdiag<T> (out: &mut [T], buf: &[T], out_offsets: &[usize; 4], buf_shape: &[usize; 4], out_s2ij_shape: &[usize; 3])
where
    T: Copy
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..buf_shape[0] {
                    let out_index = get_f_index_4d_s2ij(&[out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + k, out_offsets[3]], out_s2ij_shape);
                    let buf_index = get_f_index_4d(&[i, j, k, c], buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

#[inline(always)]
fn copy_4d_s2ij_diag<T> (out: &mut [T], buf: &[T], out_offsets: &[usize; 4], buf_shape: &[usize; 4], out_s2ij_shape: &[usize; 3])
where
    T: Copy
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..(j + 1) {
                    let out_index = get_f_index_4d_s2ij(&[out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + k, out_offsets[3]], out_s2ij_shape);
                    let buf_index = get_f_index_4d(&[i, j, k, c], buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

/* #endregion */

impl CINTR2CDATA {
    /// Remove optimizer
    pub fn optimizer_destruct(&mut self) {
        unsafe { cint::CINTdel_optimizer(&mut self.c_opt); }
    }

    /// Optimizer of libcint intors.
    /// 
    /// To use optimizer, one need to take intor information (such as `int2e_ip1`) into `optimizer::<T>`.
    /// An example could be (by properly defining `c_atm`, `c_bas`, `c_env`)
    /// 
    /// ```
    /// let mut cint_data = CINTR2CDATA::new();
    /// cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
    /// cint_data.optimizer::<int2e_ip1>();
    /// ```
    pub fn optimizer<T> (&mut self)
    where
        T: IntorBase
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

    /// Range of atomic orbitals, for specified slice of shell.
    pub fn cgto_range_slice(&self, shl_slice: &[i32; 2]) -> std::ops::Range<usize> {
        let loc = self.cgto_loc_slice(shl_slice);
        return loc.first().unwrap().clone()..loc.last().unwrap().clone();
    }

    /// Range of atomic orbitals, for specified slices of shell.
    pub fn cgto_range_slices(&self, shl_slices: &Vec<[i32; 2]>) -> Vec<std::ops::Range<usize>> {
        return shl_slices.iter().map(|shl_slice| self.cgto_range_slice(shl_slice)).collect();
    }

    pub fn check_shl_slices<T> (&self, shl_slices: &Vec<[i32; 2]>) -> Result<(), String>
    where
        T: IntorBase
    {
        let n_center = T::n_center();
        if shl_slices.len() != n_center {
            return Err(format!("number of centers {n_center} is not the same to `shl_slices` {shl_slices:?}"));
        }
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

    /// Shape of integral (in atomic orbital basis), for specified slices of shell.
    pub fn cgto_shape<T> (&self, shl_slices: &Vec<[i32; 2]>) -> Vec<usize>
    where
        T: IntorBase
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        let shape = shl_slices.iter().map(|shl_slice| {
            let loc = self.cgto_loc_slice(shl_slice);
            loc.last().unwrap().clone() - loc.first().unwrap().clone() as usize
        }).collect::<Vec<usize>>();
        return shape;
    }

    pub fn cgto_shape_s2ij<T> (&self, shl_slices: &Vec<[i32; 2]>) -> Result<Vec<usize>, String>
    where
        T: IntorBase
    {
        let n_center = T::n_center();
        self.check_shl_slices::<T>(shl_slices)?;
        if (shl_slices[0] != shl_slices[1]) {
            return Err(format!(
                "for symmetry s2ij, last two elements of `shl_slices` {:?}, {:?} should be the same",
                shl_slices[0], shl_slices[1]));
        }
        let mut shape: Vec<usize> = vec![];
        {
            let shl_slice = shl_slices[0];
            let loc = self.cgto_loc_slice(&shl_slice);
            let d = (loc.last().unwrap().clone() - loc.first().unwrap().clone()) as usize;
            shape.push(d * (d + 1) / 2);
        }
        for n in 2..n_center {
            let shl_slice = shl_slices[n];
            let loc = self.cgto_loc_slice(&shl_slice);
            let d = (loc.last().unwrap().clone() - loc.first().unwrap().clone()) as usize;
            shape.push(d);
        }
        return Ok(shape);
    }

    /// Obtain cache size for integral.
    /// 
    /// This function should be used with the shell slice one desired.
    /// 
    /// If the shell slice is not known to you currently (molecule information has passed into
    /// `CINTR2CDATA`), just pass empty `shls_slice = vec![]`, then it should give the maximum
    /// cache size for this molecule/intor.
    /// 
    /// ```no_run
    /// let mut cint_data = CINTR2CDATA::new();
    /// cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
    /// let shls_slice = vec![[0, 2], [0, 1], [1, 3], [0, 2]];
    /// // maximum cache for the given shells
    /// println!("{:?}", cint_data.max_cache_size::<int2e>(&shls_slice));
    /// // maximum cache for the whole molecule
    /// println!("{:?}", cint_data.max_cache_size::<int2e>(&vec![]));
    /// ```
    pub fn size_of_cache<T> (&mut self, shls_slice: &Vec<[i32; 2]>) -> usize
    where
        T: IntorBase
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
        T: IntorBase
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        T::n_comp() * shl_slices.iter().map(|&[shl_0, shl_1]| {
            (shl_0..shl_1).map(|shl| self.cgto_size(shl)).max().unwrap() as usize
        }).product::<usize>()
    }

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
        T: IntorBase
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
        T: IntorBase
    {
        /* #region 1. dimension definition and sanity check */

        self.check_shl_slices::<T>(shl_slices).unwrap();

        let n_comp = T::n_comp();
        let n_center = T::n_center();
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
        let thread_cache: Vec<Vec<f64>> = vec![vec![0.; cache_size]; rayon::current_num_threads()];

        // out: enable mut vector by passing immut slice
        let out_const_slice = out.as_slice();

        // reverse iteration for f-contiguous
        let index_shape_rev = index_shape.into_iter().rev().collect_vec();
        let shl_slices_rev = shl_slices.iter().rev().collect_vec();
        let cgto_locs_rel_rev = cgto_locs_rel.iter().rev().collect_vec();
        let cgto_shape_rev = cgto_shape.iter().rev().collect_vec();

        /* #endregion */

        /* #region 3. parallel integral generation */

        (0..(index_shape_rev[0] * index_shape_rev[1])).into_par_iter().for_each(|idx_01| {
            let idx_0 = idx_01 / index_shape_rev[1];
            let idx_1 = idx_01 % index_shape_rev[1];
            let shl_0 = idx_0 as i32 + shl_slices_rev[0][0];
            let shl_1 = idx_1 as i32 + shl_slices_rev[1][0];
            let cgto_0 = cgto_locs_rel_rev[0][idx_0];
            let cgto_1 = cgto_locs_rel_rev[1][idx_1];
            
            let thread_index = current_thread_index().unwrap_or(0);
            let mut cache = unsafe { cast_mut_slice(&thread_cache[thread_index]) };
            
            match n_center {
                2 =>
                {
                    let shls = [shl_1, shl_0];
                    let offset = cgto_1 + cgto_shape_rev[1] * cgto_0;
    
                    unsafe {
                        let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                        self.integral_block::<T>(out_with_offset, &shls, &cgto_shape_i32, cache);
                    }
                },

                3 =>
                for idx_2 in 0..index_shape_rev[2] {
                    let shl_2 = idx_2 as i32 + shl_slices_rev[2][0];
                    let cgto_2 = cgto_locs_rel_rev[2][idx_2];
                    let shls = [shl_2, shl_1, shl_0];
                    let offset = cgto_2 + cgto_shape_rev[2] * (cgto_1 + cgto_shape_rev[1] * cgto_0);
    
                    unsafe {
                        let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                        self.integral_block::<T>(out_with_offset, &shls, &cgto_shape_i32, cache);
                    }
                },

                4 =>
                for idx_2 in 0..index_shape_rev[2] {
                    for idx_3 in 0..index_shape_rev[3] {
                        let shl_2 = idx_2 as i32 + shl_slices_rev[2][0];
                        let shl_3 = idx_3 as i32 + shl_slices_rev[3][0];
                        let cgto_2 = cgto_locs_rel_rev[2][idx_2];
                        let cgto_3 = cgto_locs_rel_rev[3][idx_3];
                        let shls = [shl_3, shl_2, shl_1, shl_0];
                        let offset = cgto_3 + cgto_shape_rev[3] * (cgto_2 + cgto_shape_rev[2] * (cgto_1 + cgto_shape_rev[1] * cgto_0));
        
                        unsafe {
                            let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                            self.integral_block::<T>(out_with_offset, &shls, &cgto_shape_i32, cache);
                        }
                    }
                },

                _ => panic!("Not known centers {n_center:}"),
            }
        });
        /* #endregion */

        /* #region 4. cleanup */
        self.optimizer_destruct();
        /* #endregion */
    }

    pub fn integral_s1<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> Vec<f64>
    where
        T: IntorBase
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices.clone(),
            None => vec![[0, self.c_nbas]; T::n_center()],
        };
        let mut out_shape = match T::n_comp() {
            1 => vec![],
            _ => vec![T::n_comp()],
        };
        let out_shape = [out_shape, self.cgto_shape::<T>(&shl_slices)].concat();
        let out_size = out_shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(out_size);
        unsafe { out.set_len(out_size) };
        self.integral_s1_inplace::<T>(&mut out, &shl_slices);
        return out;
    }
    
    pub fn integral_s2ij_inplace<T> (&mut self, out: &mut Vec<f64>, shl_slices: &Vec<[i32; 2]>)
    where
        T: IntorBase
    {
        /* #region 1. dimension definition and sanity check */

        self.check_shl_slices::<T>(shl_slices).unwrap();

        let n_comp = T::n_comp();
        let n_center = T::n_center();
        let cgto_s2ij_shape = self.cgto_shape_s2ij::<T>(shl_slices).unwrap();
        let index_shape = shl_slices.iter().map(|[shl_start, shl_stop]| (shl_stop - shl_start) as usize).collect_vec();
        let cgto_locs_rel = self.cgto_loc_slices_relative(shl_slices);

        /* #endregion */

        /* #region 2. preparation for integral engine */

        // optimizer (make integral faster)
        self.optimizer::<T>();

        // cache: thread-local
        let cache_size = self.size_of_cache::<T>(shl_slices);
        let buf_size = self.size_of_buffer::<T>(shl_slices);
        let thread_cache: Vec<Vec<f64>> = vec![vec![0.; cache_size]; rayon::current_num_threads()];
        let thread_buf: Vec<Vec<f64>> = vec![vec![0.; buf_size]; rayon::current_num_threads()];

        // out: enable mut vector by passing immut slice
        let out_const_slice = out.as_slice();

        /* #endregion */

        /* #region 3. parallel integral generation */
        
        match n_center {
            3 => {
                let out_s2ij_shape: [usize; 3] = [cgto_s2ij_shape, vec![n_comp]].concat().try_into().unwrap();

                (0..index_shape[2]).into_par_iter().for_each(|idx_k| {
                    // thread-local variables
                    let thread_index = current_thread_index().unwrap_or(0);
                    let mut cache = unsafe { cast_mut_slice(&thread_cache[thread_index]) };
                    let mut buf = unsafe { cast_mut_slice(&thread_buf[thread_index]) };
                    // output
                    let mut out = unsafe { cast_mut_slice(&out) };
                    // main engine
                    let shl_k = idx_k as i32 + shl_slices[2][0];
                    let cgto_k = cgto_locs_rel[2][idx_k];

                    for idx_j in 0..index_shape[1] {
                        for idx_i in 0..(idx_j + 1) {
                            let shl_i = idx_i as i32 + shl_slices[0][0];
                            let shl_j = idx_j as i32 + shl_slices[0][0];
                            let cgto_i = cgto_locs_rel[0][idx_i];
                            let cgto_j = cgto_locs_rel[0][idx_j];

                            let shls = [shl_i, shl_j, shl_k];
                            unsafe { self.integral_block::<T>(buf, &shls, &vec![], cache); }

                            let buf_shape = [self.cgto_size(shl_i), self.cgto_size(shl_j), self.cgto_size(shl_k), n_comp];
                            let out_offsets = [cgto_i, cgto_j, cgto_k, 0];
                            if idx_i != idx_j {
                                copy_4d_s2ij_offdiag(out, buf, &out_offsets, &buf_shape, &out_s2ij_shape);
                            } else {
                                copy_4d_s2ij_diag(out, buf, &out_offsets, &buf_shape, &out_s2ij_shape);
                            }
                        }
                    }
                })
            },
            _ => panic!(),
        }
        
        /* #endregion */

        /* #region 4. cleanup */
        self.optimizer_destruct();
        /* #endregion */
    }

    pub fn integral_s2ij<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> Vec<f64>
    where
        T: IntorBase
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices.clone(),
            None => vec![[0, self.c_nbas]; T::n_center()],
        };
        let mut out_shape = match T::n_comp() {
            1 => vec![],
            _ => vec![T::n_comp()],
        };
        let out_shape = [out_shape, self.cgto_shape_s2ij::<T>(&shl_slices).unwrap()].concat();
        let out_size = out_shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(out_size);
        unsafe { out.set_len(out_size) };
        self.integral_s2ij_inplace::<T>(&mut out, &shl_slices);
        return out;
    }
}
