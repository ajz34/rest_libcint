use std::slice::from_raw_parts_mut;
use std::sync::mpsc::channel;
use itertools::Itertools;
use rayon::current_thread_index;
use rayon::max_num_threads;
use crate::cint_wrapper::*;
use crate::CINTR2CDATA;
use rayon::prelude::*;

unsafe impl Send for CINTR2CDATA {}
unsafe impl Sync for CINTR2CDATA {}

unsafe fn cast_mut_slice<T> (slc: &[T]) -> &mut [T] {
    let len = slc.len();
    let ptr = slc.as_ptr() as *mut T;
    let mut mut_vector = from_raw_parts_mut(ptr, len);
    return mut_vector;
}

impl CINTR2CDATA {

    /// Main integral engine for s1 symmetry.
    /// 
    /// This function a low-level API, which is not intended to be called by user.
    /// This function only works for f-contiguous integral (PySCF convention).
    pub fn integral_s1_inplace<T> (
        &mut self,
        out: &mut Vec<f64>,
        shl_slices: &Vec<[i32; 2]>
    )
    where
        T: IntorBase
    {
        /* #region 1. dimension definition and sanity check */
        let n_comp = T::n_comp();
        let n_center = T::n_center();
        assert!(
            shl_slices.len() == n_center,
            "length of `shl_slices` {shl_slices:?} is not the same to n_center {n_center:?}");
        let out_shape = self.shape_of_integral::<T>(shl_slices);
        assert!(
            out_shape.iter().product::<usize>() >= out.len(),
            "length of `out` seems not enough");
        let out_shape_i32 = out_shape.iter().map(|&v| v as i32).collect::<Vec<i32>>();
        /* #endregion */

        /* #region 2. cgto (atomic orbital) definition */
        let index_shape = shl_slices.iter().map(|[shl_start, shl_stop]| (shl_stop - shl_start) as usize).collect_vec();
        let cgto_locs_rel = self.cgto_loc_slices_relative(shl_slices);
        /* #endregion */

        /* #region 3. final preparation for integral engine */
        let cache_size = self.size_of_cache::<T>(shl_slices);
        self.optimizer::<T>();
        let mut cgto_indices = vec![0; out_shape.len()];
        let out_const_slice = out.as_slice();
        /* #endregion */

        let thread_cache: Vec<Vec<f64>> = vec![vec![0.; cache_size]; rayon::current_num_threads()];

        // iterate first two indices
        (0..(index_shape[0] * index_shape[1])).into_par_iter().for_each(
            |idx_01| {

            let idx_0 = idx_01 % index_shape[0];
            let idx_1 = idx_01 / index_shape[1];
            let shl_0 = idx_0 as i32 + shl_slices[0][0];
            let shl_1 = idx_1 as i32 + shl_slices[1][0];
            let cgto_0 = cgto_locs_rel[0][idx_0];
            let cgto_1 = cgto_locs_rel[1][idx_1];
            
            let thread_index = current_thread_index().unwrap();
            let mut cache = unsafe { cast_mut_slice(&thread_cache[thread_index]) };
            
            for idx_2 in 0..index_shape[2] {

                let shl_2 = idx_2 as i32 + shl_slices[2][0];
                let cgto_2 = cgto_locs_rel[2][idx_2];
                let shls = [shl_2, shl_1, shl_0];

                let offset = cgto_2 + out_shape[0] * (cgto_1 + out_shape[1] * cgto_0);

                unsafe {
                    let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                    self.integral_block::<T>(out_with_offset, &shls, &out_shape_i32, cache);
                }
            }
        })
    }

    pub fn integral_s1<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> Vec<f64>
    where
        T: IntorBase
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices.clone(),
            None => vec![[0, self.c_nbas]; T::n_center()],
        };
        let out_shape = self.shape_of_integral::<T>(&shl_slices);
        let out_size = out_shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(out_size);
        unsafe { out.set_len(out_size) };
        self.integral_s1_inplace::<T>(&mut out, &shl_slices);
        return out;
    }
    
    
}
