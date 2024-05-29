use std::slice::from_raw_parts_mut;
use itertools::Itertools;
use ndarray::parallel::prelude::IntoParallelIterator;
use ndarray::parallel::prelude::ParallelIterator;
use crate::cint_wrapper::*;
use crate::CINTR2CDATA;
use rayon;

/* #region indices computation */

#[inline(always)]
fn comput_index(offset: usize, indices: &[usize], strides: &[usize]) -> usize {
    debug_assert_eq!(indices.len(), strides.len());
    return offset + indices.iter().zip(strides).map(|(&i, &s)| i * s).sum::<usize>();
}

fn get_f_strides_from_shape(shape: &Vec<usize>) -> Vec<usize> {
    let mut strides = vec![1];
    for d in &shape[0..shape.len()-1] {
        strides.push(strides.last().unwrap() * d);
    }
    return strides;
}

fn get_c_strides_from_shape(shape: &Vec<usize>) -> Vec<usize> {
    let mut l = shape.iter().product::<usize>();
    let mut strides = vec![];
    for d in &shape[0..shape.len()] {
        l /= d;
        strides.push(l);
    }
    return strides;
}

fn get_indices_from_c_stride(index: usize, c_strides: &Vec<usize>) -> Vec<usize> {
    let mut index = index;
    c_strides.iter().map(|&d| {
        let div = index / d;
        index %= d;
        div
    }).collect_vec()
}

fn get_indices_from_f_stride(index: usize, f_strides: &Vec<usize>) -> Vec<usize> {
    let mut index = index;
    f_strides.iter().rev().map(|&d| {
        let div = index / d;
        index %= d;
        div
    }).collect_vec().into_iter().rev().collect_vec()
}

/* #endregion */

unsafe impl Send for CINTR2CDATA {}
unsafe impl Sync for CINTR2CDATA {}

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

        // iterate first two indices
        (0..(index_shape[0] * index_shape[1])).into_par_iter().for_each_with(
            {
                // thread-local: 
                let mut cache = Vec::<f64>::with_capacity(cache_size);
                unsafe { cache.set_len(cache_size) };
                cache
            },
            |cache, idx_01| {

            let idx_0 = idx_01 % index_shape[0];
            let idx_1 = idx_01 / index_shape[1];
            let shl_0 = idx_0 as i32 + shl_slices[0][0];
            let shl_1 = idx_1 as i32 + shl_slices[1][0];
            let cgto_0 = cgto_locs_rel[0][idx_0];
            let cgto_1 = cgto_locs_rel[1][idx_1];
            
            for idx_2 in 0..index_shape[2] {
                let shl_2 = idx_2 as i32 + shl_slices[2][0];
                let cgto_2 = cgto_locs_rel[2][idx_2];
                let shls = [shl_2, shl_1, shl_0];

                let offset = cgto_2 + out_shape[0] * (cgto_1 + out_shape[1] * cgto_0);

                unsafe {
                    let out_offset_ptr = out_const_slice.as_ptr() as *mut f64;
                    self.integral_block::<T>(
                        from_raw_parts_mut(out_offset_ptr.add(offset), 0),
                        &shls, &out_shape_i32,
                        cache);
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
