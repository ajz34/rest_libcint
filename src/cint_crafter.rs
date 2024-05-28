use std::slice::from_raw_parts_mut;

use itertools::Itertools;
use crate::cint_wrapper::*;
use crate::CINTR2CDATA;

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

impl CINTR2CDATA {

    /// Main integral engine for s1 symmetry.
    /// 
    /// Note: This function a low-level API, which is not intended to be called by user.
    pub fn integral_s1_inplace<T> (
        &mut self,
        out: &mut Vec<f64>,
        shl_slices: &Vec<[i32; 2]>,
        out_offset: usize,
        out_strides: &Vec<usize>,
    )
    where
        T: IntorBase
    {
        /* #region 1. dimension definition and sanity check */
        let n_comp = T::n_comp();
        let n_center = T::n_center();
        assert_eq!(
            shl_slices.len(), n_center,
            "length of `shl_slices` {shl_slices:?} is not the same to n_center {n_center:?}");
        let out_shape = self.shape_of_integral::<T>(shl_slices);
        let mut out_shape_without_comp = out_shape.iter().map(|&v| v as i32).collect_vec();
        if n_comp > 1 { out_shape_without_comp.pop().unwrap(); };
        
        assert_eq!(
            out_shape.len(), out_strides.len(),
            "length of `out_shape` {out_shape:?} is not the same to `out_strides` {out_strides:?}");
        /* #endregion */

        /* #region 2. cgto (atomic orbital) definition */
        let index_shape = shl_slices.iter().map(|[shl_start, shl_stop]| (shl_stop - shl_start) as usize).collect_vec();
        let index_strides = get_f_strides_from_shape(&index_shape);
        let index_size = index_shape.iter().product::<usize>();

        let cgto_locs_rel = self.cgto_loc_slices_relative(shl_slices);
        let mut cgto_ranges = vec![vec![]; n_center];
        let mut cgto_sizes = vec![vec![]; n_center];
        for (idx_center, &[shl0, shl1]) in shl_slices.iter().enumerate() {
            for idx in 0..(shl1 - shl0) {
                let idx = idx as usize;
                let cgto_loc_rel = &cgto_locs_rel[idx_center];
                cgto_ranges[idx_center].push(cgto_loc_rel[idx]..cgto_loc_rel[idx+1]);
                cgto_sizes[idx_center].push(cgto_loc_rel[idx+1] - cgto_loc_rel[idx]);
            }
        }
        /* #endregion */

        // == buffer allocation ==
        let cache_size = self.size_of_cache::<T>(shl_slices);
        let buf_size: usize = n_comp * cgto_sizes.iter().map(|cgto_size| cgto_size.iter().max().unwrap() ).product::<usize>();
        let mut cache = Vec::<f64>::with_capacity(cache_size);
        unsafe { cache.set_len(cache_size) };
        // output buffer allocation
        let mut buf = Vec::<f64>::with_capacity(buf_size);
        unsafe { buf.set_len(buf_size) };

        assert!(n_center == 3);
        assert!(n_comp == 1);

        self.optimizer::<T>();
        let mut cgto_indices = vec![0; out_shape.len()];
        
        for a0 in 0..index_shape[0] {
            let b0 = a0 as i32 + shl_slices[0][0];
            let c0 = cgto_locs_rel[0][a0];
        for a1 in 0..index_shape[1] {
            let b1 = a1 as i32 + shl_slices[1][0];
            let c1 = cgto_locs_rel[1][a1];
        for a2 in 0..index_shape[2] {
            let b2 = a2 as i32 + shl_slices[2][0];
            let c2 = cgto_locs_rel[2][a2];
            let shls = [b2, b1, b0];
            let offset = c2 + out_shape[0] * (c1 + out_shape[1] * c0);
            unsafe {
                self.integral_block::<T>(
                    out.split_at_mut(offset).1,
                    &shls, &out_shape_without_comp,
                    cache.as_mut_slice());
            }
        }}}
    }

    pub fn integral_s1<T> (&mut self, shl_slices: Option<&Vec<[i32; 2]>>) -> Vec<f64>
    where
        T: IntorBase
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices.clone(),
            None => vec![[0, self.c_nbas]; T::n_center()],
        };
        let shape = self.shape_of_integral::<T>(&shl_slices);
        let size = shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(size);
        unsafe { out.set_len(size) };
        let out_strides = get_f_strides_from_shape(&shape);
        self.integral_s1_inplace::<T>(&mut out, &shl_slices, 0, &out_strides);
        return out;
    }
    
    
}
