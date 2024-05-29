//! # rest_libcint
//!
//! The `rest_libcint` crate provides wrappers for libcint (C).
//!
//! In order to use the crate, the libcint should be installed with the outcome library `libcint.so` stored in a reachable path by users.
//!
//! Please visit <https://github.com/sunqm/libcint> for more details about the installation and the usage of libcint
//!
//! The `CINTR2CDATA` struct groups all necessary data for using `libcint`.
//! Various kinds of analytical Gaussian-type orbital (GTO) integrals provided by `libcint` are then wrapped as the methods defined on the `CINTR2CDATA` struct.
//!
//! Currently, only four kinds of integrals are available for both spheric and Cartesian GTOs, including 
//! 1) the one-electron kinetic integral (```CINTR2CDATA::cint_ijkin```),
//! 2) the one-electron nuclear attractive integral (```CINT2CDATA::cint_ijnuc```), 
//! 3) the one-electron overlap integral (```CINTR2CDATA::cint_ijovlp```), 
//! 4) the two-electron repulsive integral (```CINTR2CDATA::cint_ijkl```).
//! The other kinds of integrals are not yet ready in the current version.
//!
//! The integrals:```Vec<f64>``` produced by the CINTR2CDATA methods aformentioned are arranged in
//! the convention of `column-major` matrices according to the definition by `libcint`. 
//!
//! # Examples
//!
//!
//! ```
//! //=============================================================================
//! // Prepare `atm`, `bas` and `env` with the same data structures of those used by `libcint`.
//! // Refer to <https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf> 
//! // for the details of these data structures.
//! //=============================================================================
//! let mut atm: Vec<Vec<i32>> = vec![];
//! atm.push(vec![2,0,0,0,0,0]);
//! atm.push(vec![4,3,0,0,0,0]);
//! let mut natm = atm.len() as i32;
//! let mut bas: Vec<Vec<i32>> = vec![];
//! bas.push(vec![0,1,1,1,1,6,7,0]);
//! bas.push(vec![0,2,1,1,1,8,9,0]);
//! let mut nbas = bas.len() as i32;
//! let mut env: Vec<f64> = vec![0.0,0.0,0.0,0.7,0.0,0.0,1.0,1.0,0.5,1.0];
//! //=============================================================================
//! // Transfer `atm`, `bas`, and `env` to the raw pointers,
//! // and organize them by the data structure of `CINTR2CDATA`.
//! //=============================================================================
//! use rust_libcint::{CINTR2CDATA,CintType};
//! let mut cint_data = CINTR2CDATA::new();
//! cint_data.initial_r2c(&atm,natm,&bas,nbas,&env);
//! //=============================================================================
//! //The 2-electron repulsive integrals (ERIs) for spheric Gaussian-type orbitals
//! //=============================================================================
//! // The GTO functions considered here are spheric.
//! // For Cartesian GTOs, replace `CintType::Spheric` by 
//! //  `CintType::Cartesian` on the following line:
//! cint_data.set_cint_type(&CintType::Spheric);
//! cint_data.cint2e_optimizer_rust();
//! let buf = cint_data.cint_ijkl(0,1,1,0);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint2e ERIs: 0.5745411555937561; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron overlap integrals for spheric Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are spheric
//! cint_data.set_cint_type(&CintType::Spheric);
//! cint_data.cint1e_ovlp_optimizer_rust();
//! let buf = cint_data.cint_ijovlp(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_ovlp: 0.7096366827378776; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron kinetic integrals for Cartesian Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are Cartesian
//! cint_data.set_cint_type(&CintType::Cartesian);
//! cint_data.cint1e_kin_optimizer_rust();
//! let buf = cint_data.cint_ijkin(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_kin : 1.5780816190296618; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron nuclear attraction integrals for Cartesian Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are Cartesian
//! cint_data.set_cint_type(&CintType::Cartesian);
//! cint_data.cint1e_nuc_optimizer_rust();
//! let buf = cint_data.cint_ijnuc(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_nuc : 4.0007622494430706; v1: {:18.16}; ",v1);
//! //=============================================================================
//! // Finally deallocate the memory by transferring the raw pointers back to RUST
//! // i.e. Vec::from_raw_parts();
//! //=============================================================================
//! cint_data.final_c2r();
//! ```

#![allow(unused)]
use core::panic;
use std::collections::btree_map::Range;
use std::process::exit;
use std::{os::raw::c_int, ptr::null, ptr::null_mut};
use std::mem::ManuallyDrop;
use itertools::Itertools;

mod cint;
pub mod cint_crafter;
use ndarray::{prelude::*, IxDynImpl, OwnedRepr, SliceArg, SliceInfoElem, SliceInfo};

use crate::cint::{CINTOpt,CINTdel_optimizer};

#[derive(Clone,Copy)]
pub enum CintType {
   Spheric,
   Cartesian,
   //Spinor,  // Not yet included
}

pub enum IJOPT {
    Ovlp,
    Kinetic,
    Nuclear,
    UNKNOWN
}

pub enum IJIPOPT {
    IPOvlp,
    IPKin,
    IPNuc,
    IPRInv
}

pub enum IP3C2E {
    IP1,
    IP2,
}

#[derive(Clone)]
pub struct CINTR2CDATA {
    c_opt: *mut CINTOpt,
    c_natm: i32,
    c_nbas: i32,
    cint_type: CintType,
    c_atm: Vec<i32>,
    c_bas: Vec<i32>,
    c_env: Vec<f64>,
}

pub mod cint_wrapper;
use crate::cint_wrapper::*;

impl CINTR2CDATA {
    /// create a new, empty CINTR2CDATA.
    pub fn new() -> CINTR2CDATA {
        CINTR2CDATA {
            c_opt: null_mut(),
            c_natm: 0,
            c_nbas: 0,
            cint_type: CintType::Spheric,
            c_atm: Vec::new(),
            c_bas: Vec::new(),
            c_env: Vec::new(),
        }
    }
    pub fn set_cint_type(&mut self, ctype: &CintType) {
        self.cint_type = *ctype;
    }
    //// 
    pub fn initial_r2c(&mut self, 
                    atm: &Vec<Vec<i32>>, natm: i32, 
                    bas: &Vec<Vec<i32>>, nbas: i32, 
                    env: &Vec<f64>) {
        self.c_atm = atm.clone().into_iter().flatten().collect::<Vec<i32>>();
        self.c_bas = bas.clone().into_iter().flatten().collect::<Vec<i32>>();
        self.c_env = env.clone();
        self.c_natm = natm;
        self.c_nbas = nbas;
        self.c_opt = null_mut();
    }

    pub fn initial_r2c_with_ecp(&mut self,
                atm: &Vec<Vec<i32>>, natm:i32,
                bas: &Vec<Vec<i32>>, nbas:i32,
                ecp: &Vec<Vec<i32>>, necp:i32,
                env: &Vec<f64>) {
        self.initial_r2c(atm, natm, bas, nbas, env);
    }

    pub fn final_c2r(&mut self) {
        self.cint_del_optimizer_rust();
    }

    pub fn cint_del_optimizer_rust(&mut self) {
        unsafe{
            CINTdel_optimizer(&mut self.c_opt);
        }
    }

    pub fn cint1e_ecp_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_2e_optimizer_rust();
        unsafe {
            self.c_opt = std::ptr::null::<CINTOpt>() as *mut CINTOpt                            
        }
    }
    pub fn cint_cgto_rust(&self, index: i32) -> i32 {
        let mut dim: i32;
        unsafe {
            dim = match self.cint_type {
                CintType::Spheric => cint::CINTcgto_spheric(index as c_int, self.c_bas.as_ptr()),
                CintType::Cartesian => cint::CINTcgto_cart(index as c_int, self.c_bas.as_ptr()),
            };
        }
        dim
    }
    pub fn cint_2c2e(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0;(di*dj) as usize];
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::int2c2e_sph(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
                CintType::Cartesian => cint::int2c2e_cart(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_ip_2c2e(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0;(di*dj*3) as usize];
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::int2c2e_ip1_sph(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
                CintType::Cartesian => cint::int2c2e_ip1_cart(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_3c2e(&mut self, i:i32,j:i32,k:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0;(di*dj*dk) as usize];
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::int3c2e_sph(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
                CintType::Cartesian => cint::int3c2e_cart(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_ijkl_by_shell(&mut self, i:i32,j:i32,k:i32,l:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut dl = self.cint_cgto_rust(l);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int,l as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = vec![0.0;(di*dj*dk*dl) as usize];
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());

        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::int2e_sph(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
                CintType::Cartesian => cint::int2e_cart(
                    c_buf, null(), c_shls,
                    self.c_atm.as_ptr(), self.c_natm,
                    self.c_bas.as_ptr(), self.c_nbas,
                    self.c_env.as_ptr(),
                    self.c_opt, null_mut()),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }

    pub fn gto_norm(n:i32,a:f64) -> f64 {
        let mut r: f64 = 0.0_f64;
        unsafe {
            r = cint::CINTgto_norm(n as c_int, a);
        }
        r
    }
    pub fn cint_ij(&mut self, i:i32,j:i32,op_name: &String) -> Vec<f64> {
        // for 1e integrals: ovlp, kinetic, and nuclear
        let op_type = if op_name.to_lowercase() ==String::from("ovlp") {
            IJOPT::Ovlp
        } else if op_name.to_lowercase() ==String::from("kinetic") {
            IJOPT::Kinetic
        } else if op_name.to_lowercase() ==String::from("nuclear") {
            IJOPT::Nuclear
        } else {
            IJOPT::UNKNOWN
            //panic!("Error:: Unknown operator for GTO-ij integrals {}", op_name)
        };
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj) as usize);
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match op_type {
                IJOPT::Ovlp => {
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_ovlp_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_ovlp_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJOPT::Kinetic => {    
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_kin_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_kin_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJOPT::Nuclear => {    
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_nuc_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_nuc_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJOPT::UNKNOWN => {1}
            };
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }

    pub fn cint_ip_ij(&mut self, i:i32,j:i32,op_name: &String) -> Vec<f64> {
        // for 1e integrals: ipovlp

        let op_type = if op_name.to_lowercase() == String::from("ipovlp") {
            IJIPOPT::IPOvlp
        } else if op_name.to_lowercase() == String::from("ipkin") {
            IJIPOPT::IPKin
        } else if op_name.to_lowercase() == String::from("ipnuc") {
            IJIPOPT::IPNuc
        } else if op_name.to_lowercase() == String::from("iprinv") {
            IJIPOPT::IPRInv
        } else {
            panic!("Error:: Unknown operator for GTO-ij-ip integrals {}", op_name)
        };
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((3*di*dj) as usize);
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match op_type {
                IJIPOPT::IPOvlp => {
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_ipovlp_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_ipovlp_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJIPOPT::IPKin => {
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_ipkin_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_ipkin_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJIPOPT::IPNuc => {
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_ipnuc_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_ipnuc_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                },
                IJIPOPT::IPRInv => {
                    match self.cint_type {
                        CintType::Spheric => cint::int1e_iprinv_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int1e_iprinv_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }
                }
            };
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
            
            //println!("i={},j={},di={},dj={}", i,j,di,dj);
            //println!("new_buf={:?}", new_buf);
            //println!("new_buf_len={}", new_buf.len());
        }
        new_buf
    }

    pub fn cint_ip_3c2e(&mut self, i:i32,j:i32,k:i32,op_name: &String) -> Vec<f64> {
        let op_type = if op_name.to_lowercase() == String::from("ip1") {
            IP3C2E::IP1
        } else if op_name.to_lowercase() == String::from("ip2") {
            IP3C2E::IP2
        } else {
            panic!("Error:: Unknown operator for GTO-3c2e-ip integrals {}", op_name)
        };

        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        println!("di,dj,dk = {} {} {}", di, dj, dk);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0_f64].repeat((3*di*dj*dk) as usize);
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match op_type { 
                IP3C2E::IP1 => {
                    match self.cint_type {
                        CintType::Spheric => cint::int3c2e_ip1_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int3c2e_ip1_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }},
                IP3C2E::IP2 => {
                    match self.cint_type {
                        CintType::Spheric => cint::int3c2e_ip2_sph(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                        CintType::Cartesian => cint::int3c2e_ip2_cart(
                            c_buf, null(), c_shls,
                            self.c_atm.as_ptr(), self.c_natm,
                            self.c_bas.as_ptr(), self.c_nbas,
                            self.c_env.as_ptr(),
                            self.c_opt, null_mut()),
                    }},
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }

}

use crate::cint_wrapper::IntorBase;

#[cfg(test)]
mod simple_test_integral {
    use ndarray::SliceInfo;

    use super::*;

    #[test]
    fn test_trait_intorbase() {
        let mut cint_data = initialize_test_h2o_631g();
        cint_data.optimizer::<int2e>();
        println!("{:?}", unsafe{*cint_data.c_opt});
        let shls_slice = vec![[0, 2], [0, 1], [1, 3], [0, 2]];
        assert_eq!(cint_data.size_of_cache::<int2e>(&shls_slice), 445);
        assert_eq!(cint_data.size_of_cache::<int2e>(&vec![]), 1341);
        assert_eq!(cint_data.cgto_size_sph(1), 1);
        assert_eq!(cint_data.cgto_size_sph(3), 3);
    }

    fn initialize_test_h2o_631g() -> CINTR2CDATA {
        // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="6-31G").build()
        let c_atm = vec![
            vec![ 8, 20,  1, 23,  0,  0],
            vec![ 1, 24,  1, 27,  0,  0],
            vec![ 1, 28,  1, 31,  0,  0]];
        let c_bas = vec![
            vec![ 0,  0,  6,  1,  0, 32, 38,  0],
            vec![ 0,  0,  3,  1,  0, 44, 47,  0],
            vec![ 0,  0,  1,  1,  0, 50, 51,  0],
            vec![ 0,  1,  3,  1,  0, 52, 55,  0],
            vec![ 0,  1,  1,  1,  0, 58, 59,  0],
            vec![ 1,  0,  3,  1,  0, 60, 63,  0],
            vec![ 1,  0,  1,  1,  0, 66, 67,  0],
            vec![ 2,  0,  3,  1,  0, 60, 63,  0],
            vec![ 2,  0,  1,  1,  0, 66, 67,  0]];
        let c_env = vec![
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            1.77634256e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
            -4.44760657e-01,  0.00000000e+00,  1.71976186e+00,  0.00000000e+00,
            5.48467170e+03,  8.25234950e+02,  1.88046960e+02,  5.29645000e+01,
            1.68975700e+01,  5.79963530e+00,  2.94842455e+00,  5.42657124e+00,
            8.78126489e+00,  1.15432129e+01,  9.90055015e+00,  3.38516599e+00,
            1.55396160e+01,  3.59993360e+00,  1.01376180e+00, -2.19051792e+00,
            -9.77405528e-01,  2.88629094e+00,  2.70005800e-01,  9.46334873e-01,
            1.55396160e+01,  3.59993360e+00,  1.01376180e+00,  6.37930734e+00,
            4.91490975e+00,  2.15791055e+00,  2.70005800e-01,  5.67807022e-01,
            1.87311370e+01,  2.82539370e+00,  6.40121700e-01,  7.61926220e-01,
            1.29237100e+00,  1.47131900e+00,  1.61277800e-01,  6.42977783e-01];
        let mut cint_data = CINTR2CDATA::new();
        cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
        cint_data
    }
}


//pub fn cint2e_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
//                   c_atm: & *mut c_int, c_natm:c_int, 
//                   c_bas: & *mut c_int, c_nbas:c_int, 
//                   c_env: & *mut f64,
//                   c_opt: & *mut CINTOpt) -> Vec<f64> {
//    //
//    buf.shrink_to_fit();
//    let mut buf = ManuallyDrop::new(buf);
//    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
//
//    shls.shrink_to_fit();
//    let mut shls = ManuallyDrop::new(shls);
//    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
//    let mut new_buf:Vec<f64>;
//    unsafe {
//        cint::cint2e_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
//        //println!("debug 1 {}", &c_buf.read());
//        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
//        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
//        //println!("debug 2");
//        //vec![0.0,0.0]
//    }
//    new_buf
//}
//pub fn cint1e_ovlp_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
//                   c_atm: & *mut c_int, c_natm:c_int, 
//                   c_bas: & *mut c_int, c_nbas:c_int, 
//                   c_env: & *mut f64,
//                   c_opt: & *mut CINTOpt) -> Vec<f64> {
//    //
//    buf.shrink_to_fit();
//    let mut buf = ManuallyDrop::new(buf);
//    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
//
//    shls.shrink_to_fit();
//    let mut shls = ManuallyDrop::new(shls);
//    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
//    let mut new_buf:Vec<f64>;
//    unsafe {
//        cint::cint1e_ovlp_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
//        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
//        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
//    }
//    new_buf
//}

#[test]
pub fn test_1() {
    //=============================================================================
    // Prepare `atm`, `bas` and `env` with the same data structures of those used by `libcint`.
    // Refer to <https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf> 
    // for the details of these data structures.
    //=============================================================================
    let mut atm: Vec<Vec<i32>> = vec![];
    let mut bas: Vec<Vec<i32>> = vec![];
    let mut env = vec![0.0;20];
    let mut ptr_env = 20;
    atm.push(vec![1,ptr_env,0,0,0,0]);
    env.extend(vec![0.0,0.0,-0.8]);
    ptr_env += 3;

    atm.push(vec![1,ptr_env,0,0,0,0]);
    env.extend(vec![0.0,0.0,0.8]);
    ptr_env += 3;

    // now for basis set 1
    // H S
    // exp  bas_1  bas_2
    // 6.0    0.7    0.4
    // 2.0    0.6    0.3
    // 0.8    0.5    0.2
    env.extend(vec![6.0,2.0,0.8]);
    [6.0,2.0,0.8].into_iter().zip([0.7,0.6,0.5]).for_each(|(exp, coeff)| {
        env.push(coeff*crate::CINTR2CDATA::gto_norm(0, exp))
    });
    [6.0,2.0,0.8].into_iter().zip([0.4,0.3,0.2]).for_each(|(exp, coeff)| {
        env.push(coeff*crate::CINTR2CDATA::gto_norm(0, exp))
    });
    bas.push(
        //   ATOM_OF, ANG_OF,NPRIM_OF, NCtr_OF, KAPPA_OF, PTR_EXP, PTR_COEFF
        vec![      0,      0,       3,       2,        0, ptr_env,ptr_env+3,0]
    );
    ptr_env += 9;
    // H P
    // exp  bas_1
    // 0.9    1.0
    env.push(0.9);
    env.push(1.0*crate::CINTR2CDATA::gto_norm(1,0.9));
    bas.push(
        //   ATOM_OF, ANG_OF,NPRIM_OF, NCtr_OF, KAPPA_OF, PTR_EXP, PTR_COEFF
        vec![      0,      1,       1,       1,        0, ptr_env,ptr_env+1,0]
    );
    ptr_env += 2;

    let mut bas_2 = bas.clone();
    bas_2[0][0] = 1;
    bas_2[1][0] = 1;
    bas.extend(bas_2);

    let mut natm = atm.len() as i32;
    let mut nbas = bas.len() as i32;
    println!("{:?}", &atm);
    println!("{:?}", &bas);
    println!("{:?}", &env);
    //let mut env: Vec<f64> = vec![0.0,0.0,0.0,0.7,0.0,0.0,1.0,1.0,0.5,1.0];
    //=============================================================================
    // Transfer `atm`, `bas`, and `env` to the raw pointers,
    // and organize them by the data structure of `CINTR2CDATA`.
    //=============================================================================
    use crate::{CINTR2CDATA,CintType};
    let mut cint_data = CINTR2CDATA::new();
    cint_data.initial_r2c(&atm,natm,&bas,nbas,&env);
    cint_data.set_cint_type(&CintType::Spheric);
    //=============================================================================
    // for int1e_nuc
    //=============================================================================
    cint_data.int1e_nuc_optimizer_rust();
    let op = String::from("nuclear");
    let buf = cint_data.cint_ij(0,1,&op);
    println!("nuc: {:?}", &buf);
    //=============================================================================
    // for cint1e_ipnuc
    //=============================================================================
    cint_data.int1e_ipnuc_optimizer_rust();
    let op = String::from("ipnuc");
    let buf = cint_data.cint_ip_ij(0,1,&op);
    println!("ipnuc: {:?}", &buf);
    //cint_data.cint2e_optimizer_rust();
    //=============================================================================
    // for cint3c2e_ip
    //=============================================================================
    cint_data.int3c2e_ip1_optimizer_rust();
    let op = String::from("ip1");
    let buf = cint_data.cint_ip_3c2e(0,1,1, &op);
    println!("3c2e_ip1: {:?}", &buf);

    cint_data.final_c2r();
}
