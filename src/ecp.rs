#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct ECPOpt {
    pub u_ecp: *mut f64,
}

#[link(name = "ecp")]
extern "C" {
    pub fn ECPscalar_sph(
        out: *mut f64,
        dims: *const ::std::os::raw::c_int,
        shls: *const ::std::os::raw::c_int,
        atm: *const ::std::os::raw::c_int,
        natm: ::std::os::raw::c_int,
        bas: *const ::std::os::raw::c_int,
        nbas: ::std::os::raw::c_int,
        env: *const f64,
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> ::std::os::raw::c_int;
    pub fn ECPscalar_cart(
        out: *mut f64,
        dims: *const ::std::os::raw::c_int,
        shls: *const ::std::os::raw::c_int,
        atm: *const ::std::os::raw::c_int,
        natm: ::std::os::raw::c_int,
        bas: *const ::std::os::raw::c_int,
        nbas: ::std::os::raw::c_int,
        env: *const f64,
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> ::std::os::raw::c_int;
}