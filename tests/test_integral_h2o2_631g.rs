use ndarray::prelude::*;
use ndarray::SliceInfo;
use rest_libcint::CINTR2CDATA;
use rest_libcint::cint_wrapper::*;

#[cfg(test)]
mod test_h2o2_631g {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_integral_s1_serial_int3c2e_ip1() {

        let shl_slices = vec![[2, 6], [4, 14], [0, 8]];

        let mut cint_data = initialize();

        let out = cint_data.integral_s1_serial::<int3c2e_ip1>(Some(&shl_slices));
        let out_multiplier = Array::<f64, _>::linspace(0., 10., out.len());
        let out_multiplier = out_multiplier.into_shape(out.shape()).unwrap();
        
        assert_abs_diff_eq!(out.sum(), -217.0493972185929, epsilon=1e-12);
        assert_abs_diff_eq!((out * out_multiplier).sum(), -1435.5638714741315, epsilon=1e-12);
        
        let out = cint_data.integral_s1_serial::<int3c2e_ip1>(None);
        let out_multiplier = Array::<f64, _>::linspace(0., 10., out.len());
        let out_multiplier = out_multiplier.into_shape(out.shape()).unwrap();
        
        assert_abs_diff_eq!(out.sum(), -60.66622177750967, epsilon=1e-12);
        assert_abs_diff_eq!((out * out_multiplier).sum(), 1478.8773869408533, epsilon=1e-12);
    }

    #[test]
    fn test_integral_s1_serial_inplace_int3c2e() {

        let shl_slices = vec![[7, 14], [0, 10], [3, 12]];
        let ao_ranges = vec![[11, 22], [0, 18], [3, 20]];
        let nao = 22;

        let mut cint_data = initialize();
        let mut out = Array::<f64, _>::zeros([nao, nao, nao].f());
        let out_multiplier = Array::<f64, _>::linspace(0., 10., nao * nao * nao).into_shape((nao, nao, nao)).unwrap();

        let slc_info = SliceInfo::<_, Ix3, Ix3>::try_from(
            ao_ranges.iter().map(|slc| (slc[0]..slc[1]).into()).collect::<Vec<_>>()
        ).unwrap();
        let mut out_view = out.slice_mut(slc_info);

        cint_data.optimizer::<int3c2e>();
        cint_data.integral_s1_serial_inplace::<int3c2e, _> (&mut out_view, &shl_slices);
        
        assert_abs_diff_eq!(out.sum(), 258.83910807865027, epsilon=1e-12);
        assert_abs_diff_eq!((out * out_multiplier).sum(), 1727.3637003278836, epsilon=1e-12);
    }

    #[test]
    fn test_integral_s1_serial_inplace_int3c2e_ip1() {

        let shl_slices = vec![[2, 6], [4, 14], [0, 8]];
        let ao_ranges = vec![[2, 10], [6, 22], [0, 12]];
        let nao = 22;

        let mut cint_data = initialize();
        let mut out = Array::<f64, _>::zeros([nao, nao, nao, 3].f()).into_dimensionality::<IxDyn>().unwrap();
        let out_multiplier = Array::<f64, _>::linspace(0., 10., nao * nao * nao * 3);
        let out_multiplier = out_multiplier.into_shape((nao, nao, nao, 3)).unwrap();
        let out_multiplier = out_multiplier.into_dimensionality::<IxDyn>().unwrap();

        let mut slc_info = ao_ranges.iter().map(|slc| (slc[0]..slc[1])).collect::<Vec<_>>();
        slc_info.push(0..3);
        let slc_info = SliceInfo::<_, Ix4, Ix4>::try_from(
            slc_info.iter().map(|range| range.clone().into()).collect::<Vec<_>>()
        ).unwrap();
        let mut out_view = out.slice_mut(slc_info);

        cint_data.optimizer::<int3c2e_ip1>();
        cint_data.integral_s1_serial_inplace::<int3c2e_ip1, _> (&mut out_view, &shl_slices);
        
        assert_abs_diff_eq!(out.sum(), -217.0493972185929, epsilon=1e-12);
        assert_abs_diff_eq!((out * out_multiplier).sum(), -733.3508220023117, epsilon=1e-12);
    }

    fn initialize() -> CINTR2CDATA {
        // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="6-31G").build()
        let c_atm = vec![
            vec![ 8, 20,  1, 23,  0,  0],
            vec![ 8, 24,  1, 27,  0,  0],
            vec![ 1, 28,  1, 31,  0,  0],
            vec![ 1, 32,  1, 35,  0,  0],
        ];
        let c_bas = vec![
            vec![ 0,  0,  6,  1,  0, 44, 50,  0],
            vec![ 0,  0,  3,  1,  0, 56, 59,  0],
            vec![ 0,  0,  1,  1,  0, 62, 63,  0],
            vec![ 0,  1,  3,  1,  0, 64, 67,  0],
            vec![ 0,  1,  1,  1,  0, 70, 71,  0],
            vec![ 1,  0,  6,  1,  0, 44, 50,  0],
            vec![ 1,  0,  3,  1,  0, 56, 59,  0],
            vec![ 1,  0,  1,  1,  0, 62, 63,  0],
            vec![ 1,  1,  3,  1,  0, 64, 67,  0],
            vec![ 1,  1,  1,  1,  0, 70, 71,  0],
            vec![ 2,  0,  3,  1,  0, 36, 39,  0],
            vec![ 2,  0,  1,  1,  0, 42, 43,  0],
            vec![ 3,  0,  3,  1,  0, 36, 39,  0],
            vec![ 3,  0,  1,  1,  0, 42, 43,  0],
        ];
        let c_env = vec![
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  1.8897261245650618e+00,
            1.8897261245650618e+00,  0.0000000000000000e+00,
            0.0000000000000000e+00,  3.7794522491301237e+00,
            0.0000000000000000e+00,  0.0000000000000000e+00,
            3.7794522491301237e+00,  0.0000000000000000e+00,
            3.7794522491301237e+00,  0.0000000000000000e+00,
            1.8731137000000000e+01,  2.8253936999999998e+00,
            6.4012170000000002e-01,  7.6192621981653486e-01,
            1.2923709967021120e+00,  1.4713190025449527e+00,
            1.6127780000000000e-01,  6.4297778342414225e-01,
            5.4846716999999999e+03,  8.2523495000000003e+02,
            1.8804696000000001e+02,  5.2964500000000001e+01,
            1.6897570000000002e+01,  5.7996353000000003e+00,
            2.9484245477144135e+00,  5.4265712361854543e+00,
            8.7812648932200261e+00,  1.1543212923989751e+01,
            9.9005501519376935e+00,  3.3851659942154413e+00,
            1.5539616000000001e+01,  3.5999336000000000e+00,
            1.0137617999999999e+00, -2.1905179239492316e+00,
           -9.7740552772687828e-01,  2.8862909435948061e+00,
            2.7000580000000002e-01,  9.4633487268304284e-01,
            1.5539616000000001e+01,  3.5999336000000000e+00,
            1.0137617999999999e+00,  6.3793073366236994e+00,
            4.9149097473185055e+00,  2.1579105461599575e+00,
            2.7000580000000002e-01,  5.6780702217958712e-01];
        let mut cint_data = CINTR2CDATA::new();
        cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
        return cint_data;
    }
}