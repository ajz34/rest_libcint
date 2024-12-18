#[allow(non_snake_case)]
#[cfg(test)]
mod valid_ecp_sb2me4_tzvp {
    use std::time::Instant;
    use itertools::Itertools;
    use rest_libcint::prelude::*;
    use ndarray::prelude::*;
    use approx::*;

    #[test]
    fn test_ECPscalar() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let (out, _) = cint_data.integral_ecp_s1::<ECPscalar>(None);
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(
            out.sum(), 2844.388262741379, max_relative=1e-10);
        assert_relative_eq!(
            (out * scale).sum(), -2316.524453049404, max_relative=1e-10);
    }

    #[test]
    fn test_ECPscalar_ipnuc() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[7, 93], [14, 121]];
        let (out, _) = cint_data.integral_ecp_s1::<ECPscalar_ipnuc>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(
            out.sum(), -2332.269421722938, max_relative=1e-10);
        assert_relative_eq!(
            (out * scale).sum(), 702.9250770858847, max_relative=1e-10);
    }

    #[test]
    fn test_ECPscalar_ipiprinv() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[7, 93], [14, 121]];
        let (out, _) = cint_data.integral_ecp_s1::<ECPscalar_ipiprinv>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(
            out.sum(), 101.618985700688, max_relative=1e-10);
        assert_relative_eq!(
            (out * scale).sum(), 16.880187076578714, max_relative=1e-10);
    }

    #[test]
    fn test_ECPscalar_ignuc() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[7, 93], [14, 121]];
        let (out, _) = cint_data.integral_ecp_s1::<ECPscalar_ignuc>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(
            out.sum(), 45.903159687209964, max_relative=1e-10);
        assert_relative_eq!(
            (out * scale).sum(), 0.21557173801319873, max_relative=1e-10);
    }

    #[test]
    fn test_ECPscalar_ignuc_cart() {
        let mut cint_data = initialize();
        cint_data.set_cint_type(&CintType::Cartesian);
        let now = Instant::now();
        let shl_slices = vec![[7, 93], [14, 121]];
        let (out, _) = cint_data.integral_ecp_s1::<ECPscalar_ignuc>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(
            out.sum(), 111.56926252707224, max_relative=1e-10);
        assert_relative_eq!(
            (out * scale).sum(), -77.65227858432043, max_relative=1e-10);
    }

    fn initialize() -> CINTR2CDATA {
        let c_atm = vec![
            [23, 20,  4, 23,  0,  0],
            [23, 24,  4, 27,  0,  0],
            [ 6, 28,  1, 31,  0,  0],
            [ 6, 32,  1, 35,  0,  0],
            [ 6, 36,  1, 39,  0,  0],
            [ 6, 40,  1, 43,  0,  0],
            [ 1, 44,  1, 47,  0,  0],
            [ 1, 48,  1, 51,  0,  0],
            [ 1, 52,  1, 55,  0,  0],
            [ 1, 56,  1, 59,  0,  0],
            [ 1, 60,  1, 63,  0,  0],
            [ 1, 64,  1, 67,  0,  0],
            [ 1, 68,  1, 71,  0,  0],
            [ 1, 72,  1, 75,  0,  0],
            [ 1, 76,  1, 79,  0,  0],
            [ 1, 80,  1, 83,  0,  0],
            [ 1, 84,  1, 87,  0,  0],
            [ 1, 88,  1, 91,  0,  0],
        ];
        let c_bas = vec![
            [  0,   0,   4,   1,   0, 144, 148,   0],
            [  0,   0,   2,   1,   0, 152, 154,   0],
            [  0,   0,   1,   1,   0, 156, 157,   0],
            [  0,   0,   1,   1,   0, 158, 159,   0],
            [  0,   0,   1,   1,   0, 160, 161,   0],
            [  0,   0,   1,   1,   0, 162, 163,   0],
            [  0,   1,   3,   1,   0, 164, 167,   0],
            [  0,   1,   3,   1,   0, 170, 173,   0],
            [  0,   1,   1,   1,   0, 176, 177,   0],
            [  0,   1,   1,   1,   0, 178, 179,   0],
            [  0,   1,   1,   1,   0, 180, 181,   0],
            [  0,   2,   6,   1,   0, 182, 188,   0],
            [  0,   2,   1,   1,   0, 194, 195,   0],
            [  0,   2,   1,   1,   0, 196, 197,   0],
            [  0,   3,   1,   1,   0, 198, 199,   0],
            [  0,   3,   1,   1,   0, 200, 201,   0],
            [  1,   0,   4,   1,   0, 144, 148,   0],
            [  1,   0,   2,   1,   0, 152, 154,   0],
            [  1,   0,   1,   1,   0, 156, 157,   0],
            [  1,   0,   1,   1,   0, 158, 159,   0],
            [  1,   0,   1,   1,   0, 160, 161,   0],
            [  1,   0,   1,   1,   0, 162, 163,   0],
            [  1,   1,   3,   1,   0, 164, 167,   0],
            [  1,   1,   3,   1,   0, 170, 173,   0],
            [  1,   1,   1,   1,   0, 176, 177,   0],
            [  1,   1,   1,   1,   0, 178, 179,   0],
            [  1,   1,   1,   1,   0, 180, 181,   0],
            [  1,   2,   6,   1,   0, 182, 188,   0],
            [  1,   2,   1,   1,   0, 194, 195,   0],
            [  1,   2,   1,   1,   0, 196, 197,   0],
            [  1,   3,   1,   1,   0, 198, 199,   0],
            [  1,   3,   1,   1,   0, 200, 201,   0],
            [  2,   0,   6,   1,   0, 104, 110,   0],
            [  2,   0,   2,   1,   0, 116, 118,   0],
            [  2,   0,   1,   1,   0, 120, 121,   0],
            [  2,   0,   1,   1,   0, 122, 123,   0],
            [  2,   0,   1,   1,   0, 124, 125,   0],
            [  2,   1,   4,   1,   0, 126, 130,   0],
            [  2,   1,   1,   1,   0, 134, 135,   0],
            [  2,   1,   1,   1,   0, 136, 137,   0],
            [  2,   2,   1,   1,   0, 138, 139,   0],
            [  2,   2,   1,   1,   0, 140, 141,   0],
            [  2,   3,   1,   1,   0, 142, 143,   0],
            [  3,   0,   6,   1,   0, 104, 110,   0],
            [  3,   0,   2,   1,   0, 116, 118,   0],
            [  3,   0,   1,   1,   0, 120, 121,   0],
            [  3,   0,   1,   1,   0, 122, 123,   0],
            [  3,   0,   1,   1,   0, 124, 125,   0],
            [  3,   1,   4,   1,   0, 126, 130,   0],
            [  3,   1,   1,   1,   0, 134, 135,   0],
            [  3,   1,   1,   1,   0, 136, 137,   0],
            [  3,   2,   1,   1,   0, 138, 139,   0],
            [  3,   2,   1,   1,   0, 140, 141,   0],
            [  3,   3,   1,   1,   0, 142, 143,   0],
            [  4,   0,   6,   1,   0, 104, 110,   0],
            [  4,   0,   2,   1,   0, 116, 118,   0],
            [  4,   0,   1,   1,   0, 120, 121,   0],
            [  4,   0,   1,   1,   0, 122, 123,   0],
            [  4,   0,   1,   1,   0, 124, 125,   0],
            [  4,   1,   4,   1,   0, 126, 130,   0],
            [  4,   1,   1,   1,   0, 134, 135,   0],
            [  4,   1,   1,   1,   0, 136, 137,   0],
            [  4,   2,   1,   1,   0, 138, 139,   0],
            [  4,   2,   1,   1,   0, 140, 141,   0],
            [  4,   3,   1,   1,   0, 142, 143,   0],
            [  5,   0,   6,   1,   0, 104, 110,   0],
            [  5,   0,   2,   1,   0, 116, 118,   0],
            [  5,   0,   1,   1,   0, 120, 121,   0],
            [  5,   0,   1,   1,   0, 122, 123,   0],
            [  5,   0,   1,   1,   0, 124, 125,   0],
            [  5,   1,   4,   1,   0, 126, 130,   0],
            [  5,   1,   1,   1,   0, 134, 135,   0],
            [  5,   1,   1,   1,   0, 136, 137,   0],
            [  5,   2,   1,   1,   0, 138, 139,   0],
            [  5,   2,   1,   1,   0, 140, 141,   0],
            [  5,   3,   1,   1,   0, 142, 143,   0],
            [  6,   0,   3,   1,   0,  92,  95,   0],
            [  6,   0,   1,   1,   0,  98,  99,   0],
            [  6,   0,   1,   1,   0, 100, 101,   0],
            [  6,   1,   1,   1,   0, 102, 103,   0],
            [  7,   0,   3,   1,   0,  92,  95,   0],
            [  7,   0,   1,   1,   0,  98,  99,   0],
            [  7,   0,   1,   1,   0, 100, 101,   0],
            [  7,   1,   1,   1,   0, 102, 103,   0],
            [  8,   0,   3,   1,   0,  92,  95,   0],
            [  8,   0,   1,   1,   0,  98,  99,   0],
            [  8,   0,   1,   1,   0, 100, 101,   0],
            [  8,   1,   1,   1,   0, 102, 103,   0],
            [  9,   0,   3,   1,   0,  92,  95,   0],
            [  9,   0,   1,   1,   0,  98,  99,   0],
            [  9,   0,   1,   1,   0, 100, 101,   0],
            [  9,   1,   1,   1,   0, 102, 103,   0],
            [ 10,   0,   3,   1,   0,  92,  95,   0],
            [ 10,   0,   1,   1,   0,  98,  99,   0],
            [ 10,   0,   1,   1,   0, 100, 101,   0],
            [ 10,   1,   1,   1,   0, 102, 103,   0],
            [ 11,   0,   3,   1,   0,  92,  95,   0],
            [ 11,   0,   1,   1,   0,  98,  99,   0],
            [ 11,   0,   1,   1,   0, 100, 101,   0],
            [ 11,   1,   1,   1,   0, 102, 103,   0],
            [ 12,   0,   3,   1,   0,  92,  95,   0],
            [ 12,   0,   1,   1,   0,  98,  99,   0],
            [ 12,   0,   1,   1,   0, 100, 101,   0],
            [ 12,   1,   1,   1,   0, 102, 103,   0],
            [ 13,   0,   3,   1,   0,  92,  95,   0],
            [ 13,   0,   1,   1,   0,  98,  99,   0],
            [ 13,   0,   1,   1,   0, 100, 101,   0],
            [ 13,   1,   1,   1,   0, 102, 103,   0],
            [ 14,   0,   3,   1,   0,  92,  95,   0],
            [ 14,   0,   1,   1,   0,  98,  99,   0],
            [ 14,   0,   1,   1,   0, 100, 101,   0],
            [ 14,   1,   1,   1,   0, 102, 103,   0],
            [ 15,   0,   3,   1,   0,  92,  95,   0],
            [ 15,   0,   1,   1,   0,  98,  99,   0],
            [ 15,   0,   1,   1,   0, 100, 101,   0],
            [ 15,   1,   1,   1,   0, 102, 103,   0],
            [ 16,   0,   3,   1,   0,  92,  95,   0],
            [ 16,   0,   1,   1,   0,  98,  99,   0],
            [ 16,   0,   1,   1,   0, 100, 101,   0],
            [ 16,   1,   1,   1,   0, 102, 103,   0],
            [ 17,   0,   3,   1,   0,  92,  95,   0],
            [ 17,   0,   1,   1,   0,  98,  99,   0],
            [ 17,   0,   1,   1,   0, 100, 101,   0],
            [ 17,   1,   1,   1,   0, 102, 103,   0],
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
            1.2400000000000000e+02,  8.0000000000000000e+00,
           -2.5310584098499369e+00,  8.4277726023886190e-01,
           -2.4052374398118572e+00,  0.0000000000000000e+00,
            2.5310584098499369e+00, -8.4277726023886190e-01,
           -2.4052374398118572e+00,  0.0000000000000000e+00,
           -2.6537334016303631e+00,  2.0870488887454450e+00,
            1.5773204755905210e+00,  0.0000000000000000e+00,
           -4.0857793105660818e+00, -2.9504747139158982e+00,
           -1.6009297689277748e+00,  0.0000000000000000e+00,
            4.0857793105660818e+00,  2.9504747139158982e+00,
           -1.6009297689277748e+00,  0.0000000000000000e+00,
            2.6537334016303631e+00, -2.0870488887454450e+00,
            1.5773204755905210e+00,  0.0000000000000000e+00,
           -1.3212707871233360e+00,  3.6280404189405715e+00,
            1.9062048765219717e+00,  0.0000000000000000e+00,
           -2.1941889144873530e+00,  5.4859913467416477e-01,
            2.8699842862840113e+00,  0.0000000000000000e+00,
           -4.5376960126935924e+00,  2.7823488797641889e+00,
            2.0506711317242523e+00,  0.0000000000000000e+00,
           -3.8172854354179337e+00, -4.2123750116234504e+00,
           -3.2104045959319141e+00,  0.0000000000000000e+00,
           -3.1946252316137058e+00, -3.8099321520192873e+00,
            4.8712944299290760e-02,  0.0000000000000000e+00,
           -6.1123352881667969e+00, -2.8249478230573231e+00,
           -1.2363219286457627e+00,  0.0000000000000000e+00,
            3.8172854354179337e+00,  4.2123750116234504e+00,
           -3.2104045959319141e+00,  0.0000000000000000e+00,
            6.1123352881667969e+00,  2.8249478230573231e+00,
           -1.2363219286457627e+00,  0.0000000000000000e+00,
            3.1946252316137058e+00,  3.8099321520192873e+00,
            4.8712944299290760e-02,  0.0000000000000000e+00,
            1.3212707871233360e+00, -3.6280404189405715e+00,
            1.9062048765219717e+00,  0.0000000000000000e+00,
            4.5376960126935924e+00, -2.7823488797641889e+00,
            2.0506711317242523e+00,  0.0000000000000000e+00,
            2.1941889144873530e+00, -5.4859913467416477e-01,
            2.8699842862840113e+00,  0.0000000000000000e+00,
            3.4061340999999999e+01,  5.1235746000000004e+00,
            1.1646626000000000e+00,  9.0618446120248586e-01,
            1.6354784928239057e+00,  2.4145128304249659e+00,
            3.2723041000000003e-01,  1.0930883523645869e+00,
            1.0307241000000000e-01,  4.5959135109675275e-01,
            8.0000000000000004e-01,  2.2072263710762661e+00,
            1.3575349682000000e+04,  2.0352333679999999e+03,
            4.6322562359000000e+02,  1.3120019597999999e+02,
            4.2853015890999998e+01,  1.5584185765999999e+01,
            1.9269786369447457e+00,  3.5965017404150745e+00,
            6.1382671840596235e+00,  9.5394800177297636e+00,
            1.2777459273222963e+01,  1.3125218210409033e+01,
            6.2067138507999999e+00,  2.5764896526999999e+00,
            6.5167955479977069e+00,  1.9311152716735143e+00,
            5.7696339418999998e-01,  1.6725388683875464e+00,
            2.2972831358000001e-01,  8.3835105961496525e-01,
            9.5164440027999994e-02,  4.3288351417905013e-01,
            3.4697232243999999e+01,  7.9582622825999998e+00,
            2.3780826883000001e+00,  8.1433208183000005e-01,
            2.7827633154407163e+00,  2.9702240843169960e+00,
            2.6011391897611520e+00,  1.6425869687454335e+00,
            2.8887547253000001e-01,  6.1783531504681743e-01,
            1.0056823670999999e-01,  1.6521917035633993e-01,
            1.0970000000000000e+00,  3.0682517054988114e+00,
            3.1800000000000000e-01,  3.5137984426073676e-01,
            7.6100000000000001e-01,  1.0669052186285761e+00,
            1.6124199933000000e+03,  2.3884452096999999e+02,
            2.3998118809000001e+01,  1.5193124213000001e+01,
            4.7398852731266844e-01,  5.3111555217407824e-01,
           -3.4950640089470562e+00,  2.1794489183304044e+01,
            1.1736409733000000e+01,  6.5259774793999998e+00,
            8.7478174468651169e+00,  4.8447468436724117e+00,
            2.0247451872000002e+00,  4.2883754956682498e+00,
            9.7113418587000000e-01,  2.4715787696517073e+00,
            2.4254333997999999e-01,  8.7318655419367219e-01,
            9.2206608176999993e-02,  4.2275286483061619e-01,
            2.1568393354000000e+02,  1.6374479088000001e+01,
            9.7216283345000001e+00,  3.0415339314391066e+00,
            3.4302369393967169e+01, -6.6024353304020806e+01,
            2.7982643154000000e+00,  1.4711045033000001e+00,
            7.5165385300999998e-01,  4.5057661963693887e+00,
            2.1859318084540815e+00,  3.5073463932127641e-01,
            3.3168699849000000e-01,  7.3433622558653044e-01,
            1.3931606365999999e-01,  2.4830524102829885e-01,
            5.6399307526000003e-02,  8.0181948753247551e-02,
            1.1590312253000000e+02,  3.0474233720000001e+01,
            1.8228418239000000e+01,  4.3291456646000004e+00,
            2.1294818495999999e+00,  9.9682636692000004e-01,
            6.2185747276141159e+00,  6.7119020458836074e+00,
           -4.8553708683987953e+00,  7.5565698755090089e+00,
            4.5866372169759106e+00,  1.0954341006038497e+00,
            4.3347239862999998e-01,  6.0424340880994176e-01,
            1.4000000000000001e-01,  8.3609080734344948e-02,
            3.2369705999999998e-01,  1.5589154464098992e-01,
            1.1000000000000001e+00,  2.4442403165721562e+00,
            1.4449294999999999e+01,  1.4444978000000001e+01,
           -2.0296137999999999e+01, -1.5366801000000001e+01,
            1.6330864999999999e+01,  1.4449294999999999e+01,
            1.4444978000000001e+01,  8.5565420000000003e+00,
            2.8107158099999998e+02,  2.0296137999999999e+01,
            1.5366801000000001e+01,  6.1716603999999997e+01,
            1.4470337000000001e+01,  1.4449294999999999e+01,
            1.4444978000000001e+01,  1.3816193999999999e+01,
            8.4249240000000007e+00,  8.0927279999999993e+00,
            6.7457380000000001e+01,  2.0296137999999999e+01,
            1.5366801000000001e+01,  1.3493350300000000e+02,
            1.4716343999999999e+01,  2.9518512000000001e+01,
            1.5146319000000000e+01,  1.4886331000000000e+01,
            1.4449294999999999e+01,  1.4444978000000001e+01,
            5.9082670000000004e+00,  5.5943220000000000e+00,
            5.3143465999999997e+01,  3.5447814999999999e+01,
            2.0296137999999999e+01,  1.5366801000000001e+01,
            9.1792230000000004e+00,  1.3240252999999999e+01,
        ];
        let c_ecp = vec![
            [  0,  -1,   2,   2,   0, 202, 204,   0],
            [  0,   0,   4,   2,   0, 206, 210,   0],
            [  0,   1,   6,   2,   0, 214, 220,   0],
            [  0,   2,   6,   2,   0, 226, 232,   0],
            [  1,  -1,   2,   2,   0, 202, 204,   0],
            [  1,   0,   4,   2,   0, 206, 210,   0],
            [  1,   1,   6,   2,   0, 214, 220,   0],
            [  1,   2,   6,   2,   0, 226, 232,   0],
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