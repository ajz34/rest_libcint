#[link(name = "ecp")]
extern "C" {
    fn distance_square(r1: *const f64, r2: *const f64) -> f64;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_square() {
        let r1 = [0.0, 0.0, 0.0];
        let r2 = [1.0, 0.0, 0.0];
        let result = unsafe { distance_square(r1.as_ptr(), r2.as_ptr()) };
        assert_eq!(result, 1.0);
        let r1 = [5.0, 1.0, 0.0];
        let r2 = [1.0, 0.0, 1.0];
        let result = unsafe { distance_square(r1.as_ptr(), r2.as_ptr()) };
        println!("{:?}", result);
    }
}