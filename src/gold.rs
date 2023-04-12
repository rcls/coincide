
use crate::matrix::Matrix;
use crate::triple::Triple;
use crate::vector::Vector;

pub const GOLD: f64 = 1.618033988749895;

const fn t<T>(x: T, y: T, z: T) -> Triple<T> { Triple::new(x, y, z) }

pub const ID: Matrix = t(t(1., 0., 0.), t(0., 1., 0.), t(0., 0., 1.));

/// Clockwise around [1,0,gold].
pub const ROT5: Matrix = t(
    t(0.5, 0.5*GOLD, 0.5*GOLD-0.5),
    t(-0.5*GOLD, 0.5*GOLD-0.5, 0.5),
    t(0.5*GOLD-0.5, -0.5, 0.5*GOLD));

/// Clockwise around [1,1,1].
pub const ROT3: Matrix = t(t(0., 1., 0.), t(0., 0., 1.), t(1., 0., 0.));

pub struct IcoGroup {
    /// The 120 members of the ico. symmetry group, with the identity at
    /// index 0.
    pub matrices: [Matrix; 120],
    /// Map a matrix to an index.
    pub mat_to_index: std::collections::BTreeMap<Matrix, u8>,
    /// Matrix multiplication as a 2-d array of indexes.
    pub mult_index: [[u8; 120]; 120],
    /// Matrix inverse as an array of indexes.
    pub inv_index: [u8; 120],
    /// A basic triangle, the 120 images of this tile the sphere.
    pub corners: Matrix,
    /// Weighted sum of corners to give direction of [1,1,4*GOLD].
    pub trunc_ico_dodec: Vector,
}

impl IcoGroup {
    pub fn new() -> IcoGroup {
        let mut r = IcoGroup{
            matrices: ico_group(), mat_to_index: Default::default(),
            mult_index: [[0; 120]; 120], inv_index: [0; 120],
            corners: Matrix::default(), trunc_ico_dodec: Vector::default()
        };
        for (index, &matrix) in r.matrices.iter().enumerate() {
            r.mat_to_index.insert(matrix, index as u8);
        }
        for (li, lm) in r.mult_index.iter_mut().zip(r.matrices.iter()) {
            for (ri, rm) in li.iter_mut().zip(r.matrices.iter()) {
                *ri = r.mat_to_index[&(lm * rm).normalize()];
            }
        }
        for (i, m) in r.inv_index.iter_mut().zip(r.matrices.iter()) {
            *i = r.mat_to_index[&m.transpose()];
        }

        // Norm 1
        r.corners.x = Vector::new(0., 0., 1.).unit();
        // Norm² = 1 + GOLD² = GOLD + 1.
        r.corners.y = Vector::new(1., 0., GOLD).unit();
        // Norm² = (GOLD-1)² + GOLD² = 2GOLD² - 2GOLD + 1 = 3
        r.corners.z = Vector::new(0., GOLD - 1., GOLD).unit();

        let tid = r.corners.transpose().invert() *
            Vector::new(1., 1., 4. * GOLD);
        let tid_sum = tid.x + tid.y + tid.z;
        r.trunc_ico_dodec = tid / tid_sum;
        r
    }

    pub fn map(&self, i: u8) -> &Matrix { &self.matrices[i as usize] }
}

pub fn ico_group() -> [Matrix; 120] {
    let result = &mut [ID; 120];
    fn step(r: &mut [Matrix; 120], m: Matrix,
            start: usize, multiply: usize) -> usize {
        for i in start..start * multiply {
            r[i] = (r[i-start] * m).normalize();
        }
        start * multiply
    }
    let s = 1;
    let s = step(result, ROT5, s, 5);
    let s = step(result, ROT3, s, 3);
    let rx = t(t(1., 0., 0.), t(0., -1., 0.), t(0., 0., -1.));
    let s = step(result, rx, s, 2);
    let s = step(result, ROT3 * rx * ROT3 * ROT3, s, 2);
    let s = step(result, -ID, s, 2);
    // let s = step(result, ROT3 * ROT3 * rx * ROT3, s, 2);
    assert_eq!(s, 120);
    *result
}


pub fn integerise(x: f64) -> Option<(i32, i32)> {
    let r = x.round();
    let f = x - r;
    for i in 0..=10 {
        let g = i as f64 * GOLD;
        let fg = g - g.round();
        if (f - fg).abs() < 1e-7 {
            // f = i·ϕ - round(i·ϕ)
            // x = r+f = r - round(i·ϕ) + i·ϕ
            return Some((r as i32 - g.round() as i32, i));
        }
        if (f + fg).abs() < 1e-7 {
            return Some((r as i32 + g.round() as i32, -i));
        }
    }
    None
}

pub fn rationalise(x: f64) -> (f64, f64) {
    fn rat0(x: f64, up: f64, down: f64) -> Option<(f64, f64)> {
        let xx = x * up;
        let r = xx.round();
        let f = xx - r;
        let (i, q) = integerise(f)?;
        Some(((i as f64 + r) * down, q as f64 * down))
    }
    // Deal to integers and half integers.
    let f = x - x.round();
    if f.abs() < 1e-7 {
        (x.round(), 0.0)
    }
    else {
        rat0(x, 1.0, 1.0).unwrap_or_else(
            || rat0(x, 2.0, 0.5).unwrap_or_else(
                || rat0(x, 5.0, 0.2).unwrap_or(
                (x, 0.0))))
    }
}

pub trait Normalize {
    fn normalize(&self) -> Self;
}

impl Normalize for f64 {
    fn normalize(&self) -> f64 {
        let (i, q) = rationalise(*self);
        i + GOLD * q
    }
}

impl<T: Normalize> Normalize for Triple<T> {
    fn normalize(&self) -> Self {
        Triple::new(self.x.normalize(), self.y.normalize(), self.z.normalize())
    }
}

#[test]
fn check_gold_const()
{
    assert_eq!(GOLD, (5f64.sqrt() + 1.) / 2.);
}

#[test]
fn check_norm() {
    assert_eq!((1.0 + GOLD + 1e-10).normalize(), 1.0 + GOLD);
    assert_eq!(2.50000001.normalize(), 2.5.normalize());
    assert_eq!(rationalise(std::f64::consts::PI.normalize()),
               (std::f64::consts::PI, 0.));
    assert_eq!((0.2000000001 * (GOLD + 1.)).normalize(), 0.2 * GOLD + 0.2);
}

#[test]
fn test_rot5() {
    let r2 = &ROT5 * &ROT5;
    let r4 = &r2 * &r2;
    let i = &r4 * &ROT5;
    let z = i - ID;
    for v in &z {
        for x in v {
            assert!(x.abs() < 1e-10);
        }
    }
}

#[test]
fn test_corners() {
    let g = IcoGroup::new();
    assert!(g.corners.y.dist(ROT5 * g.corners.y) < 1e-7);
    assert!((ROT5 * g.corners.z).dist(Vector::new(1., 1., 1.).unit()) < 1e-7);
    assert!((ROT5 * g.corners.z).cross([1., 1., 1.].into()).norm() < 1e-7);

    assert_eq!(g.map(0) * &Vector::new(1., 2., 3.), [1., 2., 3.].into());

    assert!(g.trunc_ico_dodec.x > 0.);
    assert!(g.trunc_ico_dodec.y > 0.);
    assert!(g.trunc_ico_dodec.z > 0.);
    assert_eq!(g.trunc_ico_dodec * Vector::new(1., 1., 1.), 1.);

    assert!((g.corners.transpose() * g.trunc_ico_dodec).unit().dist(
        Vector::new(1., 1., 4. * GOLD).unit()) < 2e-16);
}
