
use std::ops::*;

use crate::cubic::{CubicSolution, cubic_solve};
use crate::triple::Triple;
use crate::vector::Vector;

#[cfg(test)]
use crate::test::Tolerate;

pub type Matrix = Triple<Vector>;

#[inline]
const fn t<T>(x: T, y: T, z: T) -> Triple<T> { Triple{x, y, z} }

impl Matrix {
    pub const fn transpose(self: &Matrix) -> Matrix {
        let Matrix{x, y, z} = self;
        t(t(x.x, y.x, z.x), t(x.y, y.y, z.y), t(x.z, y.z, z.z))
    }

    pub fn determinant(self: &Matrix) -> f64 {
        let Matrix{x, y, z} = self;
        x.x * y.y * z.z + x.y * y.z * z.x + x.z * y.x * z.y
            - x.z * y.y * z.x - x.x * y.z * z.y - x.y * y.x * z.z
    }

    pub fn invert(self: &Matrix) -> Matrix {
        let d = self.determinant();
        let a: [[f64; 3]; 3] = (*self).into();
        let mut inv = [[0.; 3]; 3];
        for i in 0..3 {
            let i1 = if i == 2 {0} else {i + 1};
            let i2 = if i == 0 {2} else {i - 1};
            for j in 0..3 {
                let j1 = if j == 2 {0} else {j + 1};
                let j2 = if j == 0 {2} else {j - 1};
                inv[j][i] = (a[i1][j1] * a[i2][j2] - a[i1][j2] * a[i2][j1]) / d;
            }
        }
        inv.into()
    }

    pub fn eigenvalues(self: &Matrix) -> CubicSolution {
        let Matrix{x, y, z} = self;
        let a = -1.;
        let b = x.x + y.y + z.z;
        let c = z.y * y.z - y.y * z.z
            + x.z * z.x - x.x * z.z
            + x.y * y.x - x.x * y.y;
        let d = self.determinant();
        cubic_solve(a, b, c, d)
    }

    // Least-absolute eigenvalue.  Our use-case is for a covariance matrix,
    // so only postive real values are expected.  Negative values are handled
    // ad. hoc.
    pub fn least_abs_eigenvalue(self: &Matrix) -> f64 {
        let Matrix{x, y, z} = self;
        let a = -1.;
        let b = x.x + y.y + z.z;
        let c = z.y * y.z - y.y * z.z
            + x.z * z.x - x.x * z.z
            + x.y * y.x - x.x * y.y;
        let d = self.determinant();
        // If we are singular, return zero...
        if d == 0. {
            return 0.;
        }
        // Solve for 1/eigenvalues and then invert the largest.
        match cubic_solve(d, c, b, a) {
            CubicSolution::Real(_, _, w) => 1. / w,
            CubicSolution::Mixed(x, r, i) => {
                let n = r*r + i*i;
                if x*x < n { 1. / n.sqrt() } else { 1. / x }
            }
        }
    }
}

impl Mul<&Matrix> for &Matrix {
    type Output = Matrix;
    fn mul(self, r: &Matrix) -> Matrix {
        let Matrix{x, y, z} = r.transpose();
        t(t(self.x * x, self.x * y, self.x * z),
          t(self.y * x, self.y * y, self.y * z),
          t(self.z * x, self.z * y, self.z * z))
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;
    fn mul(self, r: Matrix) -> Matrix { &self * &r }
}

impl Mul<&Vector> for &Matrix {
    type Output = Vector;
    fn mul(self, &v: &Vector) -> Vector {
        Triple{x: self.x * v, y: self.y * v, z: self.z * v}
    }
}

impl Mul<Vector> for Matrix {
    type Output = Vector;
    #[inline]
    fn mul(self, v: Vector) -> Vector { &self * &v }
}

#[test]
fn test_transpose() {
    let m: Matrix = [[0,1,2], [3,4,5], [6,7,8]].into();
    let mm: [[f64; 3]; 3] = m.into();
    let nn: [[f64; 3]; 3] = m.transpose().into();
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(mm[i][j], nn[j][i]);
        }
    }
}

#[test]
fn test_eigen() {
    let m: Matrix = [[2, 0, 0], [0, 3, 0], [0, 0, 4]].into();
    assert_eq!(m.eigenvalues(), CubicSolution::Real(2., 3., 4.));
    assert_eq!(m.least_abs_eigenvalue(), 2.);

    let m: Matrix = [[1, 0, 0], [1, 2, 0], [2, 3, 3]].into();
    assert_eq!(m.eigenvalues(), CubicSolution::Real(1., 2., 3.));
    assert_eq!(m.least_abs_eigenvalue(), 1.);

    let m: Matrix = [[2, 0, 0], [0, 3, 4], [0, 4, 9]].into();
    let CubicSolution::Real(u, v, w) = m.eigenvalues() else { panic!() };
    assert_eq!(u.tolerate(1e-10), 1.);
    assert_eq!(v, 2.);
    assert_eq!(w, 11.);
    assert_eq!(m.least_abs_eigenvalue().tolerate(1e-10), 1.);

    let m: Matrix = [[0, 1, 0], [0, 0, 1], [1, 0, 0]].into();
    let CubicSolution::Mixed(u, v, w) = m.eigenvalues() else { panic!() };
    assert_eq!(u.tolerate(1e-10), 1.);
    assert_eq!(v.tolerate(1e-10), -0.5);
    assert_eq!(w.tolerate(1e-10), 3f64.sqrt() / 2.);
    assert_eq!(m.least_abs_eigenvalue(), 1.);
}

#[test]
fn test_least_abs() {
    // Singular matix, make sure we don't crash.
    let m: Matrix = [[1., 2., 3.], [3., 2., 1.], [2., 2., 2.]].into();
    assert_eq!(m.least_abs_eigenvalue(), 0.);

    // Check a couple of cases with complex eigenvalues.
    let m: Matrix = [[10., 0., 0.], [0., 1., -1.], [0., 1., 1.]].into();
    assert_eq!(m.least_abs_eigenvalue(), 2f64.sqrt());

    let m: Matrix = [[2., 0., 0.], [0., 10., -10.], [0., 10., 10.]].into();
    assert_eq!(m.least_abs_eigenvalue().tolerate(1e-10), 2.);
}

#[test]
fn test_invert() {
    let m: Matrix = [[1, 2, 3], [2, 9, 4], [4, 5, 6]].into();
    let n = m.invert();
    let i = m * n;
    let e1 = Vector::new(1.,0.,0.);
    let e2 = Vector::new(0.,1.,0.);
    let e3 = Vector::new(0.,0.,1.);
    assert!((i * e1).dist(e1) < 1e-16);
    assert!((i * e2).dist(e2) < 1e-16);
    assert!((i * e3).dist(e3) < 1e-16);
}
