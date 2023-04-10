
use std::ops::*;

use crate::matrix::Matrix;
use crate::triple::Triple;

pub type Vector = Triple<f64>;

impl Vector {
    pub fn cross(self: Vector, r: Vector) -> Vector {
        Vector::new(
            self.y * r.z - self.z * r.y,
            self.z * r.x - self.x * r.z,
            self.x * r.y - self.y * r.x)
    }

    pub fn normsq(self) -> f64 { self * self }
    pub fn norm(self) -> f64 { self.normsq().sqrt() }
    pub fn unit(self) -> Vector { self / self.norm() }
    pub fn dist(self, other: Vector) -> f64 { (self - other).norm() }

    pub fn verticate(self) -> Matrix {
        let Vector{x, y, z} = self.unit();
        let d = x * x + y * y;
        let scale = if d > 1e-6 {(z.abs() - 1.0) / d} else {-0.5 - 0.125 * d};
        let a = x * x * scale + 1.0;
        let b = y * y * scale + 1.0;
        let c = x * y * scale;
        if z >= 0.0 {
            [[a, c, -x], [c, b, -y], [x, y, z]].into()
        }
        else {
            [[-a, -c, -x], [c, b, y], [x, y, z]].into()
        }
    }
}

/// Dot product.
impl Mul<Vector> for Vector {
    type Output = f64;
    fn mul(self: Vector, r: Vector) -> f64 {
        self.x * r.x + self.y * r.y + self.z * r.z
    }
}

#[cfg(test)]
fn test_verticate1(v : Vector) {
    let m = v.verticate();
    assert!((m.x.norm() - 1.0).abs() < 1e-15);
    assert!((m.y.norm() - 1.0).abs() < 1e-15);
    assert!((m.z.norm() - 1.0).abs() < 1e-15);
    assert!((m.x * m.y).abs() < 1e-15);
    assert!((m.y * m.z).abs() < 1e-15);
    assert!((m.z * m.x).abs() < 1e-15);
    assert!((m * v - Vector::new(0.0, 0.0, v.norm())).norm() < 1e-15);
}

#[test]
fn test_verticate() {
    let pos = [2., 1., 0.5, 1.1e-3, 1e-3, 1e-4, 5e-4, 1e-10, 1e-12];
    let mut values = vec![0.];
    values.extend(pos);
    values.extend(pos.iter().map(Neg::neg));
    for &x in &values {
        for &y in &values {
            for &z in &values {
                if x != 0. || y != 0. || z != 0. {
                    test_verticate1(Vector::new(x, y, z));
                }
            }
        }
    }
}

/*
impl Mul<f64> for Vector {
    type Output = Vector;
    fn mul(self: Vector, r: f64) -> Vector {
        Vector::new(self.x * r, self.y * r, self.z * r)
    }
}
 */
/*
impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self: f64, r: Vector) -> Vector {
        Vector::new(self * r.x, self * r.y, self * r.z)
    }
}
*/
