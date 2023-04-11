
use crate::matrix::Matrix;
use crate::vector::Vector;

#[cfg(test)]
use crate::test::Tolerate;

pub fn jacobi(m: &Matrix) -> (Matrix, Vector) {
    // m should be symmetric.  Calculate:
    // j * m * j.transpose()

    let mut diag = *m;
    let mut rot: Matrix = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into();

    for _i in 0..9 {
        // Stay symmetric!  Shouldn't be necessary, but just in case...
        diag.y.x = diag.x.y;
        diag.z.x = diag.x.z;
        diag.z.y = diag.y.z;

        let axy = diag.x.y.abs();
        let axz = diag.x.z.abs();
        let ayz = diag.y.z.abs();

        if axy > axz && axy > ayz {
            let (c, s) = jacobi2(diag.x.x, diag.y.y, diag.x.y);
            let r: Matrix = [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]].into();
            diag = r * diag * r.transpose();
            rot = r * rot;
        }
        else if !(axy > axz) && axz > ayz {
            let (c, s) = jacobi2(diag.x.x, diag.z.z, diag.x.z);
            let r: Matrix = [[c, 0., s], [0., 1., 0.], [-s, 0., c]].into();
            diag = r * diag * r.transpose();
            rot = r * rot;
        }
        else {
            let (c, s) = jacobi2(diag.y.y, diag.z.z, diag.y.z);
            let r: Matrix = [[1., 0., 0.], [0., c, s], [0., -s, c]].into();
            diag = r * diag * r.transpose();
            rot = r * rot;
        }
    }
    (rot, Vector::new(diag.x.x, diag.y.y, diag.z.z))
}

fn jacobi2(a: f64, d: f64, b: f64) -> (f64, f64) {
    // Sym. matrix [a,b;b,d], find c,s = cos θ, sin θ such that
    //
    // [c,s;-s,c] * [a,b;b,d] * [c,-s;s,c]
    //
    // is diagonal.

    // 2·c·s·(a-d)/2 + (c·c - s·s) b = 0.
    // 2·c·s = sin 2θ, c·c - s·s = (2 c·c - 1) = cos 2θ

    let h = (a - d).hypot(2. * b);
    if h == 0. {
        // It doesn't matter what we do :-)
        return (1., 0.);
    }
    let cos2 = (a - d) / h;
    let sin2by2 = b / h;
    // c² = 1/2 + cos2/2
    // s = c·s / c = sin2by2 / c
    //
    // An overall sign flip of both c,s doesn't matter so we take c to be
    // positive.
    let c = (0.5 + 0.5 * cos2).sqrt();
    // If the test fails then we have s = ±1 and zeros or rounding errors.
    // The sign of s doesn't matter in that case.
    if sin2by2.abs() < c { (c, sin2by2 / c) } else { (1., 0.) }
}

#[cfg(test)]
fn test_jacobi2_single(a: f64, d: f64, b: f64) {
    let (c, s) = jacobi2(a, d, b);

    let aa: Matrix = [[a, b, 0.], [b, d, 0.], [0., 0., 1.]].into();
    let cs: Matrix = [[c, s, 0.], [-s, c, 0.], [0., 0., 1.]].into();

    let product = cs * aa * cs.transpose();
    assert!(product.x.y.abs() < 1e-10);
    assert!(product.y.x.abs() < 1e-10);

    assert!((product.x.x * cs.x).dist(aa * cs.x).abs() < 1e-10);
    assert!((product.y.y * cs.y).dist(aa * cs.y).abs() < 1e-10);
}

#[test]
fn test_jacobi2() {
    test_jacobi2_single(1., 1., 2.);
    test_jacobi2_single(1., 2., 0.);
    test_jacobi2_single(1., 1., -2.);
    test_jacobi2_single(1., -2., 3.);
    test_jacobi2_single(3.0, -0.4082482904638631, 1.4993054670864636);
}

#[cfg(test)]
fn test_jacobi_single(m: &Matrix) {
    let (rot, eigen) = jacobi(&m);
    let diag = &rot * m * rot.transpose();
    println!("{:?}", eigen);
    const EPSILON: f64 = 1e-14;
    assert!(diag.x.y.abs() < EPSILON);
    assert!(diag.y.x.abs() < EPSILON);
    assert!(diag.x.z.abs() < EPSILON);
    assert!(diag.z.x.abs() < EPSILON);
    assert!(diag.y.z.abs() < EPSILON);
    assert!(diag.z.y.abs() < EPSILON);
    assert_eq!(eigen.x.tolerate(EPSILON), diag.x.x);
    assert_eq!(eigen.y.tolerate(EPSILON), diag.y.y);
    assert_eq!(eigen.z.tolerate(EPSILON), diag.z.z);
    assert!((diag.x.x * rot.x).dist(m * &rot.x) < EPSILON);
    assert!((diag.y.y * rot.y).dist(m * &rot.y) < EPSILON);
    assert!((diag.z.z * rot.z).dist(m * &rot.z) < EPSILON);

    let id = rot * rot.transpose();
    assert!(id.x.dist([1,0,0].into()) < EPSILON);
    assert!(id.y.dist([0,1,0].into()) < EPSILON);
    assert!(id.z.dist([0,0,1].into()) < EPSILON);
}

#[test]
fn test_jacobi() {
    test_jacobi_single(&[[1, 2, 4], [2, 1, 3], [4, 3, 1]].into());
    test_jacobi_single(&[[10, 2, 0], [2, 0, 2], [0, 2, 10]].into());

    test_jacobi_single(&[[0.6375, 0., -1.7835], [0., 7.1568, 0.],
                         [-1.7835, 0., 1.4508]].into());

    test_jacobi_single(&[[0., -1.0298, 1.0792], [-1.0298, 1.2554, 0.],
                         [1.0792, 0., 0.1547]].into());
    test_jacobi_single(&[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());
}
