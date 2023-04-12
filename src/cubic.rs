
#[cfg(test)]
use crate::test::Tolerate;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CubicSolution {
    Real(f64, f64, f64),
    Mixed(f64, f64, f64),
}

const COS60: f64 = 0.8660254037844386;

pub fn cubic_solve(a: f64, b: f64, c: f64, d: f64) -> CubicSolution {
    let delta0 = b * b - c * a * 3.0;
    let delta1 = b * b * b - c * a * b * 4.5 + d * a * a * 13.5;

    let p = delta1 * delta1 - delta0 * delta0 * delta0;

    let scale = 1.0 / -3.0 / a;

    // Solution is:
    // -(b + cubert(delta1 + sqrt(p)) + cubert_bar(delta1 - sqrt(p))) / 3a
    // where cubert & cubert_bar cycle through the cube roots in opposite
    // directions.
    if p > 0. {
        // One real, plus complex pair.
        let sqrtp = p.sqrt();
        // We can take either sign for the square root - algebraically, this
        // just swaps pls and mns.  We get the best numerical stability (and
        // avoid 0/0) if we take the one with the same sign as delta1.
        let pls = (delta1 + sqrtp.copysign(delta1)).cbrt();
        let mns = delta0 / pls;
        CubicSolution::Mixed(
            (b + (pls + mns)) * scale,
            (b - (pls + mns) * 0.5) * scale,
            ((pls - mns) * COS60 * scale).abs(),
        )
    }
    else {
        let sqrtp = (-p).sqrt();
        // Three real.  We want the cuberoot of delta1 + i * sqrtp.
        let angle = sqrtp.atan2(delta1) * (1. / 3.);
        let hypot = (delta1 * delta1 - p).sqrt().cbrt();
        let cs = angle.cos() * hypot;
        let sn = angle.sin() * hypot;
        // First root, the cube roots are (cs ± i sn), so their
        // sum is 2 * cs.
        let r1 = (b + 2.0 * cs) * scale;
        // Second root, the cube roots are (cs + i sn) * (i * cos60 - 0.5)
        // & complex conj.  For the third, use -cos60.
        let c2 = COS60 * 2. * sn;
        let r2 = (b - cs + c2) * scale;
        let r3 = (b - cs - c2) * scale;

        // Order the roots.
        let (u, v, w) = sort3a(r1, r2, r3);
        CubicSolution::Real(u, v, w)
    }
}

fn sort3a(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let (xa, ya, za) = (x.abs(), y.abs(), z.abs());
    if xa <= ya {
        if ya <= za {(x, y, z)} else if xa <= za {(x, z, y)} else {(z, x, y)}
    }
    else {
        if xa <= za {(y, x, z)} else if ya <= za {(y, z, x)} else {(z, y, x)}
    }
}


#[test]
fn test_cubert_1() {
    assert_eq!(COS60, 3f64.sqrt() / 2.);

    let c = cubic_solve(0., 0., 0., 0.);
    println!("{:?}", c.clone());

    let CubicSolution::Mixed(x, r, i) = cubic_solve(1., 0., 0., 1.)
    else { panic!() };
    assert_eq!(x.tolerate(1e-10), -1.);
    assert_eq!(r.tolerate(1e-10), 0.5);
    assert_eq!(i.tolerate(1e-10), COS60);

    let CubicSolution::Mixed(x2, r2, i2) = cubic_solve(0.5, 0., 0., 4.)
    else { panic!() };
    assert_eq!(x2, x * 2.);
    assert_eq!(r2, r * 2.);
    assert_eq!(i2, i * 2.);

    // (x - 1)(x + i)(x - i) = (x - 1)(x²+1) = x³ - x² + x - 1.
    let CubicSolution::Mixed(x, r, i) = cubic_solve(1., -1., 1., -1.)
    else { panic!() };
    assert_eq!(x.tolerate(1e-10), 1.);
    assert_eq!(r.tolerate(1e-10), 0.);
    assert_eq!(i.tolerate(1e-10), 1.);
}

#[test]
fn test_linear() {
    // (x - 1)(x - 2)(x - 3) = x³ - 6x² + 11x - 6
    let CubicSolution::Real(a, b, c) = cubic_solve(1., -6., 11., -6.)
    else { panic!() };
    assert_eq!(a.tolerate(1e-10), 1.);
    assert_eq!(b.tolerate(1e-10), 2.);
    assert_eq!(c.tolerate(1e-10), 3.);
}

#[test]
fn test_repeat() {
    // (x - 1)(x - 2)(x - 2) = x³ - 5x² + 8x - 4.
    let CubicSolution::Real(a, b, c) = cubic_solve(1., -5., 8., -4.)
    else { panic!() };
    assert_eq!(a, 1.);
    assert_eq!(b, 2.);
    assert_eq!(c, 2.);

    let CubicSolution::Real(a, b, c) = cubic_solve(1., -9., 27., -27.)
    else { panic!() };
    assert_eq!(a, 3.);
    assert_eq!(b, 3.);
    assert_eq!(c, 3.);
}

#[test]
fn test_sort3() {
    for x in -3..4i8 {
        for y in -3..4 {
            for z in -3..4 {
                let mut a = [x, y, z];
                a.sort_by_key(|p| p.abs());
                let s = sort3a(x as f64, y as f64, z as f64);
                assert_eq!((a[0].into(), a[1].into(), a[2].into()), s);
            }
        }
    }
}
