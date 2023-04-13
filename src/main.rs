#![deny(warnings)]
#![allow(mixed_script_confusables)]
#![feature(const_trait_impl)]
#![feature(type_alias_impl_trait)]

mod cubic;
mod gold;
mod jacobi;
mod matrix;
mod quad;
mod sortn;
mod test;
mod triple;
mod vector;

use gold::IcoGroup;
use matrix::Matrix;
use vector::Vector;

#[cfg(test)]
use test::Tolerate;

use std::collections::HashMap;

fn plane(vv: &[Vector]) -> (Vector, f64, f64, u8) {
    let (evectors, evalues) = jacobi::jacobi(&covar(vv));
    // Find the smallest eigen value and it's eigen vector.
    let (normal, residual) =
        if evalues.x < evalues.y && evalues.x < evalues.z {
            (evectors.x, evalues.x)
        } else if !(evalues.x < evalues.y) && evalues.y < evalues.z {
            (evectors.y, evalues.y)
        } else {
            (evectors.z, evalues.z)
        };
    let mean = vv.iter().sum::<Vector>() / vv.len() as f64;
    let distance = normal * mean;
    if distance >= 0. {
        (normal, distance, residual, 1)
    }
    else {
        (-normal, -distance, residual, 1)
    }
}


fn is_done(group: &IcoGroup, stuff: &mut HashMap<u32, (f64, u8)>,
           [a, b, c, d]: [u8; 4]) -> bool
{
    use quad::join_ascend;
    let mut one = |a, b, c, d| {
        let z = group.inv_index[a as usize];
        let p = |x| group.mult_index[z as usize][x as usize];
        let key = join_ascend(&[z, p(b), p(c), p(d)]);
        if let Some(q) = stuff.get_mut(&key) {
            q.1 += 1;
            true
        }
        else {
            false
        }
    };
    one(a, b, c, d) || one(b, a, c, d) || one(c, a, b, d) || one(d, a, b, c)
}

fn five_points(group: &IcoGroup, p: &Vector, a: u8, b: u8, c: u8, d: u8)
               -> [Vector; 5] {
    let f = |i| group.map(i) * p;
    [*p, f(a), f(b), f(c), f(d)]
}

fn evaluate(group: &IcoGroup, f: impl FnOnce(&Matrix) -> f64,
            a: u8, b: u8, c: u8, d: u8, x: f64, y: f64) -> f64 {
    let ε = 0.1;
    let z = 1. - (x + y);
    // It doesn't really matter if we normalise to a unit vector.  Not
    // normalising keeps the covariance matrix as a quad. form in x and y.
    let p = group.corners.transpose() * Vector{x, y, z};
    let mut r = f(&covar(&five_points(group, &p, a, b, c, d)));
    if x < -ε {
        r += x * x - ε * ε;
    }
    if y < -ε {
        r += y * y - ε * ε;
    }
    if z < -ε {
        r += z * z - ε * ε;
    }
    r
}

fn evaluate_eigen(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8,
                  x: f64, y: f64) -> f64 {
    evaluate(group, Matrix::least_abs_eigenvalue, a, b, c, d, x, y)
}

fn evaluate_determinant(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8,
                        x: f64, y: f64) -> f64 {
    evaluate(group, Matrix::determinant, a, b, c, d, x, y)
}

fn derivs(f: impl Fn(f64, f64) -> f64, x:f64, y: f64, ε: f64) -> [f64; 5] {
    let c = 1. / ε;
    let mut m = [[0.; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[j][i] = f(x + (i as f64 - 1.) * ε, y + (j as f64 - 1.) * ε);
        }
    }
    let c6 = c * 1./6.;
    let c3 = 1./3. * c * c;
    let c4 = 0.25 * c * c;
    let m: Matrix = m.into();
    let dx = m * Vector::new(-1., 0., 1.) * Vector::new(1., 1., 1.) * c6;
    let dy = m * Vector::new(1., 1., 1.) * Vector::new(-1., 0., 1.) * c6;
    let dxx = m * Vector::new(1., -2., 1.) * Vector::new(1., 1., 1.) * c3;
    let dyy = m * Vector::new(1., 1., 1.) * Vector::new(1., -2., 1.) * c3;
    let dxy = m * Vector::new(-1., 0., 1.) * Vector::new(-1., 0., 1.) * c4;
    [dx, dy, dxx, dyy, dxy]
}

pub fn step(f: impl Fn(f64, f64) -> f64, x: f64, y: f64, ε: f64)
            -> Option<(f64, f64)> {
    let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y, ε);

    // We want a, b to make f(x+a, f(x+b) stationary.
    //
    // f_x(x+a, y+b) ~ dx + a·dxx + b·dxy
    // f_y(x+a, y+b) ~ dy + a·dxy + b·dyy
    //
    // so we want [dxx,dxy;dxy,dyy][a;b] + [dx;dy] = 0.
    //
    // This is only valid if the matrix is positive-definite.  It is symmetric
    // so we just need positive determinant.
    let det = dxx * dyy - dxy * dxy;
    let vx = dxy * dy - dyy * dx;
    let vy = dxy * dx - dxx * dy;
    // Don't do big jumps.  For us, we start in the middle of [1,0],[0,1],[0,0]
    // so bounding the step 1 on each coordinate makes sense.
    if vx.abs() < det && vy.abs() < det && dxx + dyy >= 0. {
        // println!("Easy!");
        return Some((vx / det, vy / det));
    }
    // Ok, we need more examination.  Do an eigenvalue decomposition of the
    // Hessian.  In some circumstances one eigenvalue may be small (or negative)
    // and unusable, but with the other usable.
    let (c, s) = jacobi::jacobi2(dxx, dyy, dxy);
    // rᵀ H r should now be diagonal with r = [c,-s; s,c].
    //
    // This gives the eigenvalues:
    //
    //     c²·dxx + 2·c·s·dxy + s²·dyy  and  s²·dxx - 2·c·s·dxy c²·dyy
    //
    // with eigenvectors: [c;s] and [-s;c] respectively.
    //
    // We limit our attention to the larger eigenvalue and go towards the
    // minimum in the direction of its eigenvector.
    let λ1 = c * c * dxx + c * s * dxy * 2.0 + s * s * dyy;
    let λ2 = s * s * dxx - c * s * dxy * 2.0 + c * c * dyy;
    let (λ, x, y) = if λ1 > λ2 {(λ1, c, s)} else {(λ2, -s, c)};
    let d = x * dx + y * dy;
    if d.abs() < λ {
        Some((-d / λ * x, -d / λ * y))
    }
    else {
        None
    }
}

fn covar(vv: &[Vector]) -> Matrix {
    let mean: Vector = vv.iter().sum::<Vector>() / vv.len() as f64;
    let mut covar = [[0.; 3]; 3];
    for v in vv {
        let v: [f64; 3] = (v - mean).into();
        for (i, x) in v.iter().enumerate() {
            for (j, y) in v.iter().enumerate() {
                covar[i][j] += x * y;
            }
        }
    }
    covar.into()
}

fn main() {
    let group = IcoGroup::new();
    let point = Vector::new(1., 1., 4. * gold::GOLD + 1.);

    let mut residuals = HashMap::new();
    for a in 1..120 {
        for b in 1..a {
            for c in 1..b {
                for d in 1..c {
                    if is_done(&group, &mut residuals, [d, c, b, a]) {
                        continue;
                    }
                    let p = five_points(&group, &point, d, c, b, a);
                    //let fit = plane::<false>(&p);
                    let u = covar(&p).least_abs_eigenvalue();
                    residuals.insert(quad::join([d, c, b, a]), (u, 1));
                    if true &&
                        (u > 5. && u < 5.00005 || u > 20. && u < 20.00006)
                    {
                        println!("------------------------------------");
                        println!("Points: {:?}", p);
                        println!("Eigen: {:?}", covar(&p).eigenvalues());
                        println!("Fit: {:?}", plane(&p));
                    }
                }
            }
        }
    }
    let mut max = 0.;
    for (_, &(r, _)) in &residuals {
        if r > max {
            max = r;
        }
    }
    println!("{} {:.6} {:.6}", residuals.len(), max, (max * 0.2).sqrt());
    let mut bins = Vec::new();
    let mut counts = [0; 6];
    bins.resize((max * 0.2e12).sqrt() as usize + 1, 0u16);
    for (k, &(r, n)) in &residuals {
        bins[(r * 0.2e12).sqrt() as usize] += n as u16;
        counts[n as usize] += 1;
        if r > 180. {
            println!("Bad {:#08x} {}", k, (r * 0.2).sqrt());
        }
    }
    println!("{:?}", counts);
    println!("{:?}", &bins[0..10]);
    if false {
        for (bin, &count) in bins.iter().enumerate() {
            if count > 0 {
                println!("{:.6} {}", bin as f64 * 1e-6, count);
            }
        }
    }
}

#[test]
fn test_derivs_1() {
    fn f(x: f64, y: f64) -> f64 { x * y + 1000. }
    for x in [-2., -1., 0., 1., 2.] {
        for y in [-2., -1., 0., 1., 2.] {
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y, 0.1);
            assert_eq!(dx.tolerate(1e-12), y);
            assert_eq!(dy.tolerate(2e-12), x);
            assert_eq!(dxx, 0.);
            assert_eq!(dyy.tolerate(3e-11), 0.);
            assert_eq!(dxy.tolerate(1e-12), 1.);
        }
    }
}

#[test]
fn test_derivs_2() {
    fn f(x: f64, y: f64) -> f64 { x * x + y * y }
    for x in [-2., -1., 0., 1., 2.] {
        for y in [-2., -1., 0., 1., 2.] {
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y, 0.1);
            assert_eq!(dx.tolerate(1e-12), 2. * x);
            assert_eq!(dy.tolerate(1e-12), 2. * y);
            assert_eq!(dxx.tolerate(2e-13), 2.);
            assert_eq!(dyy.tolerate(2e-13), 2.);
            assert_eq!(dxy.tolerate(3e-14), 0.);
        }
    }
}

#[test]
fn test_derivs_3() {
    fn f(x: f64, y: f64) -> f64 { x * x - y * y }
    for x in [-2., -1., 0., 1., 2.] {
        for y in [-2., -1., 0., 1., 2.] {
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y, 0.1);
            assert_eq!(dx.tolerate(3e-15), 2. * x);
            assert_eq!(dy.tolerate(4e-15), -2. * y);
            assert_eq!(dxx.tolerate(6e-14), 2.);
            assert_eq!(dyy.tolerate(6e-14), -2.);
            assert_eq!(dxy.tolerate(2e-14), 0.);
        }
    }
}

#[test]
fn test_step() {
    // (x-y)² + 2(x+y)²
    fn f(x: f64, y: f64) -> f64 { 3. * (x * x + y * y) + 2. * x * y }
    let (dx, dy) = step(f, 0.3, 0.2, 0.1).unwrap();
    assert_eq!(dx.tolerate(7e-16), -0.3);
    assert_eq!(dy.tolerate(5e-16), -0.2);

    assert_eq!(step(|x, y| -x*x - y*y, 0., 0., 0.1), None);
}

#[test]
fn test_step_saddle() {
    let sq = |x| x * x;
    let (dx, dy) = step(|x, y| sq(x) - 20. * sq(y), 1., 1., 0.1).unwrap();
    assert_eq!(dx.tolerate(5e-14), -1.);
    assert_eq!(dy.tolerate(3e-15), 0.);

    let (dx, dy) = step(|x, y| sq(x - y) - 20. * sq(x + y), 1., 0., 0.1).unwrap();
    assert_eq!(dx.tolerate(4e-14), -0.5);
    assert_eq!(dy.tolerate(4e-14), 0.5);
}
