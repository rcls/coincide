#![deny(warnings)]
#![feature(const_trait_impl)]
#![feature(type_alias_impl_trait)]

mod cubic;
mod gold;
mod jacobi;
mod matrix;
mod quad;
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
        let key = join_ascend([z, p(b), p(c), p(d)]);
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

pub fn evaluate(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8,
            x: f64, y: f64) -> f64 {
    let epsilon = 0.1;
    let z = 1. - (x + y);
    let p = group.corners * Vector{x, y, z};
    //let (_, _, mut r, _) = plane(&five_points(group, &p, a, b, c, d));
    let mut r = covar(&five_points(group, &p, a, b, c, d))
        .least_abs_eigenvalue();
    if x < -epsilon {
        r += x * x - epsilon * epsilon;
    }
    if y < -epsilon {
        r += y * y - epsilon * epsilon;
    }
    if z < -epsilon {
        r += z * z - epsilon * epsilon;
    }
    r
}

fn derivs(f: impl Fn(f64, f64) -> f64, x:f64, y: f64) -> [f64; 5] {
    #[allow(mixed_script_confusables)]
    let ε = 0.1;
    let mut m = [[0.; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[j][i] = f(x + (i as f64 - 1.) * 0.1, y + (j as f64 - 1.) * 0.1);
        }
    }
    // println!("{:?}", m);
    let m: crate::matrix::Matrix = m.into();
    let dx = m * Vector::new(-5./3., 0., 5./3.) * Vector::new(1., 1., 1.);
    let dy = m * Vector::new(1., 1., 1.) * Vector::new(-5./3., 0., 5./3.);
    let dxx = m * Vector::new(100./3., -200./3., 100./3.) * Vector::new(1., 1., 1.);
    // println!("{:?}", m * Vector::new(1., 1., 1.));
    let dyy = m * Vector::new(1., 1., 1.) * Vector::new(100./3., -200./3., 100./3.);
    let dxy = m * Vector::new(-25., 0., 25.) * Vector::new(-1., 0., 1.);
    [dx, dy, dxx, dyy, dxy]
}

pub fn step(f: impl Fn(f64, f64) -> f64, x: f64, y: f64)
        -> Option<(f64, f64)>
{
    let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y);

    // We will move in the -(dx,dy) direction
    // f(x - t dx, y - t dy)
    // = f(x, y) - t dx² - t dy²
    // + 0.5 t² dx² dxx + t² dx dy dxy + 0.5 t² dyy
    //
    // Find t to minimize the approximation....
    let numerator = dx * dx + dy * dy;
    let denominator = dx * dx * dxx + 2. * dx * dy * dxy + dy * dy * dyy;
    // We will have delta_x = -dx * num / den, delta_y = -dy * num / den

    // If delta_x, delta_y are very small, we terminate.
    // If denominator is negative or blows up the fraction, then fail.
    let num_norm = numerator * numerator.sqrt();
    if num_norm <= denominator {
        Some((-numerator / denominator * dx, -numerator / denominator * dy))
    }
    else {
        None
    }
}

fn covar(vv: &[Vector]) -> Matrix {
    //let mean: Vector = vv.iter().sum() / vv.len();
    let mut sum = Vector::default();
    for &v in vv {
        sum += v;
    }
    let mean = sum / vv.len() as f64;
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
                        //plane::<true>(&p);
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
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y);
            assert_eq!(dx.tolerate(1e-12), y);
            assert_eq!(dy.tolerate(1e-12), x);
            assert_eq!(dxx, 0.);
            assert_eq!(dyy.tolerate(3e-11), 0.);
            assert_eq!(dxy, 1.);
        }
    }
}

#[test]
fn test_derivs_2() {
    fn f(x: f64, y: f64) -> f64 { x * x + y * y }
    for x in [-2., -1., 0., 1., 2.] {
        for y in [-2., -1., 0., 1., 2.] {
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y);
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
            let [dx, dy, dxx, dyy, dxy] = derivs(f, x, y);
            assert_eq!(dx.tolerate(3e-15), 2. * x);
            assert_eq!(dy.tolerate(4e-15), -2. * y);
            assert_eq!(dxx.tolerate(6e-14), 2.);
            assert_eq!(dyy.tolerate(6e-14), -2.);
            assert_eq!(dxy.tolerate(2e-14), 0.);
        }
    }
}
