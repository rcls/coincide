
use std::cmp::Ordering;
use std::iter::*;
use std::ops::*;

#[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
pub struct Triple<T>{pub x: T, pub y: T, pub z: T}

impl<T: PartialEq<T> + PartialOrd<T>> Ord for Triple<T> {
    #[inline]
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
impl<T: PartialEq> Eq for Triple<T> { }

impl<T> Triple<T> {
    pub const fn new(x: T, y: T, z: T) -> Triple<T> { Triple{x, y, z} }
}

impl<U, T> From<[U; 3]> for Triple<T> where U: Into<T> {
    #[inline]
    fn from(r: [U; 3]) -> Triple<T> {
        let [x, y, z] = r;
        Triple::new(x.into(), y.into(), z.into())
    }
}

impl<T, U> Into<[U; 3]> for Triple<T> where T: Into<U> {
    #[inline]
    fn into(self) -> [U; 3] {[self.x.into(), self.y.into(), self.z.into()]}
}

impl<'a, T: 'a> IntoIterator for &'a Triple<T> {
    type Item = &'a T;
    type IntoIter = impl Iterator<Item = Self::Item>;
    fn into_iter(self) -> Self::IntoIter {
        [&self.x, &self.y, &self.z].into_iter()
    }
}

// Forward non-reference versions to the reference version.
macro_rules! binop_element {
    ($Trait:ident, $method:ident) => {
        impl<'a, T: 'a> const $Trait<&'a Triple<T>> for &'a Triple<T> where for<'u> &'u T: ~const $Trait<&'u T, Output=T> {
            type Output = Triple<T>;
            fn $method(self, r: &'a Triple<T>) -> Self::Output {
                Triple::new($Trait::$method(&self.x, &r.x),
                            $Trait::$method(&self.y, &r.y),
                            $Trait::$method(&self.z, &r.z))
            }
        }
        impl<'a, T> $Trait<Triple<T>> for &'a Triple<T>
            where for<'u> &'u T: $Trait<&'u T, Output=T>
        {
            type Output = Triple<T>;
            fn $method(self, r: Triple<T>) -> Triple<T> {
                $Trait::$method(self, &r)
            }
        }
        impl<'a, T> $Trait<&'a Triple<T>> for Triple<T> where for<'u> &'u T: $Trait<&'u T, Output=T> {
            type Output = Triple<T>;
            fn $method(self, r: &'a Triple<T>) -> Triple<T> {
                $Trait::$method(&self, r)
            }
        }
        impl<T> $Trait<Triple<T>> for Triple<T> where for<'u> &'u T: $Trait<&'u T, Output=T> {
            type Output = Triple<T>;
            fn $method(self, r: Triple<T>) -> Triple<T> {
                $Trait::$method(&self, &r)
            }
        }
    }
}

binop_element!{Add, add}
binop_element!{Sub, sub}

impl<U, T: AddAssign<U>> AddAssign<Triple<U>> for Triple<T> {
    fn add_assign(&mut self, r: Triple<U>) {
        self.x += r.x;
        self.y += r.y;
        self.z += r.z;
    }
}

// Multiply on left by scalar.
impl<T: Copy> Mul<Triple<T>> for f64 where f64: Mul<T> {
    type Output = Triple<<f64 as Mul<T>>::Output>;
    fn mul(self, t: Triple<T>) -> Self::Output {
        Triple::new(self * t.x, self * t.y, self * t.z)
    }
}

// Multiply on right by scalar.
impl<T: Mul<f64>> Mul<f64> for Triple<T> {
    type Output = Triple<<T as Mul<f64>>::Output>;
    fn mul(self, r: f64) -> Self::Output {
        Triple::new(self.x * r, self.y * r, self.z * r)
    }
}

impl<T> Div<f64> for Triple<T> where T: Div<f64> {
    type Output = Triple<<T as Div<f64>>::Output>;
    fn div(self, n: f64) -> Self::Output {
        Triple::new(self.x / n, self.y / n, self.z / n)
    }
}

impl<'a, T> Neg for &'a Triple<T> where &'a T: Neg {
    type Output = Triple<<&'a T as Neg>::Output>;
    fn neg(self) -> Self::Output {
        Triple::new(-&self.x, -&self.y, -&self.z)
    }
}

impl<T> Neg for Triple<T> where for<'a> &'a T: Neg<Output=T> {
    type Output = Triple<T>;
    fn neg(self) -> Self::Output { -&self }
}

impl<T> Sum for Triple<T> where T: Default, for<'a> &'a T: Add<&'a T, Output=T> {
    fn sum<I: Iterator<Item=Self>>(iter: I) -> Triple<T> {
        //fn sum<I: Iterator<Item=Self>>(iter: I) -> Triple<T> {
        let mut t = Default::default();
        for v in iter {
            t = t + v;
        }
        t
    }
}

impl<'a, T: 'a> Sum<&'a Triple<T>> for Triple<T> where
    T: 'a + Copy + Default,
    for<'u> &'u T: Add<&'u T, Output=T>
{
    fn sum<I: Iterator<Item=&'a Self>>(iter: I) -> Triple<T> {
        let mut t = Default::default();
        for v in iter {
            t = t + v;
        }
        t
    }
}
