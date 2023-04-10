
use std::ops::{Add, Sub};

pub trait Tolerate {
    fn tolerate<T>(self, tolerance: T) -> Tolerance<Self, T> where Self: Sized {
        Tolerance{value: self, tolerance}
    }
}

impl<T> Tolerate for T { }

#[derive(Debug)]
pub struct Tolerance<V, T> {
    value: V,
    tolerance: T
}

impl<V, T: Copy, U> PartialEq<U> for Tolerance<V, T> where
    V: Add<T> + Sub<T> + Copy,
    <V as std::ops::Add<T>>::Output: PartialOrd<U>,
    <V as std::ops::Sub<T>>::Output: PartialOrd<U>,
{
    fn eq(&self, other: &U) -> bool {
        self.value - self.tolerance <= *other
            && self.value + self.tolerance >= *other
    }
}
