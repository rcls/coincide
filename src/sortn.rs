
pub trait SortN {
    type E;
    fn sortn_by_key<K: PartialOrd>(&self, f: impl Fn(&Self::E) -> K) -> Self;
    fn sortn(&self) -> Self where Self: Sized, Self::E: Copy + PartialOrd {
        self.sortn_by_key(|x| *x)
    }
}

impl<T: Copy> SortN for (T, T, T) {
    type E = T;
    fn sortn_by_key<K: PartialOrd>(&self, f: impl Fn(&Self::E) -> K) -> Self {
        let (a, b, c) = self;
        if f(a) <= f(b) {
            if f(b) <= f(c) {(*a, *b, *c)}
            else if f(a) <= f(c) {(*a, *c, *b)}
            else {(*c, *a, *b)}
        }
        else {
            if f(a) <= f(c) {(*b, *a, *c)}
            else if f(b) <= f(c) {(*b, *c, *a)}
            else {(*c, *b, *a)}
        }
    }
}

impl<T: Copy> SortN for [T; 4] {
    type E = T;
    fn sortn_by_key<K: PartialOrd>(&self, f: impl Fn(&Self::E) -> K) -> Self {
        let [a, b, c, d] = self;
        let (a, b) = if f(a) <= f(b) {(*a, *b)} else {(*b, *a)};
        let (c, d) = if f(c) <= f(d) {(*c, *d)} else {(*d, *c)};
        if f(&a) <= f(&c) {
            if f(&b) <= f(&c) {
                [a, b, c, d]
            } else if f(&b) <= f(&d) {
                [a, c, b, d]
            } else {
                [a, c, d, b]
            }
        }
        else {
            if f(&b) <= f(&d) {
                [c, a, b, d]
            }
            else if f(&a) <= f(&d) {
                [c, a, d, b]
            }
            else {
                [c, d, a, b]
            }
        }
    }
}
