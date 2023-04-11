
pub fn ascend([a, b, c, d]: [u8; 4]) -> [u8; 4] {
    let (a, b) = if a < b { (a, b) } else { (b, a) };
    let (c, d) = if c < d { (c, d) } else { (d, c) };
    if a < c {
        if b < c {
            [a, b, c, d]
        } else if b < d {
            [a, c, b, d]
        } else {
            [a, c, d, b]
        }
    }
    else {
        if b < d {
            [c, a, b, d]
        }
        else if a < d {
            [c, a, d, b]
        }
        else {
            [c, d, a, b]
        }
    }
}

pub fn join(u: [u8; 4]) -> u32 { u32::from_ne_bytes(u) }
#[cfg(test)]
pub fn split(u: u32) -> [u8; 4] { u.to_ne_bytes() }

pub fn join_ascend(v: [u8; 4]) -> u32 { join(ascend(v)) }

#[test]
fn test_ascend() {
    for a in 0..10 {
        for b in 0..10 {
            for c in 0..10 {
                for d in 0..10 {
                    let u = [a, b, c, d];
                    let v = ascend(u);
                    let mut s = u;
                    s.sort();
                    assert_eq!(s, v);
                    assert_eq!(join_ascend(u), join(s));
                }
            }
        }
    }
}

#[test]
fn test_join_split() {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    for _ in 0..1000 {
        let v: u32 = rng.gen();
        assert_eq!(join(split(v)), v);
    }
}

#[test]
fn test_split_join() {
    use rand::Rng;
    let mut rng = rand::thread_rng();
    for _ in 0..1000 {
        let v: [u8; 4] = rng.gen();
        assert_eq!(split(join(v)), v);
    }
}
