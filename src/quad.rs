
use crate::sortn::SortN;

pub fn join(u: [u8; 4]) -> u32 { u32::from_ne_bytes(u) }
pub fn split(u: u32) -> [u8; 4] { u.to_ne_bytes() }

pub fn join_ascend(v: &[u8; 4]) -> u32 { join(v.sortn()) }

#[test]
fn test_join_ascend() {
    for a in 0..10 {
        for b in 0..10 {
            for c in 0..10 {
                for d in 0..10 {
                    let u = [a, b, c, d];
                    let v = u.sortn();
                    let mut s = u;
                    s.sort();
                    assert_eq!(s, v);
                    assert_eq!(join_ascend(&u), join(s));
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
