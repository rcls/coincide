
use crate::gold::IcoGroup;

type Pentagon = [usize; 5];

const PENTAGONS: [[usize; 5]; 12] = [
    [0, 1, 2, 3, 4],
    [0, 1, 2, 4, 3],
    [0, 1, 3, 2, 4],
    [0, 1, 3, 4, 2],
    [0, 1, 4, 2, 3],
    [0, 1, 4, 3, 2],
    [0, 2, 1, 3, 4],
    [0, 2, 1, 4, 3],
    [0, 2, 3, 1, 4],
    // [0, 2, 3, 4, 1] is reverse of [0,1,4,3,2]
    [0, 2, 4, 1, 3],
    // [0, 2, 4, 3, 1] is reverse of [0,1,3,4,2]
    [0, 3, 1, 2, 4],
    [0, 3, 2, 1, 4],
    // [0, 3, 4, x, y] is reverse of [0,y,x,4,3]
    // [0, 4, x, y, z] is reverse of [0,z,y,x,4]
];

fn edges(group: &IcoGroup, pentagon: &Pentagon,
         auto: usize, a: u8, b: u8, c: u8, d: u8) -> [(u8, u8); 5] {
    let perm = |n| group.mult_index[auto][n as usize];
    let indexes = [perm(0), perm(a), perm(b), perm(c), perm(d)];
    [
        (indexes[pentagon[0]], indexes[pentagon[1]]),
        (indexes[pentagon[1]], indexes[pentagon[2]]),
        (indexes[pentagon[2]], indexes[pentagon[3]]),
        (indexes[pentagon[3]], indexes[pentagon[4]]),
        (indexes[pentagon[4]], indexes[pentagon[0]])
    ]
}

fn edge_closes(group: &IcoGroup, pentagon: &Pentagon,
               a: u8, b: u8, c: u8, d: u8, i: u8, j: u8) -> bool {
    for p in 1..120 {
        for (r, s) in edges(group, pentagon, p, a, b, c, d) {
            if (i == r && j == s) || (j == s && i == r) {
                return true;
            }
        }
    }
    false
}

fn closes(group: &IcoGroup, pentagon: &Pentagon,
          a: u8, b: u8, c: u8, d: u8) -> bool {
    for (i, j) in edges(group, pentagon, 0, a, b, c, d) {
        if !edge_closes(group, pentagon, a, b, c, d, i, j) {
            return false;
        }
    }
    true
}

pub fn closures(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8) -> usize {
    let mut count = 0;
    for p in &PENTAGONS {
        if closes(group, p, a, b, c, d) {
            println!("Pentagon {:?}", p);
            count += 1;
        }
    }
    count
}
