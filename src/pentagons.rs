
use crate::gold::IcoGroup;
use crate::sortn::SortN;

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
    // [0, 3, x, 4, y] is reverse of [0,y,4,x,3]
    // [0, 3, 4, x, y] is reverse of [0,y,x,4,3]
    // [0, 4, x, y, z] is reverse of [0,z,y,x,4]
];

fn edges(pentagon: &Pentagon, a: u8, b: u8, c: u8, d: u8) -> [(u8, u8); 5] {
    let indexes = [0, a, b, c, d];
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
    let to = |x: u8, y: u8| group.mult_index[y as usize]
        [group.inv_index[x as usize] as usize];
    for (r, s) in edges(pentagon, a, b, c, d) {
        if i != r && to(i, r) == to(j, s) {
            return true;
        }
        if i != s && to(i, s) == to(j, r) {
            return true;
        }
    }
    false
}

fn closes(group: &IcoGroup, pentagon: &Pentagon,
          a: u8, b: u8, c: u8, d: u8) -> bool {
    for (i, j) in edges(pentagon, a, b, c, d) {
        if !edge_closes(group, pentagon, a, b, c, d, i, j) {
            return false;
        }
    }
    true
}

pub fn closures(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8) -> usize {
    // Reject items with a non trivial permutation from the group.  The
    // permutation has to map 0 to a point, so the permutation is given
    // by one of the points.
    let orig = [a, b, c, d].sortn();
    for p in orig {
        let mut hit = false;
        let mut image = orig;
        for i in 0..4 {
            let x = &mut image[i];
            *x = group.mult_index[p as usize][*x as usize];
            if *x == 0 {
                hit = true;
                *x = p;
            }
        }
        if hit && orig == image.sortn() {
            // println!("Automorphism");
            return 0;
        }
    }
    // Check that each edge is the non-trivial group image of an edge.
    let mut count = 0;
    for p in &PENTAGONS {
        if closes(group, p, a, b, c, d) {
            println!("Pentagon {:?}", p);
            count += 1;
        }
    }
    // Check for compounds....
    if count != 0 {
        connected(group, a, b, c, d);
    }
    count
}


fn connected(group: &IcoGroup, a: u8, b: u8, c: u8, d: u8) -> bool {
    let mut queue = Vec::new();
    let mut seen = std::collections::HashSet::new();
    queue.push(0);
    seen.insert(0);
    let p = [a, b, c, d];
    while let Some(i) = queue.pop() {
        for j in p {
            let jp = group.mult_index[j as usize][i as usize];
            if seen.insert(jp) {
                queue.push(jp);
            }
        }
    }
    println!("Component {}", seen.len());
    seen.len() == 120
}
