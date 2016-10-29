use steinertree::{SteinerTree};

extern crate disjoint_set;
use self::disjoint_set::DisjointSet;

#[derive(Debug)]
pub struct Kruskal { }

#[derive(Eq, PartialEq, PartialOrd)]
struct Edge<R: Real> {
    len: R,
    from: usize,
    to: usize
}


impl Kruskal {
    pub fn new() -> Self {
        Kruskal {}
    }
}

impl<P: Point, M: MinkowskiSpace<P>> MST<P, M> for Kruskal {
    fn find(&mut self, terminals: &[P], geo: &M) -> SteinerTree<P> {
        let mut forest = DisjointSet::new();
        for i in 0..terminals.len() {
            forest.make_set(i);
        }

        let mut edges = Vec::with_capacity(terminals.len()*terminals.len());
        let mut st = SteinerTree::new(terminals, &[], &[]);

        for i in 0..terminals.len() - 1 {
            for n in i+1..terminals.len() {
                edges.push(Edge::<P::R> {
                    len: geo.dist(&terminals[i], &terminals[n]),
                    from: i,
                    to: n
                });
            }
        }

        edges.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mut c = 0;
        for e in edges.iter() {
            if forest.find(e.from) != forest.find(e.to) {
                st.add_edge((e.from, e.to));
                forest.union(e.from, e.to).unwrap();
                c += 1;
                if c == terminals.len() - 1 {
                    return st;
                }
            }
        };

        unreachable!("!")
    }
}