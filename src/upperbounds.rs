use traits::*;
use steinertree::{SteinerTree};
use algorithms::mst::Kruskal;

#[derive(Debug)]
pub struct LineTree;
impl<P, M> UpperBound<P, M> for LineTree
    where P: Point, M: MinkowskiSpace<P> {

    fn bound(&self, t: Vec<P>, _: &M) -> SteinerTree<P> {
        let mut edges = Vec::new();
        let mut i = 0;
        while i < t.len() - 1 {
            edges.push((i, i+1));
            i += 1;
        }

        SteinerTree::new(&t[..], &[], &edges[..])
    }
}

impl Default for LineTree {
    fn default() -> Self {
        LineTree
    }
}

#[derive(Debug)]
pub struct MSTBound;
impl<P: Point, M: MinkowskiSpace<P>> UpperBound<P, M> for MSTBound {
    fn bound(&self, t: Vec<P>, geo: &M) -> SteinerTree<P> {
        Kruskal::new().find(&t[..], geo)
    }
}
