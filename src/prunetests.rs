use traits::{Point, PruneTest};
use steinertree::{SteinerTree};

#[derive(Debug)]
pub struct UpperBoundPruning;
impl PruneTest for UpperBoundPruning {
    fn prunetest<P: Point>(&self, _: &SteinerTree<P>, len: P::R, ub: P::R)-> bool {
        len >= ub
    }
}
