use traits::{Point, PruneTest};
use steinertree::{SteinerTree};


pub trait PruneTest {
    fn prunetest<P: Point>(&self, &SteinerTree<P>, P::R, P::R) -> bool;
}

#[derive(Debug)]
pub struct UpperBoundPruning;
impl PruneTest for UpperBoundPruning {
    fn prunetest<P: Point>(&self, _: &SteinerTree<P>, len: P::R, ub: P::R)-> bool {
        len >= ub
    }
}
