
pub trait UpperBound<P: Point, M: MinkowskiSpace<P>> {
    fn bound(&self, Vec<P>, geo: &M) -> SteinerTree<P>;
}