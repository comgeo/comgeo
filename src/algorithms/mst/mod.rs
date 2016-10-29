pub mod kruskal;

pub trait MST<P: Point, M: MinkowskiSpace<P>> {
    fn find(&mut self, &[P], geo: &M) -> SteinerTree<P>;
}