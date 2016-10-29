use traits::*;
use steinertree::{SteinerTree};
use algorithms::mst::Kruskal;




#[derive(Debug)]
pub struct MSTBound;
impl<P: Point, M: MinkowskiSpace<P>> UpperBound<P, M> for MSTBound {
    fn bound(&self, t: Vec<P>, geo: &M) -> SteinerTree<P> {
        Kruskal::new().find(&t[..], geo)
    }
}
