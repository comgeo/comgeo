use std::fmt;



#[derive(Debug)]
pub struct Euclidean;
impl EuclideanSpace {
    pub fn new() -> Self {
        EuclideanSpace
    }
}
impl<P: Point> MinkowskiSpace<P> for Euclidean {
    fn norm(&self, p: &P) -> P::R {
        p.iter().fold(P::R::zero(), |sum, &c| sum + c * c).sqrt()
    }
    fn dist(&self, p1: &P, p2: &P) -> P::R {
        p1.iter()
            .zip(p2.iter())
            .fold(P::R::zero(), |sum, (&c1, &c2)|
                sum + (c1 - c2) * (c1 - c2)
            ).sqrt()
    }
}
impl fmt::Display for Euclidean {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Euclidean space")
    }
}
