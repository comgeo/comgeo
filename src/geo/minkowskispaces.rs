use traits::*;

use std::fmt;

#[derive(Debug)]
pub struct LpSpace<R> {
    p: R
}
impl<R: Real> LpSpace<R> {
    pub fn new(p: R) -> Self {
        LpSpace { p: p }
    }

    pub fn p(&self) -> R {
        self.p
    }
}
impl<R: Real, P: Point<R=R>> MinkowskiSpace<P> for LpSpace<R> {
    fn norm(&self, p: &P) -> R {
        p.iter()
            .fold(P::R::zero(), |sum, &c|
                sum + c.abs().pow(self.p))
            .pow(self.p.recip())
    }
    fn dist(&self, p1: &P, p2: &P) -> R {
         p1.iter()
             .zip(p2.iter())
             .fold(P::R::zero(), |sum, (&c1, &c2)|
                sum + (c1 - c2).abs().pow(self.p))
             .pow(self.p.recip())
    }
}
impl<R: Real> fmt::Display for LpSpace<R> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Lp space with p={}", self.p)
    }
}
#[derive(Debug)]
pub struct LInfinity;
impl<P: Point> MinkowskiSpace<P> for LpSpace<LInfinity> {
    fn norm(&self, p: &P) -> P::R {
        p.iter().fold(P::R::zero(), |a, &c| a.abs().max(c))
    }
    fn dist(&self, p1: &P, p2: &P) -> P::R {
         p1.iter()
             .zip(p2.iter())
             .fold(P::R::zero(), |a, (&c1, &c2)|
                a.abs().max(c1-c2))
    }

}
impl fmt::Display for LpSpace<LInfinity> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Lp space with p=∞")
    }
}

#[derive(Debug)]
pub struct EuclideanSpace;
impl EuclideanSpace {
    pub fn new() -> Self {
        EuclideanSpace
    }
}
impl<P: Point> MinkowskiSpace<P> for EuclideanSpace {
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
impl fmt::Display for EuclideanSpace {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Euclidean space")
    }
}
