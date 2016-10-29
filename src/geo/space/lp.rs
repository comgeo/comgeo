#[derive(Debug)]
pub struct Lp<R> {
    p: R
}
impl<R: Real> Lp<R> {
    pub fn new(p: R) -> Self {
        Lp { p: p }
    }

    pub fn p(&self) -> R {
        self.p
    }
}
impl<R: Real, P: Point<R=R>> MinkowskiSpace<P> for Lp<R> {
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
impl<R: Real> fmt::Display for Lp<R> {
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