mod lp;
mod euclidean;
mod hyper_ellipsoid;

pub use self::euclidean::Euclidean;
pub use self::lp::Lp;
pub use self::hyper_ellipsoid::HyperEllipsoid;

pub trait MinkowskiSpace<P: Point>: fmt::Display {
    fn norm(&self, &P) -> P::R;
    fn dist(&self, p1: &P, p2: &P) -> P::R {
        self.norm(p1.clone().sub(p2))
    }
}