use geo::point::{Point2d, Point3d, PointNd};
use geo::hyperplane::{Line, Plane, HyperplaneNd};

use std::slice::Iter;
use std::fmt;


pub trait HyperEllipsoid<P: Point> : MinkowskiSpace<P> {
    type H: Hyperplane<P>;

    fn comps(&self) -> &[P::R];
    fn tangent(&self, o: &P, p: &P) -> Self::H;

    fn bd_intersect<'a>(&self, p: &'a mut P) -> &'a mut P {
        let norm = self.norm(p);
        p.div(norm)
    }

    fn inner_product(&self, p1: &P, p2: &P) -> P::R {
        self.comps().iter().zip(p1.iter().zip(p2.iter()))
            .fold(P::R::zero(), |sum, (&a, (&x, &y))|
                sum + (a*a).recip() * x * y
            )
    }
}

fn ellipsoid_norm<R: Real>(e: Iter<R>, p: Iter<R>) -> R {
    e.zip(p)
     .fold(R::zero(), |sum, (&a, &x)|
         sum + ((x*x) / (a*a))
     ).sqrt()
}

#[derive(Debug)]
pub struct HyperEllipsoidSpaceNd<R> {
    c: Vec<R>
}
impl<R> HyperEllipsoidSpaceNd<R> {
    pub fn new(c: Vec<R>) -> Self {
        HyperEllipsoidSpaceNd { c: c }
    }
}
impl<R: Real> HyperEllipsoidSpace<PointNd<R>> for HyperEllipsoidSpaceNd<R> {
    type H = HyperplaneNd<R>;

    fn tangent(&self, o: &PointNd<R>, p: &PointNd<R>) -> HyperplaneNd<R> {
        let mut plane = HyperplaneNd::new(
            self.c.iter().zip(self.bd_intersect(p.clone().sub(o)).iter())
                .map(|(&a, &bdx)|
                    bdx/(a*a)
                ).collect());

        plane.translate(o);
        plane
    }

    fn comps(&self) -> &[R] {
        &self.c[..]
    }
}
impl<R: Real> MinkowskiSpace<PointNd<R>> for HyperEllipsoidSpaceNd<R> {
    fn norm(&self, p: &PointNd<R>) -> R {
        ellipsoid_norm(self.c.iter(), p.iter())
    }
}
impl<R: Real> fmt::Display for HyperEllipsoidSpaceNd<R> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "Hyperellipsoid space with unit ball {}", self.c[0]));
        for c in self.c[1..].iter() {
            try!(write!(f, " + {}", c));
        }
        write!(f, " = 1")
    }
}

#[derive(Debug)]
pub struct EllipseSpace<R> {
    comp: [R; 2]
}
impl<R> EllipseSpace<R> {
    pub fn new(a: R, b: R) -> Self {
        EllipseSpace { comp: [a, b] }
    }
}
impl<R: Real> HyperEllipsoidSpace<Point2d<R>> for EllipseSpace<R> {
    type H = Line<R>;

    fn tangent(&self, o: &Point2d<R>, p: &Point2d<R>) -> Line<R> {
        let mut line = Line::new(
            self.bd_intersect(p.clone().sub(o))
                .modify(&Point2d::new(self.comp.clone()), &|x, a|
                    x / (a*a))
                .arr().clone());
        line.translate(o);
        line
    }

    fn comps(&self) -> &[R] {
        &self.comp
    }
}
impl<R: Real> MinkowskiSpace<Point2d<R>> for EllipseSpace<R> {
    fn norm(&self, p: &Point2d<R>) -> R {
        ellipsoid_norm([self.comp[0], self.comp[1]].iter(), p.iter())
    }
}
impl<R: Real> fmt::Display for EllipseSpace<R> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "Ellipse space with unit ball {}", self.comp[0]));
        for c in self.comp[1..].iter() {
            try!(write!(f, " + {}", c));
        }
        write!(f, " = 1")
    }
}

#[derive(Debug)]
pub struct EllipsoidSpace<R> {
    comp: [R; 3]
}
impl<R> EllipsoidSpace<R> {
    pub fn new(a: R, b: R, c: R) -> Self {
        EllipsoidSpace { comp: [a, b, c] }
    }
}
impl<R: Real> HyperEllipsoidSpace<Point3d<R>> for EllipsoidSpace<R> {
    type H = Plane<R>;

    fn tangent(&self, o: &Point3d<R>, p: &Point3d<R>) -> Plane<R> {
        let mut plane = Plane::new(
            self.bd_intersect(p.clone().sub(o)).modify(
                &Point3d::new(self.comp.clone()), &|x, a|
                x / (a*a)
            ).arr().clone());
        plane.translate(o);
        plane
    }

    fn comps(&self) -> &[R] {
        &self.comp
    }
}
impl<R: Real> MinkowskiSpace<Point3d<R>> for EllipsoidSpace<R> {
    fn norm(&self, p: &Point3d<R>) -> R {
        ellipsoid_norm([self.comp[0], self.comp[1], self.comp[2]].iter(), p.iter())
    }
}
impl<R: Real> fmt::Display for EllipsoidSpace<R> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(write!(f, "Ellipsoid space with unit ball {}", self.comp[0]));
        for c in self.comp[1..].iter() {
            try!(write!(f, " + {}", c));
        }
        write!(f, " = 1")
    }
}
