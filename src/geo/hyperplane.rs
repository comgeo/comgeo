use std::slice::Iter;


pub trait Hyperplane<P: Point> {
    fn iter(&self) -> Iter<P::R>;
    fn comps(&self) -> &[P::R];
    fn comps_mut(&mut self) -> &mut [P::R];

    fn translate(&mut self, p: &P) {
        let d = self.iter().zip(p.iter())
            .fold(P::R::zero(), |sum, (&a, &x)|
                sum + a*x
            );

        for c in self.comps_mut().iter_mut() {
            *c = (*c)/d
        }
    }

    fn eucl_dist(&self, p: &P) -> P::R {
        let (dotm1, sqsum) = self.iter().zip(p.iter())
            .fold((-P::R::one(), P::R::zero()), |(dotm1, sqsum), (&a, &c)|
                (dotm1 + a*c, sqsum + a*a));

        dotm1.abs() / sqsum.sqrt()
    }

    fn normal<'a>(&self, p: &'a mut P) -> &'a mut P {
        for i in 0..p.coords().len() {
            p.coords_mut()[i] = self.comps()[i];
        }
        let norm = EuclideanSpace::new().norm(p);
        p.div(norm)
    }
}

#[derive(Debug)]
pub struct HyperplaneNd<R> {
    comp: Vec<R>
}
impl<R: Real> HyperplaneNd<R> {
    pub fn new(c: Vec<R>) -> Self {
        HyperplaneNd { comp: c }
    }
}
impl<R: Real, P: Point<R=R>> Hyperplane<P> for HyperplaneNd<R> {
    fn iter(&self) -> Iter<R> {
        self.comp[..].iter()
    }

    fn comps(&self) -> &[R] {
        &self.comp[..]
    }

    fn comps_mut(&mut self) -> &mut [R] {
        &mut self.comp[..]
    }
}
//impl<R: Real + fmt::Display> fmt::Display for HyperplaneNd<R> {
//    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//        try!(self.iter().fold(Ok(()), |_, n: &R| write!(f, "{} ", n)));
//        write!(f, "= 1")
//    }
//}

#[derive(Debug)]
pub struct Line<R> {
    comp: [R; 2]
}
impl<R: Real> Line<R> {
    pub fn new(c: [R; 2]) -> Self {
        Line { comp: c }
    }
}
impl<R: Real, P: Point<R=R>> Hyperplane<P> for Line<R> {
    fn iter(&self) -> Iter<R> {
        self.comp.iter()
    }

    fn comps(&self) -> &[R] {
        &self.comp
    }

    fn comps_mut(&mut self) -> &mut [R] {
        &mut self.comp
    }
}
//impl<E: Real> fmt::Display for Line<E> {
//    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//        print_hyperplane(self.iter(), f)
//    }
//}

#[derive(Debug)]
pub struct Plane<R> {
    comp: [R; 3]
}
impl<R: Real> Plane<R> {
    pub fn new(c: [R; 3]) -> Self {
        Plane { comp: c }
    }
}
impl<R: Real, P: Point<R=R>> Hyperplane<P> for Plane<R> {
    fn iter(&self) -> Iter<R> {
        self.comp.iter()
    }

    fn comps(&self) -> &[R] {
        &self.comp
    }

    fn comps_mut(&mut self) -> &mut [R] {
        &mut self.comp
    }
}
//impl<E: Real> fmt::Display for Plane<E> {
//    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//        print_hyperplane(self.iter(), f)
//    }
//}
