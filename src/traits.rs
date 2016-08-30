use std::fmt;
use std::slice::{Iter, IterMut};
use std::ops::*;
use steinertree::{SteinerTree, Node};
use std::time::{Duration};
use std::io::{self, Write};

use geo::spaces::{EuclideanSpace};

pub trait Point : PartialEq + Clone + fmt::Display + fmt::Debug {
    type R : Real;

    fn dim(&self) -> usize;
    fn iter(&self) -> Iter<Self::R>;
    fn iter_mut(&mut self) -> IterMut<Self::R>;
    fn coords(&self) -> &[Self::R];
    fn coords_mut(&mut self) -> &mut [Self::R];
    fn id(&self) -> usize;
    fn set_id(&mut self, usize);

    fn dot(&self, o: &Self) -> Self::R {
        debug_assert_eq!(self.dim(), o.dim());
        self.iter().zip(o.iter()).fold(Self::R::zero(), |dot, (&c, &oc)| dot + c*oc)
    }

    fn is_zero(&self) -> bool {
        for &c in self.coords() {
            if c != Self::R::zero() {
                return false;
            }
        }
        true
    }

    fn modify<'a, F>(&mut self, p: &'a Self, f: &F) -> &mut Self
        where F: Fn(Self::R, Self::R) -> Self::R {

        debug_assert_eq!(self.dim(), p.dim());
        for (c, &oc) in self.coords_mut().iter_mut().zip(p.coords().iter()) {
            *c = f(*c, oc)
        }
        self
    }

    fn scale<'a, F>(&mut self, f: &F) -> &mut Self
        where F: Fn(Self::R) -> Self::R {

        for c in self.coords_mut().iter_mut() {
            *c = f(*c)
        }
        self
    }

    fn scale_mut<'a, F>(&mut self, f: &mut F) -> &mut Self
        where F: FnMut(Self::R) -> Self::R {

        for c in self.coords_mut().iter_mut() {
            *c = (*f)(*c)
        }
        self
    }

    fn add<'a>(&mut self, rhs: &'a Self) -> &mut Self {
        self.modify(rhs, &|c, oc| c + oc)
    }

    fn sub<'a>(&mut self, rhs: &'a Self) -> &mut Self {
        self.modify(rhs, &|c, oc| c - oc)
    }

    fn mul(&mut self, rhs: Self::R) -> &mut Self {
        self.scale(&|c| c * rhs)
    }

    fn div(&mut self, rhs: Self::R) -> &mut Self {
        self.scale(&|c| c / rhs)
    }

    fn unit<M: MinkowskiSpace<Self>>(&mut self, geo: &M) -> &mut Self {
        let norm = geo.norm(self);
        self.div(norm)
    }

    fn neg(&mut self) -> &mut Self {
        self.scale(&|c| -c)
    }

    fn copy(&mut self, o: &Self) -> &mut Self {
        self.modify(o, &|_, oc| oc)
    }
}

pub trait Real : Neg<Output=Self> + Add<Output=Self>
                + Sub<Output=Self> + Mul<Output=Self>
                + Div<Output=Self> + Rem<Output=Self>
                + AddAssign + SubAssign + MulAssign + DivAssign
                + From<f64> + Into<f64>
                + Copy + PartialEq + PartialOrd
                + fmt::Display + fmt::Debug {

    fn abs(self) -> Self;
    fn recip(self) -> Self;
    fn zero() -> Self;
    fn one() -> Self;
    fn max(self, o: Self) -> Self;
    fn min(self, o: Self) -> Self;
    fn sqrt(self) -> Self;
    fn pow(self, Self) -> Self;
    fn is_number(&self) -> bool;
}

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

pub trait MinkowskiSpace<P: Point>: fmt::Display {
    fn norm(&self, &P) -> P::R;
    fn dist(&self, p1: &P, p2: &P) -> P::R {
        self.norm(p1.clone().sub(p2))
    }
}

pub trait HyperEllipsoidSpace<P: Point> : MinkowskiSpace<P> {
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

pub trait PruneTest {
    fn prunetest<P: Point>(&self, &SteinerTree<P>, P::R, P::R) -> bool;
}

pub trait UpperBound<P: Point, M: MinkowskiSpace<P>> {
    fn bound(&self, Vec<P>, geo: &M) -> SteinerTree<P>;
}

pub trait Enumerator<P: Point>: fmt::Display {
    type D: EnumeratorData;

    fn init<M: MinkowskiSpace<P>>(&mut self, Vec<P>, &M);
    fn next<M: MinkowskiSpace<P>>(&mut self, &M) -> bool;
    fn tree(&self) -> &SteinerTree<P>;
    fn tree_mut(&mut self) -> &mut SteinerTree<P>;
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}

pub trait TerminalSorter: fmt::Display {
    fn sort<P, M>(&mut self, &mut[P], &M)
        where P: Point, M: MinkowskiSpace<P>;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
}

pub trait RMT<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: RmtData;

    fn find(&mut self, tree: &mut SteinerTree<P>, geo: &M) -> P::R;
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}

pub trait SMT<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: SmtData;

    fn find(&mut self, Vec<P>, geo: &M) -> SteinerTree<P>;
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}

pub trait MST<P: Point, M: MinkowskiSpace<P>> {
    fn find(&mut self, &[P], geo: &M) -> SteinerTree<P>;
}

pub trait GeoMedian<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: GeoMedianData;

    fn find(&mut self, &mut Node<P>, &M);
    fn init(&mut self, &mut Node<P>, &M);
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}

pub trait GeoMedianStep<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: GeoMedianStepData;

    fn step(&mut self, &mut P, &mut Node<P>, usize, &M);
    fn init(&mut self, &mut Node<P>, &M);
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}

pub trait RmtData: fmt::Display + Clone  {
    fn nodes(&self) -> usize;
    fn time(&self) -> &Duration;
}

pub trait SmtData: fmt::Display + Clone {
    fn time(&self) -> &Duration;
    fn best_updates(&self) -> u64;
}

pub trait GeoMedianData: fmt::Display + Clone {
    fn time(&self) -> &Duration;
}

pub trait GeoMedianStepData: fmt::Display + Clone {
}

pub trait EnumeratorData: fmt::Display + Clone {
    fn nodes(&self) -> usize;
    fn pruned(&self) -> usize;
    fn time(&self) -> &Duration;
}
