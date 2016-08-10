use traits::*;
use geo::spaces::{EuclideanSpace};
use algorithms::rmt::{GeoMedianIter};
use upperbounds::{LineTree};
use enumerator::*;
use steinertree::{SteinerTree};
use algorithms::geomedians::*;

use std::marker::PhantomData;
use std::time::{Duration, Instant};
use std::fmt;
use std::io::{self, BufWriter, Write};

#[derive(Debug)]
pub struct SteinerBnB<P, M, R, E, U> {
    rmt: R,
    enumerator: E,
    upperbound: U,
    data: SteinerBnBData,
    _m: PhantomData<M>,
    _p: PhantomData<P>
}

impl<P, M, K, E, U> SteinerBnB<P, M, K, E, U> {
    pub fn new(rmt: K, enumerator: E, u: U) -> Self {
        SteinerBnB {
            rmt: rmt,
            enumerator: enumerator,
            upperbound: u,
            data: SteinerBnBData::new(),
            _m: PhantomData,
            _p: PhantomData
        }
    }

    pub fn rmt_alg(&mut self) -> &mut K {
        &mut self.rmt
    }

    pub fn enumerator(&mut self) -> &mut E {
        &mut self.enumerator
    }


}

impl<P, M, K, E, U> SMT<P, M> for SteinerBnB<P, M, K, E, U>
    where P: Point, M: MinkowskiSpace<P>, K: RMT<P, M>, E: Enumerator<P>,
          U: UpperBound<P, M> {

    type D = SteinerBnBData;

    fn find(&mut self, t: Vec<P>, geo: &M) -> SteinerTree<P> {
        let start = Instant::now();
        let mut best = self.upperbound.bound(t.clone(), geo);
        let mut best_len = best.len(geo);
        self.enumerator.init(t, geo);

        while self.enumerator.next(geo) {
            if self.enumerator.tree().terminals().len() == best.terminals().len() {
                self.rmt.find(self.enumerator.tree_mut(), geo);
                let len = self.enumerator.tree().len(geo);
                if len < best_len {
                    self.data.best_updates += 1;
                    best = self.enumerator.tree().clone();
                    best_len = len;
                }
            }
        }

        //best.non_degenerate();
        self.data.time = Instant::now() - start;
        best
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        fn indent(f: &mut fmt::Formatter, indent: u32) -> fmt::Result {
            for _ in 0..indent {
                try!(write!(f, " "));
            }
            Ok(())
        }

        try!(write!(f, "Steiner branch and bound algorithm using the "));
        try!(self.enumerator.print(f, inde+4));
        try!({indent(f, inde);
            write!(f, ".\nRelatively minimal trees were found using the ");
            self.rmt.print(f, inde+4) });
        Ok(())
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        try!(writeln!(w, "{}", self.data()));
        try!(self.enumerator.print_data(w));
        self.rmt.print_data(w)
    }
}

impl<P, M, K, E, U> fmt::Display for SteinerBnB<P, M, K, E, U>
    where P: Point, M: MinkowskiSpace<P>, K: RMT<P, M>, E: Enumerator<P>,
          U: UpperBound<P, M> {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}

impl<P: Point> Default
    for SteinerBnB<P, EuclideanSpace, GeoMedianIter<P, EuclideanSpace, Uteshev>,
                   GPEnumeration<P, FurthestSiteOrdering>, LineTree> {

    fn default() -> Self {
        SteinerBnB::new(
            GeoMedianIter::default(),
            GPEnumeration::default(),
            LineTree::default())
    }
}

impl<P, H> SteinerBnB<P, H, GeoMedianIter<P, H, GeoMedianEllipsoid<P, Uteshev>>,
                      GPEnumeration<P, FurthestSiteOrdering>, LineTree>
    where P: Point, H: HyperEllipsoidSpace<P> {

    pub fn default_hyperellipsoid() -> Self {
        SteinerBnB::new(
            GeoMedianIter::default_with_geomedian(
                GeoMedianEllipsoid::default()),
            GPEnumeration::default(),
            LineTree::default())
    }
}

impl<P, M, R> SteinerBnB<P, M, R, GPEnumeration<P, FurthestSiteOrdering>, LineTree>
    where P: Point,  {

    pub fn default_rmt(rmt: R) -> Self {
        SteinerBnB::new(
            rmt,
            GPEnumeration::default(),
            LineTree::default())
    }
}

#[derive(Debug, Clone)]
pub struct SteinerBnBData {
    time: Duration,
    best_updates: u64
}

impl SteinerBnBData {
    fn new() -> Self {
        SteinerBnBData {
            time: Duration::new(0, 0),
            best_updates: 0
        }
    }
}

impl SmtData for SteinerBnBData {
    fn time(&self) -> &Duration {
        &self.time
    }

    fn best_updates(&self) -> u64  {
        self.best_updates
    }
}

impl fmt::Display for SteinerBnBData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Data for the Steiner branch and bound algorithm:"));
        try!(writeln!(f, "\tTotal time: {}", printdur(&self.time)));
        writeln!(f, "\tNumber of best updates: {}", self.best_updates)
    }
}
