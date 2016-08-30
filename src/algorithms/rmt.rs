use std::marker::PhantomData;
use std::time::{Duration, Instant};
use std::fmt;
use std::io::{self, BufWriter, Write};

use traits::*;
use geo::spaces::{EuclideanSpace};
use steinertree::{SteinerTree};
use algorithms::geomedians::*;

#[derive(Debug, Clone)]
pub struct GeoMedianIterData {
    nodes: usize,
    time: Duration,
    iterations: u64,
    selftime: Duration
}

impl GeoMedianIterData {
    pub fn iterations(&self) -> u64 {
        self.iterations
    }

    pub fn selftime(&self) -> &Duration {
        &self.selftime
    }
}

impl RmtData for GeoMedianIterData {
    fn time(&self) -> &Duration {
        &self.time
    }

    fn nodes(&self) -> usize  {
        self.nodes
    }
}

impl fmt::Display for GeoMedianIterData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Data for the geometric median iteration RMT algorithm:"));
        try!(writeln!(f, "\tNumber of topologies optimized: {}", self.nodes));
        try!(writeln!(f, "\tTotal number iterations {} (one iteration optimizes \
            all Steiner points of a topology)", self.iterations));
        try!(writeln!(f, "\tAvarage number of iterations pr. topology: {}",
            (self.iterations as f64) / (self.nodes as f64)));
        try!(writeln!(f, "\tTotal time: {}", printdur(&self.time)));
        writeln!(f, "\tAvarage time pr. topology: {}", printdur(&(self.time / (self.nodes as u32))))
    }
}

#[derive(Debug)]
pub struct GeoMedianIter<P: Point, M: MinkowskiSpace<P>, G: GeoMedian<P, M>> {
    tree_len_cutoff: P::R,
    median: G,
    data: GeoMedianIterData,
    _m: PhantomData<M>
}

impl<P, M, G> GeoMedianIter<P, M, G>
    where P: Point, M: MinkowskiSpace<P>, G: GeoMedian<P, M> {

    pub fn new(a: P::R, median: G) -> GeoMedianIter<P, M, G> {
        GeoMedianIter {
            tree_len_cutoff: a,
            median: median,
            data: GeoMedianIterData {
                nodes: 0,
                time: Duration::new(0, 0),
                selftime: Duration::new(0, 0),
                iterations: 0
            },
            _m: PhantomData
        }
    }

    pub fn default_with_geomedian(median: G) -> GeoMedianIter<P, M, G> {
        Self::new(P::R::from(0.00001), median)
    }

    pub fn tree_len_cutoff(mut self, a: P::R) -> Self {
        self.tree_len_cutoff = a;
        self
    }

    pub fn geo_median_alg(&mut self) -> &mut G {
        &mut self.median
    }
}

impl<P> Default for GeoMedianIter<P, EuclideanSpace, Uteshev>
    where P: Point {

    fn default() -> Self {
        Self::new(P::R::from(0.00001), Uteshev::default())
    }
}

impl<P, M, G> RMT<P, M> for GeoMedianIter<P, M, G>
    where P: Point, M: MinkowskiSpace<P>, G: GeoMedian<P, M> {

    type D = GeoMedianIterData;

    fn find(&mut self, stree: &mut SteinerTree<P>, geo: &M) -> P::R {
        self.data.nodes += 1;
        let start = Instant::now();

        let mut last_len = stree.len(geo);
        for s in stree.steiner_points() {
            self.median.init(s, geo);
        }
        loop {
            self.data.iterations += 1;
            for s in stree.steiner_points() {
                self.median.find(s, geo);
            }

            let len = stree.len(geo);
            if last_len - len < self.tree_len_cutoff {
                self.data.time += Instant::now() - start;
                self.data.selftime = self.data.time - *self.median.data().time();
                return len;
            }

            last_len = len;
        }
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

        try!(write!(f, "Geometric median iterator that continuously finds \
            geometric medians until the change in tree length is less than {}.",
            self.tree_len_cutoff));
        try!(write!(f, " Geometric medians were found using "));
        self.median.print(f, inde)
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        try!(writeln!(w, "{}", self.data()));
        self.median.print_data(w)
    }
}

impl<P, M, G> fmt::Display for GeoMedianIter<P, M, G>
    where P: Point, M: MinkowskiSpace<P>, G: GeoMedian<P, M> {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}

