pub mod brimberglove;
pub mod chiafranco;
pub mod ostresh;
pub mod weiszfeld;

use std::marker::PhantomData;
use std::time::{Duration, Instant};
use std::fmt;
use std::io::{self, BufWriter, Write};

use geo::spaces::*;
use steinertree::{Node};


pub trait GeoMedianStep<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: GeoMedianStepData;

    fn step(&mut self, &mut P, &mut Node<P>, usize, &M);
    fn init(&mut self, &mut Node<P>, &M);
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}


pub trait GeoMedianStepData: fmt::Display + Clone {
}

fn centroid<P: Point>(y: &Node<P>) -> P {
    let mut centroid = y.p().clone();
    centroid.mul(P::R::zero());

    for n in y.neighbours() {
        centroid.add(n.p());
    }

    centroid.div(P::R::from(y.neighbours().len() as f64));
    centroid
}



#[derive(Debug, Clone)]
pub struct GeoMedianStepPrecisionErrorData {
    precisionerrors: u64
}

impl GeoMedianStepPrecisionErrorData {
    pub fn new() -> Self {
        GeoMedianStepPrecisionErrorData {
            precisionerrors: 0
        }
    }

    pub fn precision_errors(&self) -> u64 {
        self.precisionerrors
    }
}

impl GeoMedianStepData for GeoMedianStepPrecisionErrorData { }

impl fmt::Display for GeoMedianStepPrecisionErrorData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "Data for the geometric median step function:"));
        writeln!(f, "\tTotal number of precision errors encountered: {}", self.precisionerrors)
    }
}


#[derive(Debug, Clone)]
pub struct GeoMedianStepFixedPointData {
    fixedpoints: u64
}

impl GeoMedianStepFixedPointData {
    pub fn new() -> Self {
        GeoMedianStepFixedPointData {
            fixedpoints: 0
        }
    }

    pub fn fixed_points(&self) -> u64 {
        self.fixedpoints
    }
}

impl GeoMedianStepData for GeoMedianStepFixedPointData { }

impl fmt::Display for GeoMedianStepFixedPointData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "Data for the geometric median step function:"));
        writeln!(f, "\tTotal number of fixed points encountered: {}", self.fixedpoints)
    }
}













#[derive(Debug)]
pub struct GeoMedianStepper<P: Point, I> {
    node_dist_cutoff: P::R,
    step: I,
    data: GeoMedianStepperData
}

impl<P: Point, I> GeoMedianStepper<P, I> {
    pub fn new(a: P::R, step: I) -> Self {
        GeoMedianStepper {
            node_dist_cutoff: a,
            step: step,
            data: GeoMedianStepperData::new()
        }
    }

    pub fn default_with_step(step: I) -> Self {
        Self::new(P::R::from(0.00001), step)
    }

    pub fn step_alg(&mut self) -> &mut I {
        &mut self.step
    }
}

impl<P: Point> Default for GeoMedianStepper<P, Ostresh> {
    fn default() -> Self {
        GeoMedianStepper::new(P::R::from(0.0001), Ostresh::default())
    }
}

impl<P, M, I> GeoMedian<P, M> for GeoMedianStepper<P, I>
    where P: Point, M: MinkowskiSpace<P>, I: GeoMedianStep<P, M> {

    type D = GeoMedianStepperData;

    fn init(&mut self, y: &mut Node<P>, geo: &M) {
        self.data.inits += 1;
        self.step.init(y, geo);
    }

    fn find(&mut self, node: &mut Node<P>, geo: &M) {
        let start = Instant::now();
        let mut x = node.p().clone();
        let mut s: usize = 1;

        loop {
            self.data.total_steps += 1;
            self.step.step(&mut x, node, s, geo);
            let change = geo.dist(&x, node.p());

            node.p_mut().clone_from(&x);

            debug_assert!(change.is_number());
            if change < self.node_dist_cutoff {
                self.data.time += Instant::now() - start;
                return;
            }

            s += 1;
        }
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        try!(self.step.print(f, inde));
        write!(f, " stopping when the change in position \
            (measured in the current space) got below {}", self.node_dist_cutoff)
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        try!(writeln!(w, "{}", self.data()));
        self.step.print_data(w)
    }
}

impl<P, I> fmt::Display for GeoMedianStepper<P, I>
    where P: Point  {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Geometric median stepper")
    }
}

#[derive(Debug, Clone)]
pub struct GeoMedianStepperData {
    time: Duration,
    total_steps: u64,
    inits: u64
}

impl GeoMedianStepperData {
    pub fn new() -> Self {
        GeoMedianStepperData {
            time: Duration::new(0, 0),
            total_steps: 0,
            inits: 0
        }
    }

    pub fn total_steps(&self) -> u64 {
        self.total_steps
    }

    pub fn inits(&self) -> u64 {
        self.inits
    }

    pub fn average_steps(&self) -> f64 {
        (self.total_steps as f64) / (self.inits as f64)
    }

    pub fn average_time_step(&self) -> f64 {
        let time = (self.time.as_secs() as f64) + (self.time.subsec_nanos() as f64) / 1000000000.0;
        time / (self.total_steps as f64)
    }

    pub fn average_time_problem(&self) -> f64 {
        self.average_time_step() * self.average_steps()
    }
}

impl GeoMedianData for GeoMedianStepperData {
    fn time(&self) -> &Duration {
        &self.time
    }
}

impl fmt::Display for GeoMedianStepperData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Data for the geometric median iteration:"));
        writeln!(f, "\tTotal time: {}", printdur(&self.time));
        writeln!(f, "\tTotal steps: {}", self.total_steps);
        writeln!(f, "\tAverage number of steps pr. geometric median: {}", self.average_steps());
        writeln!(f, "\tAverage time pr. step: {}", self.average_time_step());
        writeln!(f, "\tAverage time pr. geometric median: {}", self.average_time_problem())
    }
}




















