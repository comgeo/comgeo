pub mod uteshev;
pub mod iterative;


pub trait GeoMedian<P: Point, M: MinkowskiSpace<P>>: fmt::Display {
    type D: GeoMedianData;

    fn find(&mut self, &mut Node<P>, &M);
    fn init(&mut self, &mut Node<P>, &M);
    fn data(&self) -> &Self::D;
    fn print(&self, &mut fmt::Formatter, u32) -> fmt::Result;
    fn print_data<W: Write>(&self, &mut W) -> io::Result<()>;
}


pub trait GeoMedianData: fmt::Display + Clone {
    fn time(&self) -> &Duration;
}




#[derive(Debug)]
pub struct GeoMedianEllipsoid<P, E>
    where P: Point, E: GeoMedian<P, EuclideanSpace> {

    eucl_median: E,
    data: GeoMedianEllipsoidData,
    _m: PhantomData<P>
}

impl<P: Point, E: GeoMedian<P, EuclideanSpace>> GeoMedianEllipsoid<P, E> {
    pub fn new(eucl_median: E) -> Self {
        GeoMedianEllipsoid {
            eucl_median: eucl_median,
            data: GeoMedianEllipsoidData::new(),
            _m: PhantomData
        }
    }

    pub fn geo_median_alg(&mut self) -> &mut E {
        &mut self.eucl_median
    }
}

impl<P: Point> Default for GeoMedianEllipsoid<P, Uteshev> {
    fn default() -> Self {
        GeoMedianEllipsoid::new(Uteshev::default())
    }
}

impl<P, S, E> GeoMedian<P, S> for GeoMedianEllipsoid<P, E>
where P: Point, S: HyperEllipsoidSpace<P>,
      E: GeoMedian<P, EuclideanSpace> {

    type D = GeoMedianEllipsoidData;

    fn init(&mut self, y: &mut Node<P>, _: &S) {
        self.eucl_median.init(y, &EuclideanSpace);
    }

    fn find(&mut self, node: &mut Node<P>, geo: &S) {
        let start = Instant::now();
        for n in node.neighbours_mut() {
            for (c, d) in n.p_mut().iter_mut().zip(geo.comps().iter()) {
                *c = *c / *d;
            }
        }

        self.eucl_median.find(node, &EuclideanSpace);

        for n in node.neighbours_mut() {
            for (c, d) in n.p_mut().iter_mut().zip(geo.comps().iter()) {
                *c = *c * *d;
            }
        }
        for (c, d) in node.p_mut().iter_mut().zip(geo.comps().iter()) {
            *c = *c * *d;
        }
        self.data.time += Instant::now() - start;
        self.data.selftime = self.data.time - *self.eucl_median.data().time();
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

        try!(writeln!(f, "Ellipsoid geometric median finder that transforms the input then use the "));
        self.eucl_median.print(f, inde)
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        try!(writeln!(w, "{}", self.data));
        self.eucl_median.print_data(w)
    }
}


impl<P, E> fmt::Display for GeoMedianEllipsoid<P, E>
where P: Point, E: GeoMedian<P, EuclideanSpace> {

    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        try!(writeln!(f, "Ellipsoid geometric median finder that transforms the input then use the "));
        self.eucl_median.print(f, 0)
    }
}


#[derive(Debug, Clone)]
pub struct GeoMedianEllipsoidData {
    time: Duration,
    selftime: Duration
}

impl GeoMedianEllipsoidData {
    pub fn new() -> Self {
        GeoMedianEllipsoidData {
            time: Duration::new(0, 0),
            selftime: Duration::new(0, 0)
        }
    }
}

impl GeoMedianData for GeoMedianEllipsoidData {
    fn time(&self) -> &Duration {
        &self.time
    }
}

impl fmt::Display for GeoMedianEllipsoidData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Geometric median ellipsoid data:"));
        writeln!(f, "\tTotal time: {}", printdur(&self.time));
        writeln!(f, "\tTotal time on coordinate transformations: {}", printdur(&self.selftime))
    }
}