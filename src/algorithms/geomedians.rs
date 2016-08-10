use std::marker::PhantomData;
use std::time::{Duration, Instant};
use std::fmt;
use std::io::{self, BufWriter, Write};

use traits::*;
use geo::spaces::*;
use steinertree::{Node};


fn centroid<P: Point>(y: &Node<P>) -> P {
    let mut centroid = y.p().clone();
    centroid.mul(P::R::zero());

    for n in y.neighbours() {
        centroid.add(n.p());
    }

    centroid.div(P::R::from(y.neighbours().len() as f64));
    centroid
}



#[derive(Debug)]
pub struct ChiaFranco<P: Point, F> {
    e: F,
    _m: PhantomData<P>
}
impl<P: Point, F: Fn(usize) -> P::R> ChiaFranco<P, F> {
    pub fn new(e: F) -> Self {
        ChiaFranco {
            e: e,
            _m: PhantomData
        }
    }

    pub fn get_epsilon(&self, s: usize) -> P::R {
        (self.e)(s)
    }

    pub fn get_err_ub(&self, tlen: usize, p: P::R, s: usize) -> P::R {
        P::R::from(tlen as f64) * P::R::from(2.0).pow(P::R::one()/p) * self.get_epsilon(s).sqrt()
    }
}

impl<P: Point, F: Fn(usize) -> P::R> GeoMedianStep<P, LpSpace<P::R>> for ChiaFranco<P, F> {
    fn init(&mut self, _: &mut Node<P>, _: &LpSpace<P::R>) {
    }

    fn step(&mut self, ynext: &mut P, y: &mut Node<P>, s: usize, geo: &LpSpace<P::R>) {
        debug_assert!(P::R::from(2.0) < geo.p());

        #[inline]
        fn h<P: Point>(point: &mut P, epsilon: P::R) {
            point.scale(&|c| hk(c, epsilon));
        }

        #[inline]
        fn hk<R: Real>(c: R, epsilon: R) -> R {
            (c * c + epsilon).sqrt()
        }

        #[inline]
        fn yk<P: Point, F: Fn(usize) -> P::R>(cf: &ChiaFranco<P, F>, y: &Node<P>, k: usize,
                        lambda: P::R, ptmp: &mut P, s: usize, geo: &LpSpace<P::R>) -> P::R {
            let mut numerator = P::R::zero();
            let mut denominator = P::R::zero();
            let yk = y.p().coords()[k];

            for t in y.neighbours() {
                let tk = t.p().coords()[k];

                ptmp.clone_from(y.p());
                h(ptmp.sub(t.p()), cf.get_epsilon(s));
                let normterm = geo.norm(ptmp).pow(geo.p()-P::R::one());
                let common = hk(yk-tk, cf.get_epsilon(s)).pow(geo.p() - P::R::from(2.0)) / normterm;

                denominator += common;
                numerator += common * ((P::R::one() - lambda) * yk + lambda * tk);
            }

            numerator / denominator
        }

        let mut tmp = y.p().clone();
        let p = geo.p();
        let (one, two) = (P::R::one(), P::R::from(2.0));

        // Calculate Gamma
        let mut gamma = y.p().clone();
        let mut gamma_k = 0;
        gamma.scale_mut(&mut |_| {
            let new = yk(self, y, gamma_k, one, &mut tmp, s, geo);
            gamma_k += 1;
            new
        });

        // Calculate psi(0) and psi(f_k)
        let (psi0, psif) = {
            let (mut psi0, mut psif) = (P::R::zero(), P::R::zero());

            for (k, (&yk, &gk)) in y.p().iter().zip(gamma.iter()).enumerate() {
                let mut square = gk - yk;
                square *= square;

                for t in y.neighbours() {
                    let tk = t.p().coords()[k];

                    tmp.clone_from(y.p());
                    h(tmp.sub(t.p()), self.get_epsilon(s));
                    let normterm = geo.norm(&tmp).pow(p-one);

                    psi0 += square * (hk(yk-tk, self.get_epsilon(s)) / normterm);
                    psif += square * (hk(yk-tk+(two / (p-one))*(gk-yk), self.get_epsilon(s)) / normterm);
                }
            }
            (psi0, psif)
        };

        let lambda = if p > P::R::from(3.0) {
            let denom = -p*(p-one)*psif;
            (two / (p-one)).min((-two*psi0*p) / denom)
        } else {
            let mut sumterm = P::R::zero();
            for t in y.neighbours() {
                tmp.clone_from(y.p());
                h(tmp.sub(t.p()), self.get_epsilon(s));
                sumterm += geo.norm(&tmp).pow(one-p);
            }
            tmp.clone_from(&gamma);
            let normterm = geo.norm(tmp.sub(y.p())).pow(p);
            let denom = p*psi0 - p*p*psi0 - p*(p-one)*(two/(p-one)).pow(p-two) * normterm * sumterm;
            (-two*psi0*p) / denom
        };

        let mut ynext_k = 0;
        ynext.scale_mut(&mut |_| {
            let new = yk(self, y, ynext_k, lambda, &mut tmp, s, geo);
            ynext_k += 1;
            new
        });
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Rodríguez-chía and Valero-Franco's iteration")
    }
}

impl<P: Point, F: Fn(usize) -> P::R> fmt::Display for ChiaFranco<P, F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}

#[derive(Debug)]
pub struct BrimbergLove { }
impl BrimbergLove {
    pub fn new() -> Self {
        BrimbergLove { }
    }
}
impl Default for BrimbergLove {
    fn default() -> Self {
        BrimbergLove::new()
    }
}
impl<P: Point> GeoMedianStep<P, LpSpace<P::R>> for BrimbergLove {
    fn init(&mut self, y: &mut Node<P>, geo: &LpSpace<P::R>) {
        if y.neighbours().any(|n| geo.dist(n.p(), y.p()) == P::R::zero()) {
            let cent = centroid(y);
            y.p_mut().clone_from(&cent);
        }
    }

    fn step(&mut self, ynext: &mut P, y: &mut Node<P>, _: usize, geo: &LpSpace<P::R>) {
        debug_assert!(P::R::one() <= geo.p() && geo.p() <= P::R::from(2.0));

        ynext.mul(P::R::zero());
        let mut div = ynext.clone();
        let mut sings = ynext.clone();
        let mut singsval = ynext.clone();

        for t in y.neighbours() {
            let dpow = geo.dist(t.p(), y.p()).pow(P::R::one() - geo.p());
            if !dpow.is_number() {
                ynext.clone_from(t.p());
                return;
            }
            for (singsk, (singskval, (ynk, (divk, (&yk, &tk)))))
                in sings.iter_mut().zip(
                   singsval.iter_mut().zip(
                   ynext.iter_mut().zip(
                   div.iter_mut().zip(
                   y.p().iter().zip(
                   t.p().iter() ))))) {

                let a = (yk-tk).abs().pow(geo.p()-P::R::from(2.0)) * dpow;
                if !a.is_number() {
                    *singsk = P::R::one();
                    *singskval = tk;
                    continue;
                }
                *ynk += tk * a;
                *divk += a;
            }
        }

        ynext.modify(&div, &|c, oc| c / oc);
        for (&singsk, (&singskval, yk))
            in sings.iter().zip(
               singsval.iter().zip(
               ynext.iter_mut() )) {

            if singsk == P::R::one() {
                *yk = singskval;
            }
        }
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Brimberg and Love's iteration")
    }
}


impl fmt::Display for BrimbergLove {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Brimberg and Love's iteration")
    }
}

#[derive(Debug)]
pub struct Weiszfeld { }
impl Weiszfeld {
    pub fn new() -> Self {
        Weiszfeld { }
    }
}
impl Default for Weiszfeld {
    fn default() -> Self {
        Weiszfeld::new()
    }
}
impl<P: Point> GeoMedianStep<P, EuclideanSpace> for Weiszfeld {
    fn init(&mut self, y: &mut Node<P>, geo: &EuclideanSpace) {
        if y.neighbours().any(|n| geo.dist(n.p(), y.p()) == P::R::zero()) {
            let cent = centroid(y);
            y.p_mut().clone_from(&cent);
        }
    }

    fn step(&mut self, x: &mut P, node: &mut Node<P>, _: usize, geo: &EuclideanSpace) {
        let mut div_sum = P::R::zero();
        x.mul(P::R::zero());

        for n in node.neighbours() {
            let d = geo.dist(n.p(), node.p());
            x.modify(n.p(), &|c, nc| c + nc / d);
            div_sum = div_sum + d.recip();
        }

        x.div(div_sum);
        if !x.iter().any(|c| !c.is_number()) {
            x.clone_from(node.p());
        }
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Weiszfeld's iteration")
    }
}


impl fmt::Display for Weiszfeld {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Weiszfeld's iteration")
    }
}

#[derive(Debug)]
pub struct Ostresh { }
impl Ostresh {
    pub fn new() -> Self {
        Ostresh { }
    }
}
impl Default for Ostresh {
    fn default() -> Self {
        Ostresh::new()
    }
}
impl<P: Point> GeoMedianStep<P, EuclideanSpace> for Ostresh {
    fn init(&mut self, _: &mut Node<P>, _: &EuclideanSpace) {
    }

    fn step(&mut self, x: &mut P, node: &mut Node<P>, _: usize, geo: &EuclideanSpace) {
        x.mul(P::R::zero());
        let mut s = P::R::zero();
        let (mut g, mut tmp) = (x.clone(), x.clone());
        let mut singularity = false;

        for n in node.neighbours() {
            let d = geo.dist(n.p(), node.p());
            tmp.clone_from(node.p());
            let gadd = tmp.sub(n.p()).div(d);
            if gadd.iter().any(|c| !c.is_number()) || !(s + d.recip()).is_number() {
                singularity = true;
                continue;
            }
            g.add(gadd);
            s += d.recip();
        }

        if singularity {
            if geo.norm(&g) <= P::R::one() {
                g.mul(P::R::zero());
            } else {
                tmp.clone_from(&g);
                g.sub(tmp.unit(geo));
            }
        }

        x.clone_from(node.p());
        x.sub(g.div(s));
        if x.iter().any(|c| !c.is_number()) {
            x.clone_from(node.p());
        }
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Weiszfeld's iteration with Ostresh's modification")
    }
}


impl fmt::Display for Ostresh {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Weiszfeld's iteration with Ostresh's modification")
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
        writeln!(w, "{}", self.data())
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
        writeln!(w, "{}", self.data)
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











#[derive(Debug)]
pub struct Uteshev {
    data: UteshevData,
}

impl Uteshev {
    pub fn new() -> Self {
        Uteshev {
            data: UteshevData::new(),
        }
    }
}

impl Default for Uteshev {
    fn default() -> Self {
        Uteshev::new()
    }
}

impl<P: Point> GeoMedian<P, EuclideanSpace> for Uteshev {
    type D = UteshevData;

    fn init(&mut self, y: &mut Node<P>, geo: &EuclideanSpace) {
        if y.neighbours().any(|n| geo.dist(n.p(), y.p()) == P::R::zero()) {
            let cent = centroid(y);
            y.p_mut().clone_from(&cent);
        }
    }

    fn find(&mut self, node: &mut Node<P>, geo: &EuclideanSpace) {
        let start = Instant::now();

        #[inline]
        fn s<R: Real>(xy: R, xz: R, yz: R) -> R {
            R::from(0.5) * (
                (xy + xz + yz) * (xy + xz) * (xy + yz) * (xz + yz)
            ).sqrt()
        }

        #[inline]
        fn ca<R: Real>(ab2: R, ac2: R, bc2: R, s: R) -> R {
            (ab2 + ac2 - bc2)
                * (R::from(3.0).sqrt() / R::from(2.0))
                + s
        }

        let (mut p, mut i) = node.neighbours_data_mut();
        let (x, y, z) = (i.next().unwrap(), i.next().unwrap(), i.next().unwrap());
        let (xy, xz, yz) = (geo.dist(x, y), geo.dist(x, z), geo.dist(y, z));
        let (xy2, xz2, yz2) = (xy*xy, xz*xz, yz*yz);
        let cosx = (xy2 + xz2 - yz2) / (P::R::from(2.0) * xy * xz);
        let cosy = (xy2 + yz2 - xz2) / (P::R::from(2.0) * xy * yz);
        let cosz = (xz2 + yz2 - xy2) / (P::R::from(2.0) * xz * yz);

        if cosx <= P::R::from(-0.5) {
            p.clone_from(x);
        } else if cosy <= P::R::from(-0.5) {
            p.clone_from(y);
        } else if cosz <= P::R::from(-0.5) {
            p.clone_from(z);
        } else {
            let s = s(xy, xz, yz);
            let (cxr, cyr, czr) = (
                ca(xy2, xz2, yz2, s).recip(),
                ca(xy2, yz2, xz2, s).recip(),
                ca(xz2, yz2, xy2, s).recip()
            );

            for (pk, (&xk, (&yk, &zk))) in p.iter_mut().zip(x.iter().zip(y.iter().zip(z.iter()))) {
                let a = xk*cxr + yk*cyr + zk*czr;
                if a.is_number() {
                    *pk = a;
                }
            }

            let mul = (cxr+cyr+czr).recip();
            p.scale_mut(&mut |c| {
                let a = mul*c;
                if a.is_number() {
                    a
                } else {
                    c
                }
            });
        }

        self.data.time += Instant::now() - start;
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Analytical solution by Uteshev")
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data)
    }
}


impl fmt::Display for Uteshev {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Analytical solution by Uteshev")
    }
}


#[derive(Debug, Clone)]
pub struct UteshevData {
    time: Duration
}

impl UteshevData {
    pub fn new() -> Self {
        UteshevData {
            time: Duration::new(0, 0)
        }
    }
}

impl GeoMedianData for UteshevData {
    fn time(&self) -> &Duration {
        &self.time
    }
}

impl fmt::Display for UteshevData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        fn printdur(dur: &Duration) -> f64 {
            (dur.as_secs() as f64) + (dur.subsec_nanos() as f64) / 1000000000.0
        }

        try!(writeln!(f, "Geometric median data for Uteshev's analytical solution:"));
        writeln!(f, "\tTotal time: {}", printdur(&self.time))
    }
}
