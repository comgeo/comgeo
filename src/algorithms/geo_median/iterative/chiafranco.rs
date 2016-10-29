#[derive(Debug)]
pub struct ChiaFrancoApprox<P: Point> {
    e: P::R,
    data: GeoMedianStepPrecisionErrorData
}
impl<P: Point> ChiaFrancoApprox<P> {
    pub fn new(e: P::R) -> Self {
        ChiaFrancoApprox {
            e: e,
            data: GeoMedianStepPrecisionErrorData::new()
        }
    }

    pub fn get_epsilon(&self) -> P::R {
        self.e
    }

    pub fn set_epsilon(&mut self, e: P::R) {
        self.e = e;
    }

    pub fn get_err_ub(&self, tlen: usize, p: P::R) -> P::R {
        P::R::from(tlen as f64) * P::R::from(2.0).pow(P::R::one()/p) * self.e.sqrt()
    }
}

impl<P: Point> GeoMedianStep<P, LpSpace<P::R>> for ChiaFrancoApprox<P> {
    type D = GeoMedianStepPrecisionErrorData;

    fn init(&mut self, _: &mut Node<P>, _: &LpSpace<P::R>) {
    }

    fn step(&mut self, ynext: &mut P, y: &mut Node<P>, _: usize, geo: &LpSpace<P::R>) {
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
        fn yk<P: Point>(cf: &ChiaFrancoApprox<P>, y: &Node<P>, k: usize,
                        lambda: P::R, ptmp: &mut P, geo: &LpSpace<P::R>) -> P::R {
            let mut numerator = P::R::zero();
            let mut denominator = P::R::zero();
            let yk = y.p().coords()[k];

            for t in y.neighbours() {
                let tk = t.p().coords()[k];

                ptmp.clone_from(y.p());
                h(ptmp.sub(t.p()), cf.get_epsilon());
                let normterm = geo.norm(ptmp).pow(geo.p()-P::R::one());
                let common = hk(yk-tk, cf.get_epsilon()).pow(geo.p() - P::R::from(2.0)) / normterm;

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
            let new = yk(self, y, gamma_k, one, &mut tmp, geo);
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
                    h(tmp.sub(t.p()), self.get_epsilon());
                    let normterm = geo.norm(&tmp).pow(p-one);

                    psi0 += square * (hk(yk-tk, self.get_epsilon()) / normterm);
                    psif += square * (hk(yk-tk+(two / (p-one))*(gk-yk), self.get_epsilon()) / normterm);
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
                h(tmp.sub(t.p()), self.get_epsilon());
                sumterm += geo.norm(&tmp).pow(one-p);
            }
            tmp.clone_from(&gamma);
            let normterm = geo.norm(tmp.sub(y.p())).pow(p);
            let denom = p*psi0 - p*p*psi0 - p*(p-one)*(two/(p-one)).pow(p-two) * normterm * sumterm;
            (-two*psi0*p) / denom
        };

        let mut ynext_k = 0;
        ynext.scale_mut(&mut |_| {
            let new = yk(self, y, ynext_k, lambda, &mut tmp, geo);
            ynext_k += 1;
            new
        });
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print(&self, f: &mut fmt::Formatter, inde: u32) -> fmt::Result {
        write!(f, "Rodríguez-chía and Valero-Franco's iteration with the \
            hyberbolic approximation where epsilon={}", self.e)
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data)
    }
}

impl<P: Point> fmt::Display for ChiaFrancoApprox<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.print(f, 0)
    }
}



#[derive(Debug)]
pub struct ChiaFranco<P: Point, F> {
    e: F,
    step: ChiaFrancoApprox<P>,
}
impl<P: Point, F: Fn(usize) -> P::R> ChiaFranco<P, F> {
    pub fn new(e: F) -> Self {
        ChiaFranco {
            step: ChiaFrancoApprox::new(e(1)),
            e: e
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
    type D = GeoMedianStepPrecisionErrorData;

    fn init(&mut self, _: &mut Node<P>, _: &LpSpace<P::R>) {
    }

    fn step(&mut self, ynext: &mut P, y: &mut Node<P>, s: usize, geo: &LpSpace<P::R>) {
        self.step.set_epsilon((self.e)(s));
        self.step.step(ynext, y, s, geo);
    }

    fn data(&self) -> &Self::D {
        &self.step.data()
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data())
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