#[derive(Debug)]
pub struct BrimbergLove {
    data: GeoMedianStepFixedPointData
}
impl BrimbergLove {
    pub fn new() -> Self {
        BrimbergLove {
            data: GeoMedianStepFixedPointData::new()
        }
    }
}
impl Default for BrimbergLove {
    fn default() -> Self {
        BrimbergLove::new()
    }
}
impl<P: Point> GeoMedianStep<P, LpSpace<P::R>> for BrimbergLove {
    type D = GeoMedianStepFixedPointData;

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
                self.data.fixedpoints += 1;
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
                    self.data.fixedpoints += 1;
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

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data)
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