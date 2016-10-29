
impl Default for Ostresh {
    fn default() -> Self {
        Ostresh::new()
    }
}
impl<P: Point> GeoMedianStep<P, EuclideanSpace> for Ostresh {
    type D = GeoMedianStepPrecisionErrorData;

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
            self.data.precisionerrors += 1;
            x.clone_from(node.p());
        }
    }

    fn data(&self) -> &Self::D {
        &self.data
    }

    fn print_data<W: Write>(&self, w: &mut W) -> io::Result<()> {
        writeln!(w, "{}", self.data)
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