#[derive(Debug)]
pub struct Weiszfeld {
    data: GeoMedianStepFixedPointData
}
impl Weiszfeld {
    pub fn new() -> Self {
        Weiszfeld {
            data: GeoMedianStepFixedPointData::new()
        }
    }
}
impl Default for Weiszfeld {
    fn default() -> Self {
        Weiszfeld::new()
    }
}
impl<P: Point> GeoMedianStep<P, EuclideanSpace> for Weiszfeld {
    type D = GeoMedianStepFixedPointData;

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
            self.data.fixedpoints += 1;
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
        write!(f, "Weiszfeld's iteration")
    }
}


impl fmt::Display for Weiszfeld {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Weiszfeld's iteration")
    }
}

#[derive(Debug)]
pub struct Ostresh {
    data: GeoMedianStepPrecisionErrorData
}
impl Ostresh {
    pub fn new() -> Self {
        Ostresh {
            data: GeoMedianStepPrecisionErrorData::new()
        }
    }
}