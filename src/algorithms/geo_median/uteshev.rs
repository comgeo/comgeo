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