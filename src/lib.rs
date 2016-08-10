pub mod traits;
pub mod geo;
pub mod algorithms;
pub mod upperbounds;
pub mod enumerator;
pub mod prunetests;
pub mod steinertree;

use traits::*;
use geo::points::*;

impl Real for f64 {
    fn abs(self) -> Self {
        self.abs()
    }

    fn recip(self) -> Self {
        self.recip()
    }

    fn zero() -> Self {
        0f64
    }

    fn one() -> Self {
        1f64
    }

    fn max(self, o: Self) -> Self {
        self.max(o)
    }

    fn min(self, o: Self) -> Self {
        self.min(o)
    }

    fn sqrt(self) -> Self {
        self.sqrt()
    }

    fn pow(self, exp: Self) -> Self {
        self.powf(exp)
    }

    fn is_number(&self) -> bool {
        self.is_finite()
    }
}


#[test]
fn it_works() {
    let mut tree = SteinerTree::new(
        &[ // Terminals
            Point2d::new([0.0f64, 0.0]), // 0
            Point2d::new([1.0f64, 0.0]), // 1
            Point2d::new([0.0f64, 1.0]), // 2
            Point2d::new([1.0f64, 1.0]), // 3
        ], &[ // Steiner points
            Point2d::new([0.2f64, 0.15]), // 4
        ], &[ // Edges
            (0,4), (1,4), (2,4), (3,4)
        ]);

    //let s = EllipseSpace::new(1f64, 1f64);
    let s = EuclideanSpace;

    println!("steiner_points().len() = {}", tree.steiner_points().len());

    let rmt = GeoMedianIter::default();

    rmt.find(&mut tree, &s);

    print!("{}", tree);
    assert_eq!(2.1, tree.len(&s));
}

#[test]
fn test_enumeration() {

    let mut smt = SMT::default();

    let tree = smt.find(vec![
        Point2d::new([0.0f64, 0.0]), // 0
        Point2d::new([1.0f64, 0.0]), // 1
        //Point2d::new([0.0f64, 1.0]), // 2
        //Point2d::new([1.0f64, 1.0]), // 3
    ], &EuclideanSpace);

    println!("Best tree: {}", tree);
    assert_eq!(2.1, tree.len(&EuclideanSpace));
}
