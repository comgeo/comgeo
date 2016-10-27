#![deny(missing_docs,
   missing_debug_implementations, missing_copy_implementations,
   trivial_casts, trivial_numeric_casts,
   unsafe_code,
   unstable_features,
   unused_import_braces, unused_qualifications)]

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