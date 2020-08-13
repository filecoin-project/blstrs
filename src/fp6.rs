//! This module implements arithmetic over the quadratic extension field Fp6.

use blst::*;

use core::{
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{Fp, Fp2};

/// This represents an element $c_0 + c_1 v + c_2 v^2$ of $\mathbb{F}_{p^6} = \mathbb{F}_{p^2} / v^3 - u - 1$.
#[derive(Copy, Clone)]
pub struct Fp6(pub(crate) blst_fp6);

impl fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:?} + ({:?})*v + ({:?})*v^2",
            self.c0(),
            self.c1(),
            self.c2()
        )
    }
}

impl From<Fp> for Fp6 {
    fn from(f: Fp) -> Fp6 {
        Fp6::from_raw_unchecked(Fp2::from(f), Fp2::zero(), Fp2::zero())
    }
}

impl From<Fp2> for Fp6 {
    fn from(f: Fp2) -> Fp6 {
        Fp6::from_raw_unchecked(f, Fp2::zero(), Fp2::zero())
    }
}

impl From<blst_fp6> for Fp6 {
    fn from(val: blst_fp6) -> Fp6 {
        Fp6(val)
    }
}

impl From<Fp6> for blst_fp6 {
    fn from(val: Fp6) -> blst_fp6 {
        val.0
    }
}

impl Default for Fp6 {
    fn default() -> Self {
        Fp6::zero()
    }
}

impl Eq for Fp6 {}

impl PartialEq for Fp6 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.fp2[0].fp[0].l == other.0.fp2[0].fp[0].l
            && self.0.fp2[0].fp[1].l == other.0.fp2[0].fp[1].l
            && self.0.fp2[1].fp[0].l == other.0.fp2[1].fp[0].l
            && self.0.fp2[1].fp[1].l == other.0.fp2[1].fp[1].l
            && self.0.fp2[2].fp[0].l == other.0.fp2[2].fp[0].l
            && self.0.fp2[2].fp[1].l == other.0.fp2[2].fp[1].l
    }
}

impl<'a> Neg for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn neg(self) -> Fp6 {
        self.neg()
    }
}

impl Neg for Fp6 {
    type Output = Fp6;

    #[inline]
    fn neg(self) -> Fp6 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn sub(self, rhs: &'b Fp6) -> Fp6 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn add(self, rhs: &'b Fp6) -> Fp6 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp6> for &'a Fp6 {
    type Output = Fp6;

    #[inline]
    fn mul(self, rhs: &'b Fp6) -> Fp6 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp6, Fp6);
impl_binops_multiplicative!(Fp6, Fp6);

impl Fp6 {
    /// Returns zero, the additive identity.
    #[inline]
    pub fn zero() -> Fp6 {
        Fp6(blst_fp6::default())
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub fn one() -> Fp6 {
        Fp6(blst_fp6 {
            fp2: [Fp2::one().into(), Fp2::zero().into(), Fp2::zero().into()],
        })
    }

    pub fn is_zero(&self) -> bool {
        self == &Fp6::zero()
    }

    /// Exponentiates this element by a number represented with `u64` limbs,
    /// least significant digit first.
    pub fn pow(&self, by: &[u64; 6]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();

                if ((*e >> i) & 1) == 1 {
                    res *= self;
                }
            }
        }
        res
    }

    /// Computes a uniformly random element using rejection sampling.
    pub fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        Fp6(blst_fp6 {
            fp2: [
                Fp2::random(rng).into(),
                Fp2::random(rng).into(),
                Fp2::random(rng).into(),
            ],
        })
    }

    /// Constructs an element of `Fp6` without checking that it is canonical.
    pub fn from_raw_unchecked(c0: Fp2, c1: Fp2, c2: Fp2) -> Fp6 {
        Fp6(blst_fp6 {
            fp2: [c0.into(), c1.into(), c2.into()],
        })
    }

    #[inline]
    pub fn add(&self, rhs: &Fp6) -> Fp6 {
        todo!()
    }

    #[inline]
    pub fn neg(&self) -> Fp6 {
        todo!()
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp6) -> Fp6 {
        todo!()
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp6) -> Fp6 {
        todo!()
    }

    /// Squares this element.
    pub fn square(&self) -> Self {
        todo!()
    }

    pub fn c0(&self) -> Fp2 {
        self.0.fp2[0].into()
    }

    pub fn c1(&self) -> Fp2 {
        self.0.fp2[1].into()
    }

    pub fn c2(&self) -> Fp2 {
        self.0.fp2[2].into()
    }
}
