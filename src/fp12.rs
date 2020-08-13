//! This module implements arithmetic over the quadratic extension field Fp12.

use blst::*;

use core::{
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{Fp, Fp2, Fp6};

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
#[derive(Copy, Clone)]
pub struct Fp12(pub(crate) blst_fp12);

impl fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0(), self.c1())
    }
}

impl From<Fp> for Fp12 {
    fn from(f: Fp) -> Fp12 {
        Fp12::from_raw_unchecked(Fp6::from(f), Fp6::zero())
    }
}

impl From<Fp2> for Fp12 {
    fn from(f: Fp2) -> Fp12 {
        Fp12::from_raw_unchecked(Fp6::from(f), Fp6::zero())
    }
}

impl From<Fp6> for Fp12 {
    fn from(f: Fp6) -> Fp12 {
        Fp12::from_raw_unchecked(f, Fp6::zero())
    }
}

impl From<blst_fp12> for Fp12 {
    fn from(val: blst_fp12) -> Fp12 {
        Fp12(val)
    }
}

impl From<Fp12> for blst_fp12 {
    fn from(val: Fp12) -> blst_fp12 {
        val.0
    }
}

impl Default for Fp12 {
    fn default() -> Self {
        Fp12::zero()
    }
}

impl Eq for Fp12 {}

impl PartialEq for Fp12 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.fp6[0].fp2[0].fp[0].l == other.0.fp6[0].fp2[0].fp[0].l
            && self.0.fp6[0].fp2[0].fp[1].l == other.0.fp6[0].fp2[0].fp[1].l
            && self.0.fp6[0].fp2[1].fp[0].l == other.0.fp6[0].fp2[1].fp[0].l
            && self.0.fp6[0].fp2[1].fp[1].l == other.0.fp6[0].fp2[1].fp[1].l
            && self.0.fp6[0].fp2[2].fp[0].l == other.0.fp6[0].fp2[2].fp[0].l
            && self.0.fp6[0].fp2[2].fp[1].l == other.0.fp6[0].fp2[2].fp[1].l
            && self.0.fp6[1].fp2[0].fp[0].l == other.0.fp6[1].fp2[0].fp[0].l
            && self.0.fp6[1].fp2[0].fp[1].l == other.0.fp6[1].fp2[0].fp[1].l
            && self.0.fp6[1].fp2[1].fp[0].l == other.0.fp6[1].fp2[1].fp[0].l
            && self.0.fp6[1].fp2[1].fp[1].l == other.0.fp6[1].fp2[1].fp[1].l
            && self.0.fp6[1].fp2[2].fp[0].l == other.0.fp6[1].fp2[2].fp[0].l
            && self.0.fp6[1].fp2[2].fp[1].l == other.0.fp6[1].fp2[2].fp[1].l
    }
}

impl<'a> Neg for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        self.neg()
    }
}

impl Neg for Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn sub(self, rhs: &'b Fp12) -> Fp12 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn add(self, rhs: &'b Fp12) -> Fp12 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, rhs: &'b Fp12) -> Fp12 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp12, Fp12);
impl_binops_multiplicative!(Fp12, Fp12);

impl Fp12 {
    /// Returns zero, the additive identity.
    #[inline]
    pub fn zero() -> Fp12 {
        Fp12(blst_fp12::default())
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub fn one() -> Fp12 {
        Fp12(blst_fp12 {
            fp6: [Fp6::one().into(), Fp6::zero().into()],
        })
    }

    pub fn is_zero(&self) -> bool {
        self == &Fp12::zero()
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
        Fp12(blst_fp12 {
            fp6: [Fp6::random(rng).into(), Fp6::random(rng).into()],
        })
    }

    /// Constructs an element of `Fp12` without checking that it is canonical.
    pub fn from_raw_unchecked(c0: Fp6, c1: Fp6) -> Fp12 {
        Fp12(blst_fp12 {
            fp6: [c0.into(), c1.into()],
        })
    }

    #[inline]
    pub fn add(&self, rhs: &Fp12) -> Fp12 {
        todo!()
    }

    #[inline]
    pub fn neg(&self) -> Fp12 {
        todo!()
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp12) -> Fp12 {
        todo!()
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp12) -> Fp12 {
        let mut out = blst_fp12::default();

        unsafe { blst_fp12_mul(&mut out, &self.0, &rhs.0) };

        Fp12(out)
    }

    /// Squares this element.
    pub fn square(&self) -> Self {
        let mut out = blst_fp12::default();

        unsafe { blst_fp12_sqr(&mut out, &self.0) };

        Fp12(out)
    }

    pub fn c0(&self) -> Fp6 {
        self.0.fp6[0].into()
    }

    pub fn c1(&self) -> Fp6 {
        self.0.fp6[1].into()
    }
}
