//! This module implements arithmetic over the quadratic extension field Fp2.

use blst::*;

use core::{
    convert::TryInto,
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use fff::Field;

use crate::Fp;

#[derive(Copy, Clone)]
pub struct Fp2(pub(crate) blst_fp2);

impl fmt::Debug for Fp2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0(), self.c1())
    }
}

impl From<Fp> for Fp2 {
    fn from(f: Fp) -> Fp2 {
        Fp2::from_raw_unchecked(f, Fp::zero())
    }
}

impl From<blst_fp2> for Fp2 {
    fn from(val: blst_fp2) -> Fp2 {
        Fp2(val)
    }
}

impl From<Fp2> for blst_fp2 {
    fn from(val: Fp2) -> blst_fp2 {
        val.0
    }
}

impl Default for Fp2 {
    fn default() -> Self {
        Fp2::zero()
    }
}

impl Eq for Fp2 {}

impl PartialEq for Fp2 {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.fp[0].l == other.0.fp[0].l && self.0.fp[1].l == other.0.fp[1].l
    }
}

impl<'a> Neg for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn neg(self) -> Fp2 {
        self.neg()
    }
}

impl Neg for Fp2 {
    type Output = Fp2;

    #[inline]
    fn neg(self) -> Fp2 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn sub(self, rhs: &'b Fp2) -> Fp2 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn add(self, rhs: &'b Fp2) -> Fp2 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp2> for &'a Fp2 {
    type Output = Fp2;

    #[inline]
    fn mul(self, rhs: &'b Fp2) -> Fp2 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp2, Fp2);
impl_binops_multiplicative!(Fp2, Fp2);

impl Fp2 {
    /// Returns zero, the additive identity.
    #[inline]
    pub fn zero() -> Fp2 {
        Fp2(blst_fp2::default())
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub fn one() -> Fp2 {
        Fp2(blst_fp2 {
            fp: [Fp::one().into(), Fp::zero().into()],
        })
    }

    pub fn is_zero(&self) -> bool {
        self == &Fp2::zero()
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
        Fp2(blst_fp2 {
            fp: [Fp::random(rng).into(), Fp::random(rng).into()],
        })
    }

    /// Constructs an element of `Fp2` without checking that it is canonical.
    pub fn from_raw_unchecked(c0: Fp, c1: Fp) -> Fp2 {
        Fp2(blst_fp2 {
            fp: [c0.into(), c1.into()],
        })
    }

    #[inline]
    pub fn add(&self, rhs: &Fp2) -> Fp2 {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_add(&mut out, &self.0, &rhs.0) };

        Fp2(out)
    }

    #[inline]
    pub fn neg(&self) -> Fp2 {
        let mut out = blst_fp2::default();

        const FLAG: usize = 0x1;

        unsafe { blst_fp2_cneg(&mut out, &self.0, FLAG) };

        Fp2(out)
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp2) -> Fp2 {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_sub(&mut out, &self.0, &rhs.0) };

        Fp2(out)
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp2) -> Fp2 {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_mul(&mut out, &self.0, &rhs.0) };

        Fp2(out)
    }

    /// Squares this element.
    pub fn square(&self) -> Self {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_sqr(&mut out, &self.0) };

        Fp2(out)
    }

    /// Multiplies `self` with `3`, returning the result.
    pub fn mul3(&self) -> Self {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_mul_by_3(&mut out as _, &self.0 as _) };

        Fp2(out)
    }

    /// Multiplies `self` with `8`, returning the result.
    pub fn mul8(&self) -> Self {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_mul_by_8(&mut out as _, &self.0 as _) };

        Fp2(out)
    }

    /// Left shift `self` by `count`, returning the result.
    pub fn shl(&self, count: usize) -> Self {
        let mut out = blst_fp2::default();

        unsafe { blst_fp2_lshift(&mut out as _, &self.0 as _, count) };

        Fp2(out)
    }

    pub fn c0(&self) -> Fp {
        self.0.fp[0]
            .try_into()
            .expect("underlying fp must be valid")
    }

    pub fn c1(&self) -> Fp {
        self.0.fp[1]
            .try_into()
            .expect("underlying fp must be valid")
    }
}
