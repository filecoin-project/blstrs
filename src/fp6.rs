//! This module implements arithmetic over the quadratic extension field Fp6.

use blst::*;

use core::{
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{
    fp::{FROBENIUS_COEFF_FP6_C1, FROBENIUS_COEFF_FP6_C2},
    Fp, Fp2,
};

use fff::Field;

/// This represents an element $c_0 + c_1 v + c_2 v^2$ of $\mathbb{F}_{p^6} = \mathbb{F}_{p^2} / v^3 - u - 1$.
#[derive(Copy, Clone)]
pub struct Fp6(pub(crate) blst_fp6);

impl fmt::Debug for Fp6 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Fp6")
            .field("c0", &self.c0())
            .field("c1", &self.c1())
            .field("c2", &self.c2())
            .finish()
    }
}

impl ::std::fmt::Display for Fp6 {
    fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(
            f,
            "Fp6({} + {} * v, {} * v^2)",
            self.c0(),
            self.c1(),
            self.c2()
        )
    }
}

impl From<Fp> for Fp6 {
    fn from(f: Fp) -> Fp6 {
        Fp6::new(Fp2::from(f), Fp2::zero(), Fp2::zero())
    }
}

impl From<Fp2> for Fp6 {
    fn from(f: Fp2) -> Fp6 {
        Fp6::new(f, Fp2::zero(), Fp2::zero())
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
    /// Constructs an element of `Fp6`.
    pub const fn new(c0: Fp2, c1: Fp2, c2: Fp2) -> Fp6 {
        Fp6(blst_fp6 {
            fp2: [c0.0, c1.0, c2.0],
        })
    }

    #[inline]
    pub fn add(&self, rhs: &Fp6) -> Fp6 {
        let c0 = (self.c0() + rhs.c0()).0;
        let c1 = (self.c1() + rhs.c1()).0;
        let c2 = (self.c2() + rhs.c2()).0;

        Fp6(blst_fp6 { fp2: [c0, c1, c2] })
    }

    #[inline]
    pub fn neg(&self) -> Fp6 {
        let c0 = (-self.c0()).0;
        let c1 = (-self.c1()).0;
        let c2 = (-self.c2()).0;

        Fp6(blst_fp6 { fp2: [c0, c1, c2] })
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp6) -> Fp6 {
        let c0 = (self.c0() - rhs.c0()).0;
        let c1 = (self.c1() - rhs.c1()).0;
        let c2 = (self.c2() - rhs.c2()).0;

        Fp6(blst_fp6 { fp2: [c0, c1, c2] })
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp6) -> Fp6 {
        let mut a_a = self.c0();
        let mut b_b = self.c1();
        let mut c_c = self.c2();
        a_a *= &rhs.c0();
        b_b *= &rhs.c1();
        c_c *= &rhs.c2();

        let mut t1 = rhs.c1();
        t1 += &rhs.c2();
        {
            let mut tmp = self.c1();
            tmp += &self.c2();

            t1 *= &tmp;
            t1 -= &b_b;
            t1 -= &c_c;
            t1.mul_by_nonresidue();
            t1 += &a_a;
        }

        let mut t3 = rhs.c0();
        t3 += &rhs.c2();
        {
            let mut tmp = self.c0();
            tmp += &self.c2();

            t3 *= &tmp;
            t3 -= &a_a;
            t3 += &b_b;
            t3 -= &c_c;
        }

        let mut t2 = rhs.c0();
        t2 += &rhs.c1();
        {
            let mut tmp = self.c0();
            tmp += &self.c1();

            t2 *= &tmp;
            t2 -= &a_a;
            t2 -= &b_b;
            c_c.mul_by_nonresidue();
            t2 += &c_c;
        }

        Fp6::new(t1, t2, t3)
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

    /// Multiply by quadratic nonresidue v.
    pub fn mul_by_nonresidue(&mut self) {
        self.0.fp2.swap(0, 1);
        self.0.fp2.swap(0, 2);

        let mut c0 = self.c0();
        c0.mul_by_nonresidue();
        self.0.fp2[0] = c0.0;
    }
}

impl Field for Fp6 {
    fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        Fp6::new(Fp2::random(rng), Fp2::random(rng), Fp2::random(rng))
    }

    fn zero() -> Self {
        Fp6::new(Fp2::zero(), Fp2::zero(), Fp2::zero())
    }

    fn one() -> Self {
        Fp6::new(Fp2::one(), Fp2::zero(), Fp2::zero())
    }

    fn is_zero(&self) -> bool {
        self.c0().is_zero() && self.c1().is_zero() && self.c2().is_zero()
    }

    fn double(&mut self) {
        *self += *self;
    }

    fn negate(&mut self) {
        *self = self.neg();
    }

    fn add_assign(&mut self, other: &Self) {
        *self += other;
    }

    fn sub_assign(&mut self, other: &Self) {
        *self -= other;
    }

    fn frobenius_map(&mut self, power: usize) {
        let mut c0 = self.c0();
        c0.frobenius_map(power);
        let mut c1 = self.c1();
        c1.frobenius_map(power);
        let mut c2 = self.c2();
        c2.frobenius_map(power);

        c1 *= &FROBENIUS_COEFF_FP6_C1[power % 6];
        c2 *= &FROBENIUS_COEFF_FP6_C2[power % 6];

        self.0.fp2[0] = c0.0;
        self.0.fp2[1] = c1.0;
        self.0.fp2[2] = c2.0;
    }

    fn square(&mut self) {
        let mut s0 = self.c0();
        s0.square();
        let mut ab = self.c0();
        ab *= &self.c1();
        let mut s1 = ab;
        s1.double();
        let mut s2 = self.c0();
        s2 -= &self.c1();
        s2 += &self.c2();
        s2.square();
        let mut bc = self.c1();
        bc *= &self.c2();
        let mut s3 = bc;
        s3.double();
        let mut s4 = self.c2();
        s4.square();

        let mut c0 = s3;
        c0.mul_by_nonresidue();
        c0 += &s0;

        let mut c1 = s4;
        c1.mul_by_nonresidue();
        c1 += &s1;

        let mut c2 = s1;
        c2 += &s2;
        c2 += &s3;
        c2 -= &s0;
        c2 -= &s4;

        self.0.fp2[0] = c0.0;
        self.0.fp2[1] = c1.0;
        self.0.fp2[2] = c2.0;
    }

    fn mul_assign(&mut self, other: &Self) {
        *self *= other;
    }

    fn inverse(&self) -> Option<Self> {
        let mut c0 = self.c2();
        c0.mul_by_nonresidue();
        c0 *= &self.c1();
        c0 = -c0;
        {
            let mut c0s = self.c0();
            c0s.square();
            c0 += &c0s;
        }
        let mut c1 = self.c2();
        c1.square();
        c1.mul_by_nonresidue();
        {
            let mut c01 = self.c0();
            c01 *= &self.c1();
            c1 -= &c01;
        }
        let mut c2 = self.c1();
        c2.square();
        {
            let mut c02 = self.c0();
            c02 *= &self.c2();
            c2 -= &c02;
        }
        let mut tmp1 = self.c2();
        tmp1 *= &c1;
        let mut tmp2 = self.c1();
        tmp2 *= &c2;
        tmp1 += &tmp2;
        tmp1.mul_by_nonresidue();
        tmp2 = self.c0();
        tmp2 *= &c0;
        tmp1 += &tmp2;
        match tmp1.inverse() {
            Some(t) => Some(Fp6::new(t * c0, t * c1, t * c2)),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Fp2, Fp6};

    use fff::{Field, PrimeField};
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_fp6_mul_nonresidue() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let nqr = Fp6::new(Fp2::zero(), Fp2::one(), Fp2::zero());

        for _ in 0..1000 {
            let mut a = Fp6::random(&mut rng);
            let mut b = a;
            a.mul_by_nonresidue();
            b.mul_assign(&nqr);

            assert_eq!(a, b);
        }
    }

    #[test]
    fn fp6_random_field_tests() {
        crate::tests::field::random_field_tests::<Fp6>();
    }

    #[test]
    fn fp6_random_frobenius_tests() {
        crate::tests::field::random_frobenius_tests::<Fp6, _>(crate::Fp::char(), 13);
    }
}
