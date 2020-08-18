//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`

use blst::*;

use core::{
    convert::TryInto,
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use fff::{Field, PrimeField};

/// `Fp` values are always in
/// Montgomery form; i.e., Scalar(a) = aR mod p, with R = 2^384.
#[derive(Copy, Clone)]
pub struct Fp(pub(crate) blst_fp);

/// Representation of a `Fp`, in regular coordinates.
#[derive(Default, Clone, Copy)]
pub struct FpRepr(blst_fp);

// -((2**384) mod q) mod q
pub(crate) const NEGATIVE_ONE: Fp = Fp(blst_fp {
    l: [
        0x43f5fffffffcaaae,
        0x32b7fff2ed47fffd,
        0x7e83a49a2e99d69,
        0xeca8f3318332bb7a,
        0xef148d1ea0f4c069,
        0x40ab3263eff0206,
    ],
});

// Coefficients for the Frobenius automorphism.
pub(crate) const FROBENIUS_COEFF_FP2_C1: [Fp; 2] = [
    // Fp(-1)**(((q^0) - 1) / 2)
    Fp(blst_fp {
        l: [
            0x760900000002fffd,
            0xebf4000bc40c0002,
            0x5f48985753c758ba,
            0x77ce585370525745,
            0x5c071a97a256ec6d,
            0x15f65ec3fa80e493,
        ],
    }),
    // Fp(-1)**(((q^1) - 1) / 2)
    Fp(blst_fp {
        l: [
            0x43f5fffffffcaaae,
            0x32b7fff2ed47fffd,
            0x7e83a49a2e99d69,
            0xeca8f3318332bb7a,
            0xef148d1ea0f4c069,
            0x40ab3263eff0206,
        ],
    }),
];

impl AsRef<[u64]> for FpRepr {
    fn as_ref(&self) -> &[u64] {
        &self.0.l
    }
}

impl AsMut<[u64]> for FpRepr {
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0.l
    }
}

const LIMBS: usize = 6;
const LIMB_BITS: usize = 64;

impl fmt::Debug for FpRepr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x")?;
        for &b in self.0.l.iter().rev() {
            write!(f, "{:016x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for FpRepr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x")?;
        for &b in self.0.l.iter().rev() {
            write!(f, "{:016x}", b)?;
        }
        Ok(())
    }
}

impl From<u64> for FpRepr {
    fn from(val: u64) -> FpRepr {
        FpRepr(blst_fp {
            l: [val, 0, 0, 0, 0, 0],
        })
    }
}

impl From<u64> for Fp {
    fn from(val: u64) -> Fp {
        let mut out = blst_fp::default();
        unsafe { blst_fp_from_uint64(&mut out, &val) };
        Fp(out)
    }
}

impl Ord for FpRepr {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        for (a, b) in self.0.l.iter().rev().zip(other.0.l.iter().rev()) {
            if a < b {
                return std::cmp::Ordering::Less;
            } else if a > b {
                return std::cmp::Ordering::Greater;
            }
        }

        std::cmp::Ordering::Equal
    }
}

impl PartialOrd for FpRepr {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for FpRepr {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l == other.0.l
    }
}
impl Eq for FpRepr {}

impl fff::PrimeFieldRepr for FpRepr {
    fn sub_noborrow(&mut self, other: &Self) {
        let mut borrow = 0;

        for (a, b) in self.0.l.iter_mut().zip(other.0.l.iter()) {
            *a = fff::sbb(*a, *b, &mut borrow);
        }
    }

    fn add_nocarry(&mut self, other: &Self) {
        let mut carry = 0;

        for (a, b) in self.0.l.iter_mut().zip(other.0.l.iter()) {
            *a = fff::adc(*a, *b, &mut carry);
        }
    }

    fn num_bits(&self) -> u32 {
        let mut ret = (LIMBS as u32) * LIMB_BITS as u32;
        for i in self.0.l.iter().rev() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != LIMB_BITS as u32 {
                break;
            }
        }

        ret
    }

    fn is_zero(&self) -> bool {
        self.0.l.iter().all(|&e| e == 0)
    }

    fn is_odd(&self) -> bool {
        self.0.l[0] & 1 == 1
    }

    fn is_even(&self) -> bool {
        !self.is_odd()
    }

    fn div2(&mut self) {
        let mut t = 0;
        for i in self.0.l.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    fn shr(&mut self, mut n: u32) {
        if n as usize >= LIMB_BITS * LIMBS {
            *self = Self::from(0);
            return;
        }

        while n >= LIMB_BITS as u32 {
            let mut t = 0;
            for i in self.0.l.iter_mut().rev() {
                std::mem::swap(&mut t, i);
            }
            n -= LIMB_BITS as u32;
        }

        if n > 0 {
            let mut t = 0;
            for i in self.0.l.iter_mut().rev() {
                let t2 = *i << (LIMB_BITS as u32 - n);
                *i >>= n;
                *i |= t;
                t = t2;
            }
        }
    }

    fn mul2(&mut self) {
        let mut last = 0;
        for i in &mut self.0.l {
            let tmp = *i >> 63;
            *i <<= 1;
            *i |= last;
            last = tmp;
        }
    }

    fn shl(&mut self, mut n: u32) {
        if n as usize >= LIMB_BITS * LIMBS {
            *self = Self::from(0);
            return;
        }

        while n >= LIMB_BITS as u32 {
            let mut t = 0;
            for i in &mut self.0.l {
                std::mem::swap(&mut t, i);
            }
            n -= LIMB_BITS as u32;
        }

        if n > 0 {
            let mut t = 0;
            for i in &mut self.0.l {
                let t2 = *i >> (LIMB_BITS as u32 - n);
                *i <<= n;
                *i |= t;
                t = t2;
            }
        }
    }
}

/// Elements are ordered lexicographically.
impl Ord for Fp {
    #[inline(always)]
    fn cmp(&self, other: &Fp) -> ::std::cmp::Ordering {
        self.into_repr().cmp(&other.into_repr())
    }
}

impl PartialOrd for Fp {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fp) -> Option<::std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes_le();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes_le();
        write!(f, "Fp(0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

impl From<Fp> for FpRepr {
    fn from(val: Fp) -> Self {
        let mut out = blst_fp::default();
        unsafe { blst_fp_from(&mut out, &val.0) };
        FpRepr(out)
    }
}

impl From<Fp> for blst_fp {
    fn from(val: Fp) -> blst_fp {
        val.0
    }
}

impl From<blst_fp> for Fp {
    fn from(val: blst_fp) -> Fp {
        Fp(val)
    }
}

impl Default for Fp {
    fn default() -> Self {
        Fp::zero()
    }
}

impl Eq for Fp {}

impl PartialEq for Fp {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l[0] == other.0.l[0]
            && self.0.l[1] == other.0.l[1]
            && self.0.l[2] == other.0.l[2]
            && self.0.l[3] == other.0.l[3]
            && self.0.l[4] == other.0.l[4]
            && self.0.l[5] == other.0.l[5]
    }
}

impl<'a> Neg for &'a Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        self.neg()
    }
}

impl Neg for Fp {
    type Output = Fp;

    #[inline]
    fn neg(self) -> Fp {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn sub(self, rhs: &'b Fp) -> Fp {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn add(self, rhs: &'b Fp) -> Fp {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp> for &'a Fp {
    type Output = Fp;

    #[inline]
    fn mul(self, rhs: &'b Fp) -> Fp {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp, Fp);
impl_binops_multiplicative!(Fp, Fp);

impl fff::Field for Fp {
    fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        // The number of bits we should "shave" from a randomly sampled reputation.
        const REPR_SHAVE_BITS: usize = 384 - Fp::NUM_BITS as usize;

        loop {
            let mut raw = blst_fp::default();
            for i in 0..4 {
                raw.l[i] = rng.next_u64();
            }

            // Mask away the unused most-significant bits.
            raw.l[3] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if let Ok(valid_el) = raw.try_into() {
                return valid_el;
            }
        }
    }

    fn zero() -> Self {
        Fp(blst_fp {
            l: [0, 0, 0, 0, 0, 0],
        })
    }

    fn one() -> Self {
        Fp(blst_fp {
            l: [
                0x760900000002fffd,
                0xebf4000bc40c0002,
                0x5f48985753c758ba,
                0x77ce585370525745,
                0x5c071a97a256ec6d,
                0x15f65ec3fa80e493,
            ],
        })
    }

    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    fn square(&mut self) {
        let mut raw = blst_fp::default();
        unsafe { blst_fp_sqr(&mut raw, &self.0) }

        self.0 = raw;
    }

    fn double(&mut self) {
        *self += *self;
    }

    fn negate(&mut self) {
        *self = -&*self;
    }
    fn add_assign(&mut self, other: &Self) {
        *self += other;
    }

    fn sub_assign(&mut self, other: &Self) {
        *self -= other;
    }

    fn mul_assign(&mut self, other: &Self) {
        *self = self.mul(other)
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let mut out = blst_fp::default();

        unsafe { blst_fp_eucl_inverse(&mut out, &self.0) };

        Some(Fp(out))
    }

    fn frobenius_map(&mut self, _: usize) {
        // This has no effect in a prime field.
    }
}

const MODULUS: FpRepr = FpRepr(blst_fp {
    l: [
        0xb9feffffffffaaab,
        0x1eabfffeb153ffff,
        0x6730d2a0f6b0f624,
        0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,
        0x1a0111ea397fe69a,
    ],
});

impl FpRepr {
    pub fn new(raw: [u64; 6]) -> Self {
        FpRepr(blst_fp { l: raw })
    }
}

impl fff::PrimeField for Fp {
    type Repr = FpRepr;

    const NUM_BITS: u32 = 381;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const S: u32 = 1;

    fn from_repr(repr: Self::Repr) -> Result<Self, fff::PrimeFieldDecodingError> {
        if FpRepr(repr.0) < MODULUS {
            let mut out = blst_fp::default();
            unsafe { blst_fp_to(&mut out, &repr.0) }
            Ok(Fp(out))
        } else {
            Err(fff::PrimeFieldDecodingError::NotInField(
                "not in field".to_string(),
            ))
        }
    }

    /// Convert a biginteger representation into a prime field element, if
    /// the number is an element of the field.
    fn into_repr(&self) -> Self::Repr {
        (*self).into()
    }

    fn char() -> Self::Repr {
        MODULUS
    }

    fn multiplicative_generator() -> Self {
        Fp::from_repr(FpRepr(blst_fp {
            l: [2, 0, 0, 0, 0, 0],
        }))
        .unwrap()
    }

    fn root_of_unity() -> Self {
        Fp::from_repr(FpRepr(blst_fp {
            l: [
                13402431016077863594,
                2210141511517208575,
                7435674573564081700,
                7239337960414712511,
                5412103778470702295,
                1873798617647539866,
            ],
        }))
        .unwrap()
    }
}

impl fff::SqrtField for Fp {
    fn legendre(&self) -> fff::LegendreSymbol {
        const MOD_MINUS_1_OVER_2: [u64; 6] = [
            15924587544893707605,
            1105070755758604287,
            12941209323636816658,
            12843041017062132063,
            2706051889235351147,
            936899308823769933,
        ];
        // s = self^((modulus - 1) // 2)
        let s = self.pow(MOD_MINUS_1_OVER_2);
        if s == Self::zero() {
            ::fff::LegendreSymbol::Zero
        } else if s == Self::one() {
            ::fff::LegendreSymbol::QuadraticResidue
        } else {
            ::fff::LegendreSymbol::QuadraticNonResidue
        }
    }

    fn sqrt(&self) -> Option<Self> {
        // Shank's algorithm for q mod 4 = 3
        // https://eprint.iacr.org/2012/685.pdf (page 9, algorithm 2)

        let mut a1 = self.pow(&[
            17185665809301629610u64,
            552535377879302143u64,
            15693976698673184137u64,
            15644892545385841839u64,
            10576397981472451381u64,
            468449654411884966u64,
        ]);

        let mut a0 = a1;
        a0.square();
        a0 *= self;

        const RNEG: [u64; 6] = [
            4897101644811774638u64,
            3654671041462534141u64,
            569769440802610537u64,
            17053147383018470266u64,
            17227549637287919721u64,
            291242102765847046u64,
        ];

        if a0.0.l == RNEG {
            None
        } else {
            a1 *= self;
            Some(a1)
        }
    }
}

impl Fp {
    /// Attempts to convert a little-endian byte representation of
    /// a scalar into an `Fp`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 48]) -> Option<Fp> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_fp::default();

        unsafe {
            blst_fp_from_lendian(&mut raw, in_v.as_mut_ptr());
        }

        raw.try_into().ok()
    }

    /// Attempts to convert a big-endian byte representation of
    /// a scalar into an `Fp`, failing if the input is not canonical.
    pub fn from_bytes_be(bytes: &[u8; 48]) -> Option<Fp> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_fp::default();

        unsafe {
            blst_fp_from_bendian(&mut raw, in_v.as_mut_ptr());
        }

        raw.try_into().ok()
    }

    /// Converts an element of `Fp` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 48] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 48];

        unsafe {
            blst_lendian_from_fp(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 48];
        out.copy_from_slice(&out_v);

        out
    }

    /// Converts an element of `Fp` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes_be(&self) -> [u8; 48] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 48];

        unsafe {
            blst_bendian_from_fp(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 48];
        out.copy_from_slice(&out_v);

        out
    }

    /// Constructs an element of `Fp` without checking that it is canonical.
    pub fn from_raw_unchecked(v: [u64; 6]) -> Fp {
        let mut inner = blst_fp::default();
        inner.l.copy_from_slice(&v);
        Fp(inner)
    }

    #[inline]
    pub fn add(&self, rhs: &Fp) -> Fp {
        let mut out = blst_fp::default();

        unsafe { blst_fp_add(&mut out, &self.0, &rhs.0) };

        Fp(out)
    }

    #[inline]
    pub fn neg(&self) -> Fp {
        let mut out = blst_fp::default();

        const FLAG: usize = 0x1;

        unsafe { blst_fp_cneg(&mut out, &self.0, FLAG) };

        Fp(out)
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp) -> Fp {
        let mut out = blst_fp::default();

        unsafe { blst_fp_sub(&mut out, &self.0, &rhs.0) };

        Fp(out)
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp) -> Fp {
        let mut out = blst_fp::default();

        unsafe { blst_fp_mul(&mut out, &self.0, &rhs.0) };

        Fp(out)
    }

    /// Multiplies `self` with `3`, returning the result.
    pub fn mul3(&self) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_mul_by_3(&mut out, &self.0) };

        Fp(out)
    }

    /// Multiplies `self` with `8`, returning the result.
    pub fn mul8(&self) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_mul_by_8(&mut out, &self.0) };

        Fp(out)
    }

    /// Left shift `self` by `count`, returning the result.
    pub fn shl(&self, count: usize) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_lshift(&mut out, &self.0, count) };

        Fp(out)
    }
}

#[cfg(test)]
mod tests {
    use super::{Fp, FpRepr};

    use fff::{Field, PrimeField, PrimeFieldRepr};
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_modulus() {
        assert_eq!(
            format!("{:?}", Fp::char()), "0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab"
            );
    }

    #[test]
    fn test_neg_one() {
        let mut o = Fp::one();
        o.negate();

        assert_eq!(
            Fp(FpRepr::new([
                0x43f5fffffffcaaae,
                0x32b7fff2ed47fffd,
                0x7e83a49a2e99d69,
                0xeca8f3318332bb7a,
                0xef148d1ea0f4c069,
                0x40ab3263eff0206,
            ])
            .0),
            o
        );
    }

    #[test]
    fn test_fp_repr_ordering() {
        use std::cmp::Ordering;

        fn assert_equality(a: FpRepr, b: FpRepr) {
            assert_eq!(a, b);
            assert!(a.cmp(&b) == Ordering::Equal);
        }

        fn assert_lt(a: FpRepr, b: FpRepr) {
            assert!(a < b);
            assert!(b > a);
        }

        assert_equality(
            FpRepr::new([9999, 9999, 9999, 9999, 9999, 9999]),
            FpRepr::new([9999, 9999, 9999, 9999, 9999, 9999]),
        );
        assert_equality(
            FpRepr::new([9999, 9998, 9999, 9999, 9999, 9999]),
            FpRepr::new([9999, 9998, 9999, 9999, 9999, 9999]),
        );
        assert_equality(
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9999]),
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9999]),
        );
        assert_lt(
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9998]),
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9999]),
        );
        assert_lt(
            FpRepr::new([9999, 9999, 9999, 9997, 9998, 9999]),
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9999]),
        );
        assert_lt(
            FpRepr::new([9, 9999, 9999, 9997, 9998, 9999]),
            FpRepr::new([9999, 9999, 9999, 9997, 9999, 9999]),
        );
    }

    #[test]
    fn test_fp_repr_from() {
        assert_eq!(FpRepr::from(100), FpRepr::new([100, 0, 0, 0, 0, 0]));
    }

    #[test]
    fn test_fp_repr_is_odd() {
        assert!(!FpRepr::from(0).is_odd());
        assert!(FpRepr::from(0).is_even());
        assert!(FpRepr::from(1).is_odd());
        assert!(!FpRepr::from(1).is_even());
        assert!(!FpRepr::from(324834872).is_odd());
        assert!(FpRepr::from(324834872).is_even());
        assert!(FpRepr::from(324834873).is_odd());
        assert!(!FpRepr::from(324834873).is_even());
    }

    #[test]
    fn test_fp_repr_is_zero() {
        assert!(FpRepr::from(0).is_zero());
        assert!(!FpRepr::from(1).is_zero());
        assert!(!FpRepr::new([0, 0, 0, 0, 1, 0]).is_zero());
    }

    #[test]
    fn test_fp_repr_div2() {
        let mut a = FpRepr::new([
            0x8b0ad39f8dd7482a,
            0x147221c9a7178b69,
            0x54764cb08d8a6aa0,
            0x8519d708e1d83041,
            0x41f82777bd13fdb,
            0xf43944578f9b771b,
        ]);
        a.div2();
        assert_eq!(
            a,
            FpRepr::new([
                0xc58569cfc6eba415,
                0xa3910e4d38bc5b4,
                0xaa3b265846c53550,
                0xc28ceb8470ec1820,
                0x820fc13bbde89fed,
                0x7a1ca22bc7cdbb8d
            ])
        );
        for _ in 0..10 {
            a.div2();
        }
        assert_eq!(
            a,
            FpRepr::new([
                0x6d31615a73f1bae9,
                0x54028e443934e2f1,
                0x82a8ec99611b14d,
                0xfb70a33ae11c3b06,
                0xe36083f04eef7a27,
                0x1e87288af1f36e
            ])
        );
        for _ in 0..300 {
            a.div2();
        }
        assert_eq!(
            a,
            FpRepr::new([0x7288af1f36ee3608, 0x1e8, 0x0, 0x0, 0x0, 0x0])
        );
        for _ in 0..50 {
            a.div2();
        }
        assert_eq!(a, FpRepr::new([0x7a1ca2, 0x0, 0x0, 0x0, 0x0, 0x0]));
        for _ in 0..22 {
            a.div2();
        }
        assert_eq!(a, FpRepr::new([0x1, 0x0, 0x0, 0x0, 0x0, 0x0]));
        a.div2();
        assert!(a.is_zero());
    }

    #[test]
    fn test_fp_repr_shr() {
        let mut a = FpRepr::new([
            0xaa5cdd6172847ffd,
            0x43242c06aed55287,
            0x9ddd5b312f3dd104,
            0xc5541fd48046b7e7,
            0x16080cf4071e0b05,
            0x1225f2901aea514e,
        ]);
        a.shr(0);
        assert_eq!(
            a,
            FpRepr::new([
                0xaa5cdd6172847ffd,
                0x43242c06aed55287,
                0x9ddd5b312f3dd104,
                0xc5541fd48046b7e7,
                0x16080cf4071e0b05,
                0x1225f2901aea514e
            ])
        );
        a.shr(1);
        assert_eq!(
            a,
            FpRepr::new([
                0xd52e6eb0b9423ffe,
                0x21921603576aa943,
                0xceeead98979ee882,
                0xe2aa0fea40235bf3,
                0xb04067a038f0582,
                0x912f9480d7528a7
            ])
        );
        a.shr(50);
        assert_eq!(
            a,
            FpRepr::new([
                0x8580d5daaa50f54b,
                0xab6625e7ba208864,
                0x83fa9008d6fcf3bb,
                0x19e80e3c160b8aa,
                0xbe52035d4a29c2c1,
                0x244
            ])
        );
        a.shr(130);
        assert_eq!(
            a,
            FpRepr::new([
                0xa0fea40235bf3cee,
                0x4067a038f0582e2a,
                0x2f9480d7528a70b0,
                0x91,
                0x0,
                0x0
            ])
        );
        a.shr(64);
        assert_eq!(
            a,
            FpRepr::new([0x4067a038f0582e2a, 0x2f9480d7528a70b0, 0x91, 0x0, 0x0, 0x0])
        );
    }

    #[test]
    fn test_fp_repr_mul2() {
        let mut a = FpRepr::from(23712937547);
        a.mul2();
        assert_eq!(a, FpRepr::new([0xb0acd6c96, 0x0, 0x0, 0x0, 0x0, 0x0]));
        for _ in 0..60 {
            a.mul2();
        }
        assert_eq!(
            a,
            FpRepr::new([0x6000000000000000, 0xb0acd6c9, 0x0, 0x0, 0x0, 0x0])
        );
        for _ in 0..300 {
            a.mul2();
        }
        assert_eq!(
            a,
            FpRepr::new([0x0, 0x0, 0x0, 0x0, 0x0, 0xcd6c960000000000])
        );
        for _ in 0..17 {
            a.mul2();
        }
        assert_eq!(
            a,
            FpRepr::new([0x0, 0x0, 0x0, 0x0, 0x0, 0x2c00000000000000])
        );
        for _ in 0..6 {
            a.mul2();
        }
        assert!(a.is_zero());
    }

    #[test]
    fn test_fp_repr_num_bits() {
        let mut a = FpRepr::from(0);
        assert_eq!(0, a.num_bits());
        a = FpRepr::from(1);
        for i in 1..385 {
            assert_eq!(i, a.num_bits());
            a.mul2();
        }
        assert_eq!(0, a.num_bits());
    }

    #[test]
    fn test_fp_repr_sub_noborrow() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let mut t = FpRepr::new([
            0x827a4a08041ebd9,
            0x3c239f3dcc8f0d6b,
            0x9ab46a912d555364,
            0x196936b17b43910b,
            0xad0eb3948a5c34fd,
            0xd56f7b5ab8b5ce8,
        ]);
        t.sub_noborrow(&FpRepr::new([
            0xc7867917187ca02b,
            0x5d75679d4911ffef,
            0x8c5b3e48b1a71c15,
            0x6a427ae846fd66aa,
            0x7a37e7265ee1eaf9,
            0x7c0577a26f59d5,
        ]));
        assert!(
            t == FpRepr::new([
                0x40a12b8967c54bae,
                0xdeae37a0837d0d7b,
                0xe592c487bae374e,
                0xaf26bbc934462a61,
                0x32d6cc6e2b7a4a03,
                0xcdaf23e091c0313
            ])
        );

        for _ in 0..1000 {
            let mut a = Fp::random(&mut rng).into_repr();
            a.0.l[5] >>= 30;
            let mut b = a;
            for _ in 0..10 {
                b.mul2();
            }
            let mut c = b;
            for _ in 0..10 {
                c.mul2();
            }

            assert!(a < b);
            assert!(b < c);

            let mut csub_ba = c;
            csub_ba.sub_noborrow(&b);
            csub_ba.sub_noborrow(&a);

            let mut csub_ab = c;
            csub_ab.sub_noborrow(&a);
            csub_ab.sub_noborrow(&b);

            assert_eq!(csub_ab, csub_ba);
        }

        // Subtracting q+1 from q should produce -1 (mod 2**384)
        let mut qplusone = FpRepr::new([
            0xb9feffffffffaaab,
            0x1eabfffeb153ffff,
            0x6730d2a0f6b0f624,
            0x64774b84f38512bf,
            0x4b1ba7b6434bacd7,
            0x1a0111ea397fe69a,
        ]);
        qplusone.sub_noborrow(&FpRepr::new([
            0xb9feffffffffaaac,
            0x1eabfffeb153ffff,
            0x6730d2a0f6b0f624,
            0x64774b84f38512bf,
            0x4b1ba7b6434bacd7,
            0x1a0111ea397fe69a,
        ]));
        assert_eq!(
            qplusone,
            FpRepr::new([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff
            ])
        );
    }

    #[test]
    fn test_fp_repr_add_nocarry() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let mut t = FpRepr::new([
            0x827a4a08041ebd9,
            0x3c239f3dcc8f0d6b,
            0x9ab46a912d555364,
            0x196936b17b43910b,
            0xad0eb3948a5c34fd,
            0xd56f7b5ab8b5ce8,
        ]);
        t.add_nocarry(&FpRepr::new([
            0xc7867917187ca02b,
            0x5d75679d4911ffef,
            0x8c5b3e48b1a71c15,
            0x6a427ae846fd66aa,
            0x7a37e7265ee1eaf9,
            0x7c0577a26f59d5,
        ]));
        assert!(
            t == FpRepr::new([
                0xcfae1db798be8c04,
                0x999906db15a10d5a,
                0x270fa8d9defc6f79,
                0x83abb199c240f7b6,
                0x27469abae93e1ff6,
                0xdd2fd2d4dfab6be
            ])
        );

        // Test for the associativity of addition.
        for _ in 0..1000 {
            let mut a = Fp::random(&mut rng).into_repr();
            let mut b = Fp::random(&mut rng).into_repr();
            let mut c = Fp::random(&mut rng).into_repr();

            // Unset the first few bits, so that overflow won't occur.
            a.0.l[5] >>= 3;
            b.0.l[5] >>= 3;
            c.0.l[5] >>= 3;

            let mut abc = a;
            abc.add_nocarry(&b);
            abc.add_nocarry(&c);

            let mut acb = a;
            acb.add_nocarry(&c);
            acb.add_nocarry(&b);

            let mut bac = b;
            bac.add_nocarry(&a);
            bac.add_nocarry(&c);

            let mut bca = b;
            bca.add_nocarry(&c);
            bca.add_nocarry(&a);

            let mut cab = c;
            cab.add_nocarry(&a);
            cab.add_nocarry(&b);

            let mut cba = c;
            cba.add_nocarry(&b);
            cba.add_nocarry(&a);

            assert_eq!(abc, acb);
            assert_eq!(abc, bac);
            assert_eq!(abc, bca);
            assert_eq!(abc, cab);
            assert_eq!(abc, cba);
        }

        // Adding 1 to (2^384 - 1) should produce zero
        let mut x = FpRepr::new([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
        ]);
        x.add_nocarry(&FpRepr::from(1));
        assert!(x.is_zero());
    }

    #[test]
    fn test_fp_add_assign() {
        {
            // Random number
            let mut tmp = Fp(FpRepr::new([
                0x624434821df92b69,
                0x503260c04fd2e2ea,
                0xd9df726e0d16e8ce,
                0xfbcb39adfd5dfaeb,
                0x86b8a22b0c88b112,
                0x165a2ed809e4201b,
            ])
            .0);
            assert!(!tmp.is_zero());
            // Test that adding zero has no effect.
            tmp.add_assign(&Fp(FpRepr::from(0).0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0x624434821df92b69,
                    0x503260c04fd2e2ea,
                    0xd9df726e0d16e8ce,
                    0xfbcb39adfd5dfaeb,
                    0x86b8a22b0c88b112,
                    0x165a2ed809e4201b
                ])
                .0)
            );
            // Add one and test for the result.
            tmp.add_assign(&Fp(FpRepr::from(1).0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0x624434821df92b6a,
                    0x503260c04fd2e2ea,
                    0xd9df726e0d16e8ce,
                    0xfbcb39adfd5dfaeb,
                    0x86b8a22b0c88b112,
                    0x165a2ed809e4201b
                ])
                .0)
            );
            // Add another random number that exercises the reduction.
            tmp.add_assign(&Fp(FpRepr::new([
                0x374d8f8ea7a648d8,
                0xe318bb0ebb8bfa9b,
                0x613d996f0a95b400,
                0x9fac233cb7e4fef1,
                0x67e47552d253c52,
                0x5c31b227edf25da,
            ])
            .0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0xdf92c410c59fc997,
                    0x149f1bd05a0add85,
                    0xd3ec393c20fba6ab,
                    0x37001165c1bde71d,
                    0x421b41c9f662408e,
                    0x21c38104f435f5b
                ])
                .0)
            );
            // Add one to (q - 1) and test for the result.
            tmp = Fp(FpRepr::new([
                0xb9feffffffffaaaa,
                0x1eabfffeb153ffff,
                0x6730d2a0f6b0f624,
                0x64774b84f38512bf,
                0x4b1ba7b6434bacd7,
                0x1a0111ea397fe69a,
            ])
            .0);
            tmp.add_assign(&Fp(FpRepr::from(1).0));
            assert!(tmp.into_repr().is_zero());
            // Add a random number to another one such that the result is q - 1
            tmp = Fp(FpRepr::new([
                0x531221a410efc95b,
                0x72819306027e9717,
                0x5ecefb937068b746,
                0x97de59cd6feaefd7,
                0xdc35c51158644588,
                0xb2d176c04f2100,
            ])
            .0);
            tmp.add_assign(&Fp(FpRepr::new([
                0x66ecde5bef0fe14f,
                0xac2a6cf8aed568e8,
                0x861d70d86483edd,
                0xcc98f1b7839a22e8,
                0x6ee5e2a4eae7674e,
                0x194e40737930c599,
            ])
            .0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0xb9feffffffffaaaa,
                    0x1eabfffeb153ffff,
                    0x6730d2a0f6b0f624,
                    0x64774b84f38512bf,
                    0x4b1ba7b6434bacd7,
                    0x1a0111ea397fe69a
                ])
                .0)
            );
            // Add one to the result and test for it.
            tmp.add_assign(&Fp(FpRepr::from(1).0));
            assert!(tmp.into_repr().is_zero());
        }

        // Test associativity

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Generate a, b, c and ensure (a + b) + c == a + (b + c).
            let a = Fp::random(&mut rng);
            let b = Fp::random(&mut rng);
            let c = Fp::random(&mut rng);

            let mut tmp1 = a;
            tmp1.add_assign(&b);
            tmp1.add_assign(&c);

            let mut tmp2 = b;
            tmp2.add_assign(&c);
            tmp2.add_assign(&a);

            // assert!(tmp1.is_valid());
            // assert!(tmp2.is_valid());
            assert_eq!(tmp1, tmp2);
        }
    }

    #[test]
    fn test_fp_sub_assign() {
        {
            // Test arbitrary subtraction that tests reduction.
            let mut tmp = Fp(FpRepr::new([
                0x531221a410efc95b,
                0x72819306027e9717,
                0x5ecefb937068b746,
                0x97de59cd6feaefd7,
                0xdc35c51158644588,
                0xb2d176c04f2100,
            ])
            .0);
            tmp.sub_assign(&Fp(FpRepr::new([
                0x98910d20877e4ada,
                0x940c983013f4b8ba,
                0xf677dc9b8345ba33,
                0xbef2ce6b7f577eba,
                0xe1ae288ac3222c44,
                0x5968bb602790806,
            ])
            .0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0x748014838971292c,
                    0xfd20fad49fddde5c,
                    0xcf87f198e3d3f336,
                    0x3d62d6e6e41883db,
                    0x45a3443cd88dc61b,
                    0x151d57aaf755ff94
                ])
                .0)
            );

            // Test the opposite subtraction which doesn't test reduction.
            tmp = Fp(FpRepr::new([
                0x98910d20877e4ada,
                0x940c983013f4b8ba,
                0xf677dc9b8345ba33,
                0xbef2ce6b7f577eba,
                0xe1ae288ac3222c44,
                0x5968bb602790806,
            ])
            .0);
            tmp.sub_assign(&Fp(FpRepr::new([
                0x531221a410efc95b,
                0x72819306027e9717,
                0x5ecefb937068b746,
                0x97de59cd6feaefd7,
                0xdc35c51158644588,
                0xb2d176c04f2100,
            ])
            .0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0x457eeb7c768e817f,
                    0x218b052a117621a3,
                    0x97a8e10812dd02ed,
                    0x2714749e0f6c8ee3,
                    0x57863796abde6bc,
                    0x4e3ba3f4229e706
                ])
                .0)
            );

            // Test for sensible results with zero
            tmp = Fp(FpRepr::from(0).0);
            tmp.sub_assign(&Fp(FpRepr::from(0).0));
            assert!(tmp.is_zero());

            tmp = Fp(FpRepr::new([
                0x98910d20877e4ada,
                0x940c983013f4b8ba,
                0xf677dc9b8345ba33,
                0xbef2ce6b7f577eba,
                0xe1ae288ac3222c44,
                0x5968bb602790806,
            ])
            .0);
            tmp.sub_assign(&Fp(FpRepr::from(0).0));
            assert_eq!(
                tmp,
                Fp(FpRepr::new([
                    0x98910d20877e4ada,
                    0x940c983013f4b8ba,
                    0xf677dc9b8345ba33,
                    0xbef2ce6b7f577eba,
                    0xe1ae288ac3222c44,
                    0x5968bb602790806
                ])
                .0)
            );
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure that (a - b) + (b - a) = 0.
            let a = Fp::random(&mut rng);
            let b = Fp::random(&mut rng);

            let mut tmp1 = a;
            tmp1.sub_assign(&b);

            let mut tmp2 = b;
            tmp2.sub_assign(&a);

            tmp1.add_assign(&tmp2);
            assert!(tmp1.is_zero());
        }
    }

    #[test]
    fn test_fp_mul_assign() {
        assert_eq!(
            Fp(FpRepr::new([
                0xcc6200000020aa8a,
                0x422800801dd8001a,
                0x7f4f5e619041c62c,
                0x8a55171ac70ed2ba,
                0x3f69cc3a3d07d58b,
                0xb972455fd09b8ef,
            ])
            .0) * &Fp(FpRepr::new([
                0x329300000030ffcf,
                0x633c00c02cc40028,
                0xbef70d925862a942,
                0x4f7fa2a82a963c17,
                0xdf1eb2575b8bc051,
                0x1162b680fb8e9566,
            ])
            .0),
            Fp(FpRepr::new([
                0x9dc4000001ebfe14,
                0x2850078997b00193,
                0xa8197f1abb4d7bf,
                0xc0309573f4bfe871,
                0xf48d0923ffaf7620,
                0x11d4b58c7a926e66
            ])
            .0)
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * b) * c = a * (b * c)
            let a = Fp::random(&mut rng);
            let b = Fp::random(&mut rng);
            let c = Fp::random(&mut rng);

            let mut tmp1 = a;
            tmp1.mul_assign(&b);
            tmp1.mul_assign(&c);

            let mut tmp2 = b;
            tmp2.mul_assign(&c);
            tmp2.mul_assign(&a);

            assert_eq!(tmp1, tmp2);
        }

        for _ in 0..1000000 {
            // Ensure that r * (a + b + c) = r*a + r*b + r*c

            let r = Fp::random(&mut rng);
            let mut a = Fp::random(&mut rng);
            let mut b = Fp::random(&mut rng);
            let mut c = Fp::random(&mut rng);

            let mut tmp1 = a;
            tmp1.add_assign(&b);
            tmp1.add_assign(&c);
            tmp1.mul_assign(&r);

            a.mul_assign(&r);
            b.mul_assign(&r);
            c.mul_assign(&r);

            a.add_assign(&b);
            a.add_assign(&c);

            assert_eq!(tmp1, a);
        }
    }

    #[test]
    fn test_fp_squaring() {
        let mut a = Fp(FpRepr::new([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0x19ffffffffffffff,
        ])
        .0);
        assert!(!a.is_zero());
        a.square();
        assert_eq!(
            a,
            Fp::from_repr(FpRepr::new([
                0x1cfb28fe7dfbbb86,
                0x24cbe1731577a59,
                0xcce1d4edc120e66e,
                0xdc05c659b4e15b27,
                0x79361e5a802c6a23,
                0x24bcbe5d51b9a6f
            ]))
            .unwrap()
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * a) = a^2
            let a = Fp::random(&mut rng);

            let mut tmp = a;
            tmp.square();

            let mut tmp2 = a;
            tmp2.mul_assign(&a);

            assert_eq!(tmp, tmp2);
        }
    }

    #[test]
    fn test_fp_inverse() {
        assert!(Fp::zero().inverse().is_none());

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let one = Fp::one();

        for _ in 0..1000 {
            // Ensure that a * a^-1 = 1
            let mut a = Fp::random(&mut rng);
            let ainv = a.inverse().unwrap();
            a.mul_assign(&ainv);
            assert_eq!(a, one);
        }
    }

    #[test]
    fn test_fp_double() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure doubling a is equivalent to adding a to itself.
            let mut a = Fp::random(&mut rng);
            let mut b = a;
            b.add_assign(&a);
            a.double();
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_fp_negate() {
        {
            let mut a = Fp::zero();
            a.negate();

            assert!(a.is_zero());
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure (a - (-a)) = 0.
            let mut a = Fp::random(&mut rng);
            let mut b = a;
            b.negate();
            a.add_assign(&b);

            assert!(a.is_zero());
        }
    }

    #[test]
    fn test_fp_pow() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..1000 {
            // Exponentiate by various small numbers and ensure it consists with repeated
            // multiplication.
            let a = Fp::random(&mut rng);
            let target = a.pow(&[i]);
            let mut c = Fp::one();
            for _ in 0..i {
                c.mul_assign(&a);
            }
            assert_eq!(c, target);
        }

        for _ in 0..1000 {
            // Exponentiating by the modulus should have no effect in a prime field.
            let a = Fp::random(&mut rng);

            assert_eq!(a, a.pow(Fp::char()));
        }
    }

    #[test]
    fn test_fp_sqrt() {
        use fff::SqrtField;

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        assert_eq!(Fp::zero().sqrt().unwrap(), Fp::zero());

        for _ in 0..1000 {
            // Ensure sqrt(a^2) = a or -a
            let a = Fp::random(&mut rng);
            let mut nega = a;
            nega.negate();
            let mut b = a;
            b.square();

            let b = b.sqrt().unwrap();

            assert!(a == b || nega == b);
        }

        for _ in 0..1000 {
            // Ensure sqrt(a)^2 = a for random a
            let a = Fp::random(&mut rng);

            if let Some(mut tmp) = a.sqrt() {
                tmp.square();

                assert_eq!(a, tmp);
            }
        }
    }

    #[test]
    fn test_fp_from_into_repr() {
        // q + 1 should not be in the field
        assert!(Fp::from_repr(FpRepr::new([
            0xb9feffffffffaaac,
            0x1eabfffeb153ffff,
            0x6730d2a0f6b0f624,
            0x64774b84f38512bf,
            0x4b1ba7b6434bacd7,
            0x1a0111ea397fe69a
        ]))
        .is_err());

        // q should not be in the field
        assert!(Fp::from_repr(Fp::char()).is_err());

        // Multiply some arbitrary representations to see if the result is as expected.
        let a = FpRepr::new([
            0x4a49dad4ff6cde2d,
            0xac62a82a8f51cd50,
            0x2b1f41ab9f36d640,
            0x908a387f480735f1,
            0xae30740c08a875d7,
            0x6c80918a365ef78,
        ]);
        let mut a_fp = Fp::from_repr(a).unwrap();
        let b = FpRepr::new([
            0xbba57917c32f0cf0,
            0xe7f878cf87f05e5d,
            0x9498b4292fd27459,
            0xd59fd94ee4572cfa,
            0x1f607186d5bb0059,
            0xb13955f5ac7f6a3,
        ]);
        let b_fp = Fp::from_repr(b).unwrap();
        let c = FpRepr::new([
            0xf5f70713b717914c,
            0x355ea5ac64cbbab1,
            0xce60dd43417ec960,
            0xf16b9d77b0ad7d10,
            0xa44c204c1de7cdb7,
            0x1684487772bc9a5a,
        ]);
        a_fp.mul_assign(&b_fp);
        assert_eq!(a_fp.into_repr(), c);

        // Zero should be in the field.
        assert!(Fp::from_repr(FpRepr::from(0)).unwrap().is_zero());

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Try to turn Fp elements into representations and back again, and compare.
            let a = Fp::random(&mut rng);
            let a_repr = a.into_repr();
            let b_repr = FpRepr::from(a);
            assert_eq!(a_repr, b_repr);
            let a_again = Fp::from_repr(a_repr).unwrap();

            assert_eq!(a, a_again);
        }
    }

    #[test]
    fn test_fp_repr_display() {
        assert_eq!(
        format!("{}", FpRepr::new([0xa956babf9301ea24, 0x39a8f184f3535c7b, 0xb38d35b3f6779585, 0x676cc4eef4c46f2c, 0xb1d4aad87651e694, 0x1947f0d5f4fe325a])),
        "0x1947f0d5f4fe325ab1d4aad87651e694676cc4eef4c46f2cb38d35b3f677958539a8f184f3535c7ba956babf9301ea24".to_string()
    );
        assert_eq!(
        format!("{}", FpRepr::new([0xb4171485fd8622dd, 0x864229a6edec7ec5, 0xc57f7bdcf8dfb707, 0x6db7ff0ecea4584a, 0xf8d8578c4a57132d, 0x6eb66d42d9fcaaa])),
        "0x06eb66d42d9fcaaaf8d8578c4a57132d6db7ff0ecea4584ac57f7bdcf8dfb707864229a6edec7ec5b4171485fd8622dd".to_string()
    );
        assert_eq!(
        format!("{}", FpRepr::new([0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff])),
        "0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff".to_string()
    );
        assert_eq!(
        format!("{}", FpRepr::new([0, 0, 0, 0, 0, 0])),
        "0x000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000".to_string()
    );
    }

    #[test]
    fn test_fp_display() {
        assert_eq!(
            format!("{}", Fp::from_repr(FpRepr::new([
                0xa956babf9301ea24,
                0x39a8f184f3535c7b,
                0xb38d35b3f6779585,
                0x676cc4eef4c46f2c,
                0xb1d4aad87651e694,
                0x1947f0d5f4fe325a
            ])
            ).unwrap()),
            "Fp(0x1947f0d5f4fe325ab1d4aad87651e694676cc4eef4c46f2cb38d35b3f677958539a8f184f3535c7ba956babf9301ea24)".to_string()
    );
        assert_eq!(
        format!("{}", Fp::from_repr(FpRepr::new([0xe28e79396ac2bbf8, 0x413f6f7f06ea87eb, 0xa4b62af4a792a689, 0xb7f89f88f59c1dc5, 0x9a551859b1e43a9a, 0x6c9f5a1060de974])).unwrap()),
        "Fp(0x06c9f5a1060de9749a551859b1e43a9ab7f89f88f59c1dc5a4b62af4a792a689413f6f7f06ea87ebe28e79396ac2bbf8)".to_string()
    );
    }

    #[test]
    fn test_fp_num_bits() {
        assert_eq!(Fp::NUM_BITS, 381);
        assert_eq!(Fp::CAPACITY, 380);
    }

    #[test]
    fn test_fp_root_of_unity() {
        use fff::SqrtField;

        assert_eq!(Fp::S, 1);
        assert_eq!(
            Fp::multiplicative_generator(),
            Fp::from_repr(FpRepr::from(2)).unwrap()
        );
        assert_eq!(
            Fp::multiplicative_generator().pow(&[
                0xdcff7fffffffd555,
                0xf55ffff58a9ffff,
                0xb39869507b587b12,
                0xb23ba5c279c2895f,
                0x258dd3db21a5d66b,
                0xd0088f51cbff34d
            ]),
            Fp::root_of_unity()
        );
        assert_eq!(
            Fp::root_of_unity().pow(&[1 << Fp::S, 0, 0, 0, 0, 0]),
            Fp::one()
        );
        assert!(Fp::multiplicative_generator().sqrt().is_none());
    }

    #[test]
    fn fp_field_tests() {
        crate::tests::field::random_field_tests::<Fp>();
        crate::tests::field::random_sqrt_tests::<Fp>();
        crate::tests::field::random_frobenius_tests::<Fp, _>(Fp::char(), 13);
        crate::tests::field::from_str_tests::<Fp>();
    }

    #[test]
    fn test_fp_ordering() {
        // FpRepr's ordering is well-tested, but we still need to make sure the Fp
        // elements aren't being compared in Montgomery form.
        for i in 0..100 {
            let a = FpRepr::from(i + 1);
            let b = FpRepr::from(i);
            assert!(
                Fp::from_repr(a).unwrap() > Fp::from_repr(b).unwrap(),
                "{}: {:?} > {:?}",
                i,
                a,
                b
            );
        }
    }

    #[test]
    fn fp_repr_tests() {
        crate::tests::repr::random_repr_tests::<Fp>();
    }

    #[test]
    fn test_fp_repr_conversion() {
        let a = Fp::from(1);
        let b = Fp::from_repr(FpRepr::from(1)).unwrap();
        assert_eq!(a, b);
        assert_eq!(Fp::from(1).into_repr(), FpRepr::from(1));

        let a = Fp::from(12);
        assert_eq!(a, Fp::from_repr(a.into_repr()).unwrap());

        let a = FpRepr::from(12);
        assert_eq!(Fp::from_repr(a).unwrap().into_repr(), a);
    }
}
