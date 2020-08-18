//! An implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::{
    cmp,
    convert::TryInto,
    fmt, mem,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use fff::{Field, PrimeField};

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
///
/// The inner representation is stored in Montgomery form.
#[derive(Default, Clone, Copy)]
pub struct Scalar(blst_fr);

/// Representation of a `Scalar`, in regular coordinates.
#[derive(Default, Clone, Copy)]
pub struct ScalarRepr(pub(crate) blst_fr);

impl AsRef<[u64]> for ScalarRepr {
    fn as_ref(&self) -> &[u64] {
        &self.0.l
    }
}

impl AsMut<[u64]> for ScalarRepr {
    fn as_mut(&mut self) -> &mut [u64] {
        &mut self.0.l
    }
}

const LIMBS: usize = 4;
const LIMB_BITS: usize = 64;

const MODULUS: ScalarRepr = ScalarRepr(blst_fr {
    l: [
        0xffffffff00000001,
        0x53bda402fffe5bfe,
        0x3339d80809a1d805,
        0x73eda753299d7d48,
    ],
});

/// R = 2^256 mod q
const R: Scalar = Scalar(blst_fr {
    l: [
        0x00000001fffffffe,
        0x5884b7fa00034802,
        0x998c4fefecbc4ff5,
        0x1824b159acc5056f,
    ],
});

impl fmt::Debug for ScalarRepr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x")?;
        for &b in self.0.l.iter().rev() {
            write!(f, "{:016x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for ScalarRepr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x")?;
        for &b in self.0.l.iter().rev() {
            write!(f, "{:016x}", b)?;
        }
        Ok(())
    }
}

impl From<u64> for ScalarRepr {
    fn from(val: u64) -> ScalarRepr {
        ScalarRepr(blst_fr { l: [val, 0, 0, 0] })
    }
}

impl Ord for ScalarRepr {
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        for (a, b) in self.0.l.iter().rev().zip(other.0.l.iter().rev()) {
            match a.cmp(b) {
                cmp::Ordering::Greater => {
                    return cmp::Ordering::Greater;
                }
                cmp::Ordering::Less => {
                    return cmp::Ordering::Less;
                }
                cmp::Ordering::Equal => {}
            }
        }

        cmp::Ordering::Equal
    }
}

impl PartialOrd for ScalarRepr {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for ScalarRepr {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l == other.0.l
    }
}
impl Eq for ScalarRepr {}

impl fff::PrimeFieldRepr for ScalarRepr {
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
                mem::swap(&mut t, i);
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
                mem::swap(&mut t, i);
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

pub const S: u32 = 32;

impl fmt::Debug for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Scalar({:?})", self.into_repr())
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Scalar({})", self.into_repr())
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l[0] == other.0.l[0]
            && self.0.l[1] == other.0.l[1]
            && self.0.l[2] == other.0.l[2]
            && self.0.l[3] == other.0.l[3]
    }
}
impl Eq for Scalar {}

impl From<Scalar> for ScalarRepr {
    fn from(val: Scalar) -> Self {
        let mut out = blst_fr::default();
        unsafe { blst_fr_from(&mut out, &val.0) };
        ScalarRepr(out)
    }
}

impl From<Scalar> for blst_fr {
    fn from(val: Scalar) -> blst_fr {
        val.0
    }
}

impl From<blst_fr> for Scalar {
    fn from(val: blst_fr) -> Scalar {
        Scalar(val)
    }
}

impl From<u64> for Scalar {
    fn from(val: u64) -> Scalar {
        // safe because u64 is always in the field
        let raw = blst_fr { l: [val, 0, 0, 0] };
        let mut out = blst_fr::default();
        unsafe { blst_fr_to(&mut out, &raw) };

        Scalar(out)
    }
}

impl<'a> Neg for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        let mut out = blst_fr::default();

        const FLAG: usize = 0x1;

        unsafe { blst_fr_cneg(&mut out, &self.0, FLAG) };

        Scalar(out)
    }
}

impl Neg for Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &'b Scalar) -> Scalar {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &'b Scalar) -> Scalar {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &'b Scalar) -> Scalar {
        self.mul(rhs)
    }
}

impl_binops_additive!(Scalar, Scalar);
impl_binops_multiplicative!(Scalar, Scalar);

impl fff::Field for Scalar {
    fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        // The number of bits we should "shave" from a randomly sampled reputation.
        const REPR_SHAVE_BITS: usize = 256 - Scalar::NUM_BITS as usize;

        loop {
            let mut raw = blst_fr::default();
            for i in 0..4 {
                raw.l[i] = rng.next_u64();
            }

            // Mask away the unused most-significant bits.
            raw.l[3] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if ScalarRepr(raw) < MODULUS {
                return Scalar(raw);
            }
        }
    }

    fn zero() -> Self {
        Scalar(blst_fr::default())
    }

    fn one() -> Self {
        R
    }

    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    fn square(&mut self) {
        let mut raw = blst_fr::default();
        unsafe { blst_fr_sqr(&mut raw, &self.0) }

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
        *self *= other;
    }

    fn inverse(&self) -> Option<Self> {
        #[inline(always)]
        fn square_assign_multi(n: &mut Scalar, num_times: usize) {
            for _ in 0..num_times {
                n.square();
            }
        }
        // found using https://github.com/kwantam/addchain
        let mut t0 = *self;
        t0.square();
        let mut t1 = t0 * self;
        let mut t16 = t0;
        t16.square();
        let mut t6 = t16;
        t6.square();
        let mut t5 = t6 * t0;
        t0 = t6 * t16;
        let mut t12 = t5 * t16;
        let mut t2 = t6;
        t2.square();
        let mut t7 = t5 * t6;
        let mut t15 = t0 * t5;
        let mut t17 = t12;
        t17.square();
        t1 *= t17;
        let mut t3 = t7 * t2;
        let t8 = t1 * t17;
        let t4 = t8 * t2;
        let t9 = t8 * t7;
        t7 = t4 * t5;
        let t11 = t4 * t17;
        t5 = t9 * t17;
        let t14 = t7 * t15;
        let t13 = t11 * t12;
        t12 = t11 * t17;
        t15 *= &t12;
        t16 *= &t15;
        t3 *= &t16;
        t17 *= &t3;
        t0 *= &t17;
        t6 *= &t0;
        t2 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t17;
        square_assign_multi(&mut t0, 9);
        t0 *= &t16;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 9);
        t0 *= &t15;
        square_assign_multi(&mut t0, 7);
        t0 *= &t14;
        square_assign_multi(&mut t0, 7);
        t0 *= &t13;
        square_assign_multi(&mut t0, 10);
        t0 *= &t12;
        square_assign_multi(&mut t0, 9);
        t0 *= &t11;
        square_assign_multi(&mut t0, 8);
        t0 *= &t8;
        square_assign_multi(&mut t0, 8);
        t0 *= self;
        square_assign_multi(&mut t0, 14);
        t0 *= &t9;
        square_assign_multi(&mut t0, 10);
        t0 *= &t8;
        square_assign_multi(&mut t0, 15);
        t0 *= &t7;
        square_assign_multi(&mut t0, 10);
        t0 *= &t6;
        square_assign_multi(&mut t0, 8);
        t0 *= &t5;
        square_assign_multi(&mut t0, 16);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 7);
        t0 *= &t4;
        square_assign_multi(&mut t0, 9);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t3;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 8);
        t0 *= &t2;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;
        square_assign_multi(&mut t0, 5);
        t0 *= &t1;

        if self.is_zero() {
            None
        } else {
            Some(t0)
        }
    }

    fn frobenius_map(&mut self, _: usize) {
        // This has no effect in a prime field.
    }
}

impl ScalarRepr {
    pub fn new(raw: [u64; 4]) -> Self {
        ScalarRepr(blst_fr { l: raw })
    }
}

impl fff::PrimeField for Scalar {
    type Repr = ScalarRepr;

    const NUM_BITS: u32 = 255;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const S: u32 = S;

    fn from_repr(repr: Self::Repr) -> Result<Self, fff::PrimeFieldDecodingError> {
        if ScalarRepr(repr.0) < MODULUS {
            let mut out = blst_fr::default();
            unsafe { blst_fr_to(&mut out, &repr.0) }
            Ok(Scalar(out))
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
        Scalar::from_repr(ScalarRepr::new([7, 0, 0, 0])).unwrap()
    }

    fn root_of_unity() -> Self {
        Scalar(blst_fr {
            l: [
                0xb9b58d8c5f0e466a,
                0x5b1b4c801819d7ec,
                0x0af53ae352a31e64,
                0x5bf3adda19e9b27b,
            ],
        })
    }
}

impl fff::SqrtField for Scalar {
    fn legendre(&self) -> fff::LegendreSymbol {
        const MOD_MINUS_1_OVER_2: [u64; 4] = [
            0x7fffffff80000000,
            0xa9ded2017fff2dff,
            0x199cec0404d0ec02,
            0x39f6d3a994cebea4,
        ];
        // s = self^((modulus - 1) // 2)
        let s = self.pow(MOD_MINUS_1_OVER_2);
        if s == Self::zero() {
            fff::LegendreSymbol::Zero
        } else if s == Self::one() {
            fff::LegendreSymbol::QuadraticResidue
        } else {
            fff::LegendreSymbol::QuadraticNonResidue
        }
    }

    fn sqrt(&self) -> Option<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // https://eprint.iacr.org/2012/685.pdf (page 12, algorithm 5)
        match self.legendre() {
            fff::LegendreSymbol::Zero => Some(*self),
            fff::LegendreSymbol::QuadraticNonResidue => None,
            fff::LegendreSymbol::QuadraticResidue => {
                let mut c = Self::root_of_unity();
                let mut r = self.pow([
                    9223141137265459200,
                    347036667491570177,
                    10722717374829358084,
                    972477353,
                ]);
                let mut t = self.pow([
                    18446282274530918399,
                    694073334983140354,
                    2998690675949164552,
                    1944954707,
                ]);
                let mut m = S;

                while t != Self::one() {
                    let mut i = 1;
                    {
                        let mut t2i = t;
                        t2i.square();
                        loop {
                            if t2i == Self::one() {
                                break;
                            }
                            t2i.square();
                            i += 1;
                        }
                    }

                    for _ in 0..(m - i - 1) {
                        c.square();
                    }
                    r *= &c;
                    c.square();
                    t *= &c;
                    m = i;
                }

                Some(r)
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct NotInFieldError;

impl fmt::Display for NotInFieldError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Not in field")
    }
}

impl std::error::Error for NotInFieldError {}

impl TryInto<Scalar> for blst_scalar {
    type Error = NotInFieldError;

    fn try_into(self) -> Result<Scalar, Self::Error> {
        if !unsafe { blst_scalar_fr_check(&self) } {
            return Err(NotInFieldError);
        }

        // Safe because valid fr check was just made above.
        let fr: blst_fr = unsafe { std::mem::transmute(self) };

        let mut out = blst_fr::default();
        // convert to montgomery
        unsafe { blst_fr_to(&mut out, &fr) }

        Ok(Scalar(out))
    }
}

impl Scalar {
    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 32]) -> Option<Scalar> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_scalar::default();

        unsafe {
            blst_scalar_from_lendian(&mut raw, in_v.as_mut_ptr());
        }

        raw.try_into().ok()
    }

    /// Attempts to convert a big-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_be(bytes: &[u8; 32]) -> Option<Scalar> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_scalar::default();

        unsafe {
            blst_scalar_from_bendian(&mut raw, in_v.as_mut_ptr());
        }

        raw.try_into().ok()
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Scalar` representation.
    pub fn from_raw(val: [u64; 4]) -> Self {
        let original = blst_fr { l: val };

        let mut raw = blst_fr::default();

        // Convert to montgomery form
        unsafe { blst_fr_to(&mut raw, &original) }

        Scalar(raw)
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 32] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 32];

        let mut out = blst_fr::default();
        unsafe {
            // transform out of montgomery
            blst_fr_from(&mut out, &self.0);
        }

        // Safe because any valid blst_fr is also a valid blst_scalar.
        let scalar: blst_scalar = unsafe { std::mem::transmute(out) };

        unsafe {
            blst_lendian_from_scalar(out_v.as_mut_ptr(), &scalar);
        }

        let mut out = [0u8; 32];
        out.copy_from_slice(&out_v);

        out
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes_be(&self) -> [u8; 32] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 32];

        let mut out = blst_fr::default();
        unsafe {
            // transform out of montgomery
            blst_fr_from(&mut out, &self.0);
        }

        // Safe because any valid blst_fr is also a valid blst_scalar.
        let scalar: blst_scalar = unsafe { std::mem::transmute(out) };
        unsafe {
            blst_bendian_from_scalar(out_v.as_mut_ptr(), &scalar);
        }
        let mut out = [0u8; 32];
        out.copy_from_slice(&out_v);

        out
    }

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_mul(&mut out, &self.0, &rhs.0) };

        Scalar(out)
    }

    /// Subtracts `rhs` from `self`, returning the result.
    #[inline]
    pub fn sub(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_sub(&mut out, &self.0, &rhs.0) };

        Scalar(out)
    }

    /// Adds `rhs` to `self`, returning the result.
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_add(&mut out, &self.0, &rhs.0) };

        Scalar(out)
    }

    /// Returns true if this element is zero.
    pub fn is_zero(&self) -> bool {
        self.0.l.iter().all(|&e| e == 0)
    }

    /// Multiplies `self` with `3`, returning the result.
    pub fn mul3(&self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_mul_by_3(&mut out, &self.0) };

        Scalar(out)
    }

    /// Left shift `self` by `count`, returning the result.
    pub fn shl(&self, count: usize) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_lshift(&mut out, &self.0, count) };

        Scalar(out)
    }

    /// Right shift `self` by `count`, returning the result.
    pub fn shr(&self, count: usize) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_rshift(&mut out, &self.0, count) };

        Scalar(out)
    }
}

#[cfg(test)]
mod tests {
    use super::{Scalar, ScalarRepr, MODULUS, R};

    use fff::{Field, PrimeField, PrimeFieldRepr, SqrtField};
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    /// INV = -(q^{-1} mod 2^64) mod 2^64
    const INV: u64 = 0xfffffffeffffffff;

    /// R^2 = 2^512 mod q
    const R2: Scalar = Scalar(blst::blst_fr {
        l: [
            0xc999e990f3f29c6d,
            0x2b6cedcb87925c23,
            0x05d314967254398f,
            0x0748d9d99f59ff11,
        ],
    });

    const LARGEST: Scalar = Scalar(blst::blst_fr {
        l: [
            0xffffffff00000000,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ],
    });

    #[test]
    fn test_inv() {
        // Compute -(q^{-1} mod 2^64) mod 2^64 by exponentiating
        // by totient(2**64) - 1

        let mut inv = 1u64;
        for _ in 0..63 {
            inv = inv.wrapping_mul(inv);
            inv = inv.wrapping_mul(MODULUS.0.l[0]);
        }
        inv = inv.wrapping_neg();

        assert_eq!(inv, INV);
    }

    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Scalar::zero()),
            "Scalar(0x0000000000000000000000000000000000000000000000000000000000000000)"
        );
        assert_eq!(
            format!("{:?}", Scalar::one()),
            "Scalar(0x0000000000000000000000000000000000000000000000000000000000000001)"
        );
        assert_eq!(
            format!("{:?}", R2),
            "Scalar(0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe)"
        );
    }

    #[test]
    fn test_equality() {
        assert_eq!(Scalar::zero(), Scalar::zero());
        assert_eq!(Scalar::one(), Scalar::one());
        assert_eq!(Scalar(R2.0), Scalar(R2.0));

        assert!(Scalar::zero() != Scalar::one());
        assert!(Scalar::one() != Scalar(R2.0));
    }

    #[test]
    fn test_to_bytes() {
        assert_eq!(
            Scalar::zero().to_bytes_le(),
            [
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            Scalar::one().to_bytes_le(),
            [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ]
        );

        assert_eq!(
            R2.to_bytes_le(),
            [
                254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88, 245, 79, 188, 236,
                239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
            ]
        );

        assert_eq!(
            (-&Scalar::one()).to_bytes_le(),
            [
                0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ]
        );
    }

    #[test]
    fn test_from_bytes() {
        assert_eq!(
            Scalar::from_bytes_le(&[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::zero()
        );

        assert_eq!(
            Scalar::from_bytes_le(&[
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0
            ])
            .unwrap(),
            Scalar::one()
        );

        assert_eq!(
            Scalar::from_bytes_le(&[
                254, 255, 255, 255, 1, 0, 0, 0, 2, 72, 3, 0, 250, 183, 132, 88, 245, 79, 188, 236,
                239, 79, 140, 153, 111, 5, 197, 172, 89, 177, 36, 24
            ])
            .unwrap(),
            R2
        );

        // -1 should work
        assert!(Scalar::from_bytes_le(&[
            0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_some());

        // modulus is invalid
        assert!(Scalar::from_bytes_le(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none());

        // Anything larger than the modulus is invalid
        assert!(Scalar::from_bytes_le(&[
            2, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none());
        assert!(Scalar::from_bytes_le(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 58, 51, 72, 125, 157, 41, 83, 167, 237, 115
        ])
        .is_none());
        assert!(Scalar::from_bytes_le(&[
            1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9, 8,
            216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 116
        ])
        .is_none());
    }

    #[test]
    fn test_zero() {
        assert_eq!(Scalar::zero(), -&Scalar::zero());
        assert_eq!(Scalar::zero(), Scalar::zero() + Scalar::zero());
        assert_eq!(Scalar::zero(), Scalar::zero() - Scalar::zero());
        assert_eq!(Scalar::zero(), Scalar::zero() * Scalar::zero());
    }

    #[test]
    fn test_addition() {
        let mut tmp = LARGEST;
        tmp += &LARGEST;

        assert_eq!(
            tmp,
            Scalar(
                ScalarRepr::new([
                    0xfffffffeffffffff,
                    0x53bda402fffe5bfe,
                    0x3339d80809a1d805,
                    0x73eda753299d7d48
                ])
                .0
            )
        );

        let mut tmp = LARGEST;
        tmp += &Scalar(ScalarRepr::new([1, 0, 0, 0]).0);

        assert_eq!(tmp, Scalar::zero());
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar(ScalarRepr::new([1, 0, 0, 0]).0));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar(ScalarRepr::new([1, 0, 0, 0]).0);
        assert_eq!(tmp, LARGEST);

        {
            let mut a = Scalar::zero();
            a = -a;

            assert!(a.is_zero());
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure (a - (-a)) = 0.
            let mut a = Scalar::random(&mut rng);
            let mut b = a;
            b = -b;
            a += &b;

            assert!(a.is_zero());
        }
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Scalar::zero());

        let mut tmp = Scalar::zero();
        tmp -= &LARGEST;

        let mut tmp2 = Scalar(MODULUS.0);
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_multiplication() {
        let mut tmp = Scalar(
            ScalarRepr::new([
                0x6b7e9b8faeefc81a,
                0xe30a8463f348ba42,
                0xeff3cb67a8279c9c,
                0x3d303651bd7c774d,
            ])
            .0,
        );
        tmp *= &Scalar(
            ScalarRepr::new([
                0x13ae28e3bc35ebeb,
                0xa10f4488075cae2c,
                0x8160e95a853c3b5d,
                0x5ae3f03b561a841d,
            ])
            .0,
        );
        assert!(
            tmp == Scalar(
                ScalarRepr::new([
                    0x23717213ce710f71,
                    0xdbee1fe53a16e1af,
                    0xf565d3e1c2a48000,
                    0x4426507ee75df9d7
                ])
                .0
            )
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * b) * c = a * (b * c)
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);
            let c = Scalar::random(&mut rng);

            let mut tmp1 = a;
            tmp1 *= &b;
            tmp1 *= &c;

            let mut tmp2 = b;
            tmp2 *= &c;
            tmp2 *= &a;

            assert_eq!(tmp1, tmp2);
        }

        for _ in 0..1000000 {
            // Ensure that r * (a + b + c) = r*a + r*b + r*c

            let r = Scalar::random(&mut rng);
            let mut a = Scalar::random(&mut rng);
            let mut b = Scalar::random(&mut rng);
            let mut c = Scalar::random(&mut rng);

            let mut tmp1 = a;
            tmp1 += &b;
            tmp1 += &c;
            tmp1 *= &r;

            a *= &r;
            b *= &r;
            c *= &r;

            a += &b;
            a += &c;

            assert_eq!(tmp1, a);
        }
    }

    #[test]
    fn test_inverse_is_pow() {
        let q_minus_2 = [
            0xfffffffeffffffff,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ];

        let mut r1 = R;
        let mut r2 = R;

        for _ in 0..100 {
            r1 = r1.inverse().unwrap();
            r2 = r2.pow(&q_minus_2);

            assert_eq!(r1, r2);
            // Add R so we check something different next time around
            r1 += &R;
            r2 = r1;
        }
    }

    #[test]
    fn test_sqrt() {
        {
            assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        }

        let mut square = Scalar(
            ScalarRepr::new([
                0x46cd85a5f273077e,
                0x1d30c47dd68fc735,
                0x77f656f60beca0eb,
                0x494aa01bdf32468d,
            ])
            .0,
        );

        let mut none_count = 0;

        for _ in 0..100 {
            let square_root = square.sqrt();
            if square_root.is_none() {
                none_count += 1;
            } else {
                assert_eq!(square_root.unwrap() * square_root.unwrap(), square);
            }
            square -= Scalar::one();
        }

        assert_eq!(49, none_count);
    }

    #[test]
    fn test_from_raw() {
        assert_eq!(
            Scalar::from_raw([
                0x1fffffffd,
                0x5884b7fa00034802,
                0x998c4fefecbc4ff5,
                0x1824b159acc5056f
            ]),
            Scalar::from_raw([0xffffffffffffffff; 4])
        );

        assert_eq!(Scalar::from_raw(MODULUS.0.l), Scalar::zero());

        assert_eq!(Scalar::from_raw([1, 0, 0, 0]), R);
    }

    #[test]
    fn test_double() {
        let a = Scalar::from_raw([
            0x1fff3231233ffffd,
            0x4884b7fa00034802,
            0x998c4fefecbc4ff3,
            0x1824b159acc50562,
        ]);

        let mut b = a;
        b.double();
        assert_eq!(b, a + a);
    }

    #[test]
    fn test_scalar_repr_ordering() {
        fn assert_equality(a: ScalarRepr, b: ScalarRepr) {
            assert_eq!(a, b);
            assert!(a.cmp(&b) == ::std::cmp::Ordering::Equal);
        }

        fn assert_lt(a: ScalarRepr, b: ScalarRepr) {
            assert!(a < b);
            assert!(b > a);
        }

        assert_equality(
            ScalarRepr::new([9999, 9999, 9999, 9999]),
            ScalarRepr::new([9999, 9999, 9999, 9999]),
        );
        assert_equality(
            ScalarRepr::new([9999, 9998, 9999, 9999]),
            ScalarRepr::new([9999, 9998, 9999, 9999]),
        );
        assert_equality(
            ScalarRepr::new([9999, 9999, 9999, 9997]),
            ScalarRepr::new([9999, 9999, 9999, 9997]),
        );
        assert_lt(
            ScalarRepr::new([9999, 9997, 9999, 9998]),
            ScalarRepr::new([9999, 9997, 9999, 9999]),
        );
        assert_lt(
            ScalarRepr::new([9999, 9997, 9998, 9999]),
            ScalarRepr::new([9999, 9997, 9999, 9999]),
        );
        assert_lt(
            ScalarRepr::new([9, 9999, 9999, 9997]),
            ScalarRepr::new([9999, 9999, 9999, 9997]),
        );
    }

    #[test]
    fn test_scalar_repr_from() {
        assert_eq!(ScalarRepr::from(100), ScalarRepr::new([100, 0, 0, 0]));
    }

    #[test]
    fn test_scalar_repr_is_odd() {
        assert!(!ScalarRepr::from(0).is_odd());
        assert!(ScalarRepr::from(0).is_even());
        assert!(ScalarRepr::from(1).is_odd());
        assert!(!ScalarRepr::from(1).is_even());
        assert!(!ScalarRepr::from(324834872).is_odd());
        assert!(ScalarRepr::from(324834872).is_even());
        assert!(ScalarRepr::from(324834873).is_odd());
        assert!(!ScalarRepr::from(324834873).is_even());
    }

    #[test]
    fn test_scalar_repr_is_zero() {
        assert!(ScalarRepr::from(0).is_zero());
        assert!(!ScalarRepr::from(1).is_zero());
        assert!(!ScalarRepr::new([0, 0, 1, 0]).is_zero());
    }

    #[test]
    fn test_scalar_repr_div2() {
        let mut a = ScalarRepr::new([
            0xbd2920b19c972321,
            0x174ed0466a3be37e,
            0xd468d5e3b551f0b5,
            0xcb67c072733beefc,
        ]);
        a.div2();
        assert_eq!(
            a,
            ScalarRepr::new([
                0x5e949058ce4b9190,
                0x8ba76823351df1bf,
                0x6a346af1daa8f85a,
                0x65b3e039399df77e
            ])
        );
        for _ in 0..10 {
            a.div2();
        }
        assert_eq!(
            a,
            ScalarRepr::new([
                0x6fd7a524163392e4,
                0x16a2e9da08cd477c,
                0xdf9a8d1abc76aa3e,
                0x196cf80e4e677d
            ])
        );
        for _ in 0..200 {
            a.div2();
        }
        assert_eq!(a, ScalarRepr::new([0x196cf80e4e67, 0x0, 0x0, 0x0]));
        for _ in 0..40 {
            a.div2();
        }
        assert_eq!(a, ScalarRepr::new([0x19, 0x0, 0x0, 0x0]));
        for _ in 0..4 {
            a.div2();
        }
        assert_eq!(a, ScalarRepr::new([0x1, 0x0, 0x0, 0x0]));
        a.div2();
        assert!(a.is_zero());
    }

    #[test]
    fn test_scalar_repr_shr() {
        let mut a = ScalarRepr::new([
            0xb33fbaec482a283f,
            0x997de0d3a88cb3df,
            0x9af62d2a9a0e5525,
            0x36003ab08de70da1,
        ]);
        a.shr(0);
        assert_eq!(
            a,
            ScalarRepr::new([
                0xb33fbaec482a283f,
                0x997de0d3a88cb3df,
                0x9af62d2a9a0e5525,
                0x36003ab08de70da1
            ])
        );
        a.shr(1);
        assert_eq!(
            a,
            ScalarRepr::new([
                0xd99fdd762415141f,
                0xccbef069d44659ef,
                0xcd7b16954d072a92,
                0x1b001d5846f386d0
            ])
        );
        a.shr(50);
        assert_eq!(
            a,
            ScalarRepr::new([
                0xbc1a7511967bf667,
                0xc5a55341caa4b32f,
                0x75611bce1b4335e,
                0x6c0
            ])
        );
        a.shr(130);
        assert_eq!(a, ScalarRepr::new([0x1d5846f386d0cd7, 0x1b0, 0x0, 0x0]));
        a.shr(64);
        assert_eq!(a, ScalarRepr::new([0x1b0, 0x0, 0x0, 0x0]));
    }

    #[test]
    fn test_scalar_repr_mul2() {
        let mut a = ScalarRepr::from(23712937547);
        a.mul2();
        assert_eq!(a, ScalarRepr::new([0xb0acd6c96, 0x0, 0x0, 0x0]));
        for _ in 0..60 {
            a.mul2();
        }
        assert_eq!(
            a,
            ScalarRepr::new([0x6000000000000000, 0xb0acd6c9, 0x0, 0x0])
        );
        for _ in 0..128 {
            a.mul2();
        }
        assert_eq!(
            a,
            ScalarRepr::new([0x0, 0x0, 0x6000000000000000, 0xb0acd6c9])
        );
        for _ in 0..60 {
            a.mul2();
        }
        assert_eq!(a, ScalarRepr::new([0x0, 0x0, 0x0, 0x9600000000000000]));
        for _ in 0..7 {
            a.mul2();
        }
        assert!(a.is_zero());
    }

    #[test]
    fn test_scalar_repr_num_bits() {
        let mut a = ScalarRepr::from(0);
        assert_eq!(0, a.num_bits());
        a = ScalarRepr::from(1);
        for i in 1..257 {
            assert_eq!(i, a.num_bits());
            a.mul2();
        }
        assert_eq!(0, a.num_bits());
    }

    #[test]
    fn test_scalar_repr_sub_noborrow() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let mut t = ScalarRepr::new([
            0x8e62a7e85264e2c3,
            0xb23d34c1941d3ca,
            0x5976930b7502dd15,
            0x600f3fb517bf5495,
        ]);
        t.sub_noborrow(&ScalarRepr::new([
            0xd64f669809cbc6a4,
            0xfa76cb9d90cf7637,
            0xfefb0df9038d43b3,
            0x298a30c744b31acf,
        ]));
        assert!(
            t == ScalarRepr::new([
                0xb813415048991c1f,
                0x10ad07ae88725d92,
                0x5a7b851271759961,
                0x36850eedd30c39c5
            ])
        );

        for _ in 0..1000 {
            let mut a = Scalar::random(&mut rng).into_repr();
            a.0.l[3] >>= 30;
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

        // Subtracting r+1 from r should produce -1 (mod 2**256)
        let mut qplusone = ScalarRepr::new([
            0xffffffff00000001,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ]);
        qplusone.sub_noborrow(&ScalarRepr::new([
            0xffffffff00000002,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ]));
        assert_eq!(
            qplusone,
            ScalarRepr::new([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff
            ])
        );
    }

    #[test]
    fn test_scalar_legendre() {
        use fff::LegendreSymbol::*;
        use fff::SqrtField;

        assert_eq!(QuadraticResidue, Scalar::one().legendre());
        assert_eq!(Zero, Scalar::zero().legendre());

        let e = ScalarRepr::new([
            0x0dbc5349cd5664da,
            0x8ac5b6296e3ae29d,
            0x127cb819feceaa3b,
            0x3a6b21fb03867191,
        ]);
        assert_eq!(QuadraticResidue, Scalar::from_repr(e).unwrap().legendre());
        let e = ScalarRepr::new([
            0x96341aefd047c045,
            0x9b5f4254500a4d65,
            0x1ee08223b68ac240,
            0x31d9cd545c0ec7c6,
        ]);
        assert_eq!(
            QuadraticNonResidue,
            Scalar::from_repr(e).unwrap().legendre()
        );
    }

    #[test]
    fn test_scalar_repr_add_nocarry() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let mut t = ScalarRepr::new([
            0xd64f669809cbc6a4,
            0xfa76cb9d90cf7637,
            0xfefb0df9038d43b3,
            0x298a30c744b31acf,
        ]);
        t.add_nocarry(&ScalarRepr::new([
            0x8e62a7e85264e2c3,
            0xb23d34c1941d3ca,
            0x5976930b7502dd15,
            0x600f3fb517bf5495,
        ]));
        assert_eq!(
            t,
            ScalarRepr::new([
                0x64b20e805c30a967,
                0x59a9ee9aa114a02,
                0x5871a104789020c9,
                0x8999707c5c726f65
            ])
        );

        // Test for the associativity of addition.
        for _ in 0..1000 {
            let mut a = Scalar::random(&mut rng).into_repr();
            let mut b = Scalar::random(&mut rng).into_repr();
            let mut c = Scalar::random(&mut rng).into_repr();

            // Unset the first few bits, so that overflow won't occur.
            a.0.l[3] >>= 3;
            b.0.l[3] >>= 3;
            c.0.l[3] >>= 3;

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

        // Adding 1 to (2^256 - 1) should produce zero
        let mut x = ScalarRepr::new([
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
            0xffffffffffffffff,
        ]);
        x.add_nocarry(&ScalarRepr::from(1));
        assert!(x.is_zero());
    }

    #[test]
    fn test_scalar_add_assign() {
        {
            // Random number
            let mut tmp = Scalar(
                ScalarRepr::new([
                    0x437ce7616d580765,
                    0xd42d1ccb29d1235b,
                    0xed8f753821bd1423,
                    0x4eede1c9c89528ca,
                ])
                .0,
            );
            // assert!(tmp.is_valid());
            // Test that adding zero has no effect.
            tmp.add_assign(&Scalar(ScalarRepr::from(0).0));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ])
                    .0
                )
            );
            // Add one and test for the result.
            tmp.add_assign(&Scalar(ScalarRepr::from(1).0));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0x437ce7616d580766,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ])
                    .0
                )
            );
            // Add another random number that exercises the reduction.
            tmp.add_assign(&Scalar(
                ScalarRepr::new([
                    0x946f435944f7dc79,
                    0xb55e7ee6533a9b9b,
                    0x1e43b84c2f6194ca,
                    0x58717ab525463496,
                ])
                .0,
            ));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0xd7ec2abbb24fe3de,
                        0x35cdf7ae7d0d62f7,
                        0xd899557c477cd0e9,
                        0x3371b52bc43de018
                    ])
                    .0
                )
            );
            // Add one to (r - 1) and test for the result.
            tmp = Scalar(
                ScalarRepr::new([
                    0xffffffff00000000,
                    0x53bda402fffe5bfe,
                    0x3339d80809a1d805,
                    0x73eda753299d7d48,
                ])
                .0,
            );
            tmp.add_assign(&Scalar(ScalarRepr::from(1).0));
            assert!(tmp.into_repr().is_zero());
            // Add a random number to another one such that the result is r - 1
            tmp = Scalar(
                ScalarRepr::new([
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ])
                .0,
            );
            tmp.add_assign(&Scalar(
                ScalarRepr::new([
                    0x521a525223349e70,
                    0xa99bb5f3d8231f31,
                    0xde8e397bebe477e,
                    0x1ad08e5041d7c321,
                ])
                .0,
            ));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0xffffffff00000000,
                        0x53bda402fffe5bfe,
                        0x3339d80809a1d805,
                        0x73eda753299d7d48
                    ])
                    .0
                )
            );
            // Add one to the result and test for it.
            tmp.add_assign(&Scalar(ScalarRepr::from(1).0));
            assert!(tmp.into_repr().is_zero());
        }

        // Test associativity

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..1000 {
            // Generate a, b, c and ensure (a + b) + c == a + (b + c).
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);
            let c = Scalar::random(&mut rng);

            let mut tmp1 = a;
            tmp1.add_assign(&b);
            tmp1.add_assign(&c);

            let mut tmp2 = b;
            tmp2.add_assign(&c);
            tmp2.add_assign(&a);

            // assert!(tmp1.is_valid());
            // assert!(tmp2.is_valid());
            assert_eq!(tmp1, tmp2, "round {}", i);
        }
    }

    #[test]
    fn test_scalar_sub_assign() {
        {
            // Test arbitrary subtraction that tests reduction.
            let mut tmp = Scalar(
                ScalarRepr::new([
                    0x6a68c64b6f735a2b,
                    0xd5f4d143fe0a1972,
                    0x37c17f3829267c62,
                    0xa2f37391f30915c,
                ])
                .0,
            );
            tmp.sub_assign(&Scalar(
                ScalarRepr::new([
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ])
                .0,
            ));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0xbc83189d92a7f89c,
                        0x7f908737d62d38a3,
                        0x45aa62cfe7e4c3e1,
                        0x24ffc5896108547d
                    ])
                    .0
                )
            );

            // Test the opposite subtraction which doesn't test reduction.
            tmp = Scalar(
                ScalarRepr::new([
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ])
                .0,
            );
            tmp.sub_assign(&Scalar(
                ScalarRepr::new([
                    0x6a68c64b6f735a2b,
                    0xd5f4d143fe0a1972,
                    0x37c17f3829267c62,
                    0xa2f37391f30915c,
                ])
                .0,
            ));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ])
                    .0
                )
            );

            // Test for sensible results with zero
            tmp = Scalar(ScalarRepr::from(0).0);
            tmp.sub_assign(&Scalar(ScalarRepr::from(0).0));
            assert!(tmp.is_zero());

            tmp = Scalar(
                ScalarRepr::new([
                    0x437ce7616d580765,
                    0xd42d1ccb29d1235b,
                    0xed8f753821bd1423,
                    0x4eede1c9c89528ca,
                ])
                .0,
            );
            tmp.sub_assign(&Scalar(ScalarRepr::from(0).0));
            assert_eq!(
                tmp,
                Scalar(
                    ScalarRepr::new([
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ])
                    .0
                )
            );
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure that (a - b) + (b - a) = 0.
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);

            let mut tmp1 = a;
            tmp1.sub_assign(&b);

            let mut tmp2 = b;
            tmp2.sub_assign(&a);

            tmp1.add_assign(&tmp2);
            assert!(tmp1.is_zero());
        }
    }

    #[test]
    fn test_scalar_mul_assign() {
        let mut tmp = Scalar(
            ScalarRepr::new([
                0x6b7e9b8faeefc81a,
                0xe30a8463f348ba42,
                0xeff3cb67a8279c9c,
                0x3d303651bd7c774d,
            ])
            .0,
        );
        tmp.mul_assign(&Scalar(
            ScalarRepr::new([
                0x13ae28e3bc35ebeb,
                0xa10f4488075cae2c,
                0x8160e95a853c3b5d,
                0x5ae3f03b561a841d,
            ])
            .0,
        ));
        assert!(
            tmp == Scalar(
                ScalarRepr::new([
                    0x23717213ce710f71,
                    0xdbee1fe53a16e1af,
                    0xf565d3e1c2a48000,
                    0x4426507ee75df9d7
                ])
                .0
            )
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * b) * c = a * (b * c)
            let a = Scalar::random(&mut rng);
            let b = Scalar::random(&mut rng);
            let c = Scalar::random(&mut rng);

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

            let r = Scalar::random(&mut rng);
            let mut a = Scalar::random(&mut rng);
            let mut b = Scalar::random(&mut rng);
            let mut c = Scalar::random(&mut rng);

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
    fn test_scalar_squaring() {
        let mut a = Scalar(
            ScalarRepr::new([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x73eda753299d7d47,
            ])
            .0,
        );
        // assert!(a.is_valid());
        a.square();
        assert_eq!(
            a,
            Scalar::from_repr(ScalarRepr::new([
                0xc0d698e7bde077b8,
                0xb79a310579e76ec2,
                0xac1da8d0a9af4e5f,
                0x13f629c49bf23e97
            ]))
            .unwrap()
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * a) = a^2
            let a = Scalar::random(&mut rng);

            let mut tmp = a;
            tmp.square();

            let mut tmp2 = a;
            tmp2.mul_assign(&a);

            assert_eq!(tmp, tmp2);
        }
    }

    #[test]
    fn test_scalar_inverse() {
        assert!(Scalar::zero().inverse().is_none());

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let one = Scalar::one();

        for i in 0..1000 {
            // Ensure that a * a^-1 = 1
            let mut a = Scalar::random(&mut rng);
            let ainv = a.inverse().unwrap();
            a.mul_assign(&ainv);
            assert_eq!(a, one, "round {}", i);
        }
    }

    #[test]
    fn test_scalar_double() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure doubling a is equivalent to adding a to itself.
            let mut a = Scalar::random(&mut rng);
            let mut b = a;
            b.add_assign(&a);
            a.double();
            assert_eq!(a, b);
        }
    }

    #[test]
    fn test_scalar_negate() {
        {
            let mut a = Scalar::zero();
            a.negate();

            assert!(a.is_zero());
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure (a - (-a)) = 0.
            let mut a = Scalar::random(&mut rng);
            let mut b = a;
            b.negate();
            a.add_assign(&b);

            assert!(a.is_zero());
        }
    }

    #[test]
    fn test_scalar_pow() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..1000 {
            // Exponentiate by various small numbers and ensure it consists with repeated
            // multiplication.
            let a = Scalar::random(&mut rng);
            let target = a.pow(&[i]);
            let mut c = Scalar::one();
            for _ in 0..i {
                c.mul_assign(&a);
            }
            assert_eq!(c, target);
        }

        for _ in 0..1000 {
            // Exponentiating by the modulus should have no effect in a prime field.
            let a = Scalar::random(&mut rng);

            assert_eq!(a, a.pow(Scalar::char()));
        }
    }

    #[test]
    fn test_scalar_sqrt() {
        use fff::SqrtField;

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());

        for _ in 0..1000 {
            // Ensure sqrt(a^2) = a or -a
            let a = Scalar::random(&mut rng);
            let mut nega = a;
            nega.negate();
            let mut b = a;
            b.square();

            let b = b.sqrt().unwrap();

            assert!(a == b || nega == b);
        }

        for _ in 0..1000 {
            // Ensure sqrt(a)^2 = a for random a
            let a = Scalar::random(&mut rng);

            if let Some(mut tmp) = a.sqrt() {
                tmp.square();

                assert_eq!(a, tmp);
            }
        }
    }

    #[test]
    fn test_scalar_from_into_repr() {
        // r + 1 should not be in the field
        assert!(Scalar::from_repr(ScalarRepr::new([
            0xffffffff00000002,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48
        ]))
        .is_err());

        // r should not be in the field
        assert!(Scalar::from_repr(Scalar::char()).is_err());

        // Multiply some arbitrary representations to see if the result is as expected.
        let a = ScalarRepr::new([
            0x25ebe3a3ad3c0c6a,
            0x6990e39d092e817c,
            0x941f900d42f5658e,
            0x44f8a103b38a71e0,
        ]);
        let mut a_fr = Scalar::from_repr(a).unwrap();
        let b = ScalarRepr::new([
            0x264e9454885e2475,
            0x46f7746bb0308370,
            0x4683ef5347411f9,
            0x58838d7f208d4492,
        ]);
        let b_fr = Scalar::from_repr(b).unwrap();
        let c = ScalarRepr::new([
            0x48a09ab93cfc740d,
            0x3a6600fbfc7a671,
            0x838567017501d767,
            0x7161d6da77745512,
        ]);
        a_fr.mul_assign(&b_fr);
        assert_eq!(a_fr.into_repr(), c);

        // Zero should be in the field.
        assert!(Scalar::from_repr(ScalarRepr::from(0)).unwrap().is_zero());

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..1000 {
            // Try to turn Scalar elements into representations and back again, and compare.
            let a = Scalar::random(&mut rng);
            let a_repr = a.into_repr();
            let b_repr = ScalarRepr::from(a);
            assert_eq!(a_repr, b_repr);
            let a_again = Scalar::from_repr(a_repr).unwrap();

            assert_eq!(a, a_again, "{}", i);
        }
    }

    #[test]
    fn test_scalar_repr_display() {
        assert_eq!(
            format!(
                "{}",
                ScalarRepr::new([
                    0x2829c242fa826143,
                    0x1f32cf4dd4330917,
                    0x932e4e479d168cd9,
                    0x513c77587f563f64
                ])
            ),
            "0x513c77587f563f64932e4e479d168cd91f32cf4dd43309172829c242fa826143".to_string()
        );
        assert_eq!(
            format!(
                "{}",
                ScalarRepr::new([
                    0x25ebe3a3ad3c0c6a,
                    0x6990e39d092e817c,
                    0x941f900d42f5658e,
                    0x44f8a103b38a71e0
                ])
            ),
            "0x44f8a103b38a71e0941f900d42f5658e6990e39d092e817c25ebe3a3ad3c0c6a".to_string()
        );
        assert_eq!(
            format!(
                "{}",
                ScalarRepr::new([
                    0xffffffffffffffff,
                    0xffffffffffffffff,
                    0xffffffffffffffff,
                    0xffffffffffffffff
                ])
            ),
            "0xffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff".to_string()
        );
        assert_eq!(
            format!("{}", ScalarRepr::new([0, 0, 0, 0])),
            "0x0000000000000000000000000000000000000000000000000000000000000000".to_string()
        );
    }

    #[test]
    fn test_scalar_display() {
        assert_eq!(
            format!(
                "{}",
                Scalar::from_repr(ScalarRepr::new([
                    0xc3cae746a3b5ecc7,
                    0x185ec8eb3f5b5aee,
                    0x684499ffe4b9dd99,
                    0x7c9bba7afb68faa
                ]))
                .unwrap()
            ),
            "Scalar(0x07c9bba7afb68faa684499ffe4b9dd99185ec8eb3f5b5aeec3cae746a3b5ecc7)"
                .to_string()
        );
        assert_eq!(
            format!(
                "{}",
                Scalar::from_repr(ScalarRepr::new([
                    0x44c71298ff198106,
                    0xb0ad10817df79b6a,
                    0xd034a80a2b74132b,
                    0x41cf9a1336f50719
                ]))
                .unwrap()
            ),
            "Scalar(0x41cf9a1336f50719d034a80a2b74132bb0ad10817df79b6a44c71298ff198106)"
                .to_string()
        );
    }

    #[test]
    fn test_scalar_num_bits() {
        assert_eq!(Scalar::NUM_BITS, 255);
        assert_eq!(Scalar::CAPACITY, 254);
    }

    #[test]
    fn test_scalar_root_of_unity() {
        use fff::SqrtField;

        assert_eq!(Scalar::S, 32);
        assert_eq!(
            Scalar::multiplicative_generator(),
            Scalar::from_repr(ScalarRepr::from(7)).unwrap()
        );
        assert_eq!(
            Scalar::multiplicative_generator().pow([
                0xfffe5bfeffffffff,
                0x9a1d80553bda402,
                0x299d7d483339d808,
                0x73eda753
            ]),
            Scalar::root_of_unity()
        );
        assert_eq!(Scalar::root_of_unity().pow([1 << Scalar::S]), Scalar::one());
        assert!(Scalar::multiplicative_generator().sqrt().is_none());
    }

    #[test]
    fn scalar_field_tests() {
        crate::tests::field::random_field_tests::<Scalar>();
        crate::tests::field::random_sqrt_tests::<Scalar>();
        crate::tests::field::random_frobenius_tests::<Scalar, _>(Scalar::char(), 13);
        crate::tests::field::from_str_tests::<Scalar>();
    }

    #[test]
    fn scalar_repr_tests() {
        crate::tests::repr::random_repr_tests::<Scalar>();
    }

    #[test]
    fn test_modulus() {
        assert_eq!(
            format!("{:?}", Scalar::char()),
            "0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001"
        );
    }

    #[test]
    fn test_scalar_repr_conversion() {
        let a = ScalarRepr::from(1);
        let b = ScalarRepr(blst::blst_fr { l: [1, 0, 0, 0] });
        assert_eq!(a, b);

        let a = Scalar::from(1);
        let b = Scalar::from_repr(ScalarRepr::from(1)).unwrap();
        assert_eq!(a, b);
        assert_eq!(Scalar::from(1).into_repr(), ScalarRepr::from(1));

        let a = Scalar::from(12);
        assert_eq!(a, Scalar::from_repr(a.into_repr()).unwrap());

        let a = ScalarRepr::from(12);
        assert_eq!(Scalar::from_repr(a).unwrap().into_repr(), a);
    }
}
