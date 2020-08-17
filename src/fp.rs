//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`

use blst::*;

use core::{
    convert::TryInto,
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use fff::Field;

/// `Fp` values are always in
/// Montgomery form; i.e., Scalar(a) = aR mod p, with R = 2^384.
#[derive(Copy, Clone)]
pub struct Fp(pub(crate) blst_fp);

/// Representation of a `Fp`, in regular coordinates.
#[derive(Default, Clone, Copy)]
pub struct FpRepr(blst_fp);

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
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for FpRepr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "0x")?;
        for &b in self.0.l.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl From<u32> for FpRepr {
    fn from(val: u32) -> FpRepr {
        let mut raw = blst_fp::default();

        unsafe { blst_fp_from_uint32(&mut raw as *mut _, val as *const _) };

        FpRepr(raw)
    }
}

impl From<u64> for FpRepr {
    fn from(val: u64) -> FpRepr {
        let mut raw = blst_fp::default();

        unsafe { blst_fp_from_uint64(&mut raw as *mut _, val as *const _) };

        FpRepr(raw)
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
            *self = Self::from(0u32);
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
            *self = Self::from(0u32);
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

impl fmt::Debug for Fp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes_le();
        write!(f, "0x")?;
        for &b in tmp.iter() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Display for Fp {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes_le();
        write!(f, "0x")?;
        for &b in tmp.iter() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
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

impl TryInto<Fp> for blst_fp {
    type Error = NotInFieldError;

    fn try_into(self) -> Result<Fp, Self::Error> {
        let fp = Fp(self);

        if !fp.is_valid() {
            return Err(NotInFieldError);
        }

        Ok(fp)
    }
}

impl From<Fp> for FpRepr {
    fn from(val: Fp) -> Self {
        let raw: blst_fp = val.into();
        FpRepr(raw)
    }
}

impl From<FpRepr> for Fp {
    fn from(val: FpRepr) -> Self {
        let mut raw = blst_fp::default();
        unsafe { blst_fp_from(&mut raw, &val.0) };
        Fp(raw)
    }
}

impl From<Fp> for blst_fp {
    fn from(val: Fp) -> blst_fp {
        val.0
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
        use fff::PrimeField;

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
        *self *= other;
    }

    fn inverse(&self) -> Option<Self> {
        todo!()
    }

    fn frobenius_map(&mut self, power: usize) {
        todo!()
    }
}

impl fff::PrimeField for Fp {
    type Repr = FpRepr;

    const NUM_BITS: u32 = 381;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const S: u32 = 1;

    fn from_repr(repr: Self::Repr) -> Result<Self, fff::PrimeFieldDecodingError> {
        repr.0.try_into().map_err(|err: NotInFieldError| {
            fff::PrimeFieldDecodingError::NotInField(err.to_string())
        })
    }

    /// Convert a biginteger representation into a prime field element, if
    /// the number is an element of the field.
    fn into_repr(&self) -> Self::Repr {
        (*self).into()
    }

    fn char() -> Self::Repr {
        FpRepr(blst_fp {
            l: [
                13402431016077863595,
                2210141511517208575,
                7435674573564081700,
                7239337960414712511,
                5412103778470702295,
                1873798617647539866,
            ],
        })
        .into()
    }

    fn multiplicative_generator() -> Self {
        FpRepr(blst_fp {
            l: [2, 0, 0, 0, 0, 0],
        })
        .into()
    }

    fn root_of_unity() -> Self {
        FpRepr(blst_fp {
            l: [
                13402431016077863594,
                2210141511517208575,
                7435674573564081700,
                7239337960414712511,
                5412103778470702295,
                1873798617647539866,
            ],
        })
        .into()
    }
}

impl fff::SqrtField for Fp {
    fn legendre(&self) -> fff::LegendreSymbol {
        todo!()
    }

    fn sqrt(&self) -> Option<Self> {
        todo!()
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

    /// Computes the multiplicative inverse of this field
    /// element, returning None in the case that this element
    /// is zero.
    pub fn invert(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }

        let mut out = blst_fp::default();

        unsafe { blst_fp_eucl_inverse(&mut out, &self.0) };

        Some(Fp(out))
    }

    /// Computes a uniformly random element using rejection sampling.
    pub fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        // The number of bits we should "shave" from a randomly sampled reputation.
        const REPR_SHAVE_BITS: usize = 384 - 381;

        loop {
            let mut raw = Fp::default();
            for i in 0..6 {
                raw.0.l[i] = rng.next_u64();
            }

            // Mask away the unused most-significant bits.
            raw.0.l[5] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if raw.is_valid() {
                return raw;
            }
        }
    }

    /// Constructs an element of `Fp` without checking that it is canonical.
    pub fn from_raw_unchecked(v: [u64; 6]) -> Fp {
        let mut inner = blst_fp::default();
        inner.l.copy_from_slice(&v);
        Fp(inner)
    }

    pub fn is_valid(&self) -> bool {
        todo!()
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

    /// Squares this element.
    pub fn square(&self) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_sqr(&mut out, &self.0) };

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
