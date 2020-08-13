//! This module provides an implementation of the BLS12-381 base field `GF(p)`
//! where `p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab`

use blst::*;

use core::{
    convert::TryInto,
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// `Fp` values are always in
/// Montgomery form; i.e., Scalar(a) = aR mod p, with R = 2^384.
#[derive(Copy, Clone)]
pub struct Fp(pub(crate) blst_fp);

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

impl Fp {
    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Fp {
        Fp(blst_fp {
            l: [0, 0, 0, 0, 0, 0],
        })
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Fp {
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

    pub fn is_zero(&self) -> bool {
        self == &Fp::zero()
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into an `Fp`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 48]) -> Option<Fp> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_fp::default();

        unsafe {
            blst_fp_from_lendian(&mut raw as _, in_v.as_mut_ptr());
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
            blst_fp_from_bendian(&mut raw as _, in_v.as_mut_ptr());
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

    pub fn pow(&self, by: &[u64; 6]) -> Self {
        todo!()
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

        unsafe { blst_fp_mul_by_3(&mut out as _, &self.0 as _) };

        Fp(out)
    }

    /// Multiplies `self` with `8`, returning the result.
    pub fn mul8(&self) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_mul_by_8(&mut out as _, &self.0 as _) };

        Fp(out)
    }

    /// Left shift `self` by `count`, returning the result.
    pub fn shl(&self, count: usize) -> Self {
        let mut out = blst_fp::default();

        unsafe { blst_fp_lshift(&mut out as _, &self.0 as _, count) };

        Fp(out)
    }
}
