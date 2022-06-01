//! An implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::{
    borrow::Borrow,
    cmp,
    convert::TryInto,
    fmt,
    iter::{Product, Sum},
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use byte_slice_cast::AsByteSlice;
use ff::{Field, FieldBits, PrimeField, PrimeFieldBits};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
///
/// The inner representation `blst_fr` is stored in Montgomery form as little-endian `u64` limbs.
#[derive(Default, Clone, Copy)]
#[repr(transparent)]
pub struct Scalar(pub(crate) blst_fr);

// GENERATOR = 7 (multiplicative generator of r-1 order, that is also quadratic nonresidue)
const GENERATOR: Scalar = Scalar(blst_fr {
    l: [
        0x0000_000e_ffff_fff1,
        0x17e3_63d3_0018_9c0f,
        0xff9c_5787_6f84_57b0,
        0x3513_3220_8fc5_a8c4,
    ],
});

// Little-endian non-Montgomery form not reduced mod p.
#[allow(dead_code)]
const MODULUS: [u64; 4] = [
    0xffff_ffff_0000_0001,
    0x53bd_a402_fffe_5bfe,
    0x3339_d808_09a1_d805,
    0x73ed_a753_299d_7d48,
];

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 8] = [
    0x0000_0001,
    0xffff_ffff,
    0xfffe_5bfe,
    0x53bd_a402,
    0x09a1_d805,
    0x3339_d808,
    0x299d_7d48,
    0x73ed_a753,
];

// Little-endian non-Montgomery form not reduced mod p.
const MODULUS_REPR: [u8; 32] = [
    0x01, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0xfe, 0x5b, 0xfe, 0xff, 0x02, 0xa4, 0xbd, 0x53,
    0x05, 0xd8, 0xa1, 0x09, 0x08, 0xd8, 0x39, 0x33, 0x48, 0x7d, 0x9d, 0x29, 0x53, 0xa7, 0xed, 0x73,
];

// `2^S` root of unity in little-endian Montgomery form.
const ROOT_OF_UNITY: Scalar = Scalar(blst_fr {
    l: [
        0xb9b5_8d8c_5f0e_466a,
        0x5b1b_4c80_1819_d7ec,
        0x0af5_3ae3_52a3_1e64,
        0x5bf3_adda_19e9_b27b,
    ],
});

const ZERO: Scalar = Scalar(blst_fr { l: [0, 0, 0, 0] });

/// `R = 2^256 mod q` in little-endian Montgomery form which is equivalent to 1 in little-endian
/// non-Montgomery form.
///
/// sage> mod(2^256, 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
/// sage> 0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe
const R: Scalar = Scalar(blst_fr {
    l: [
        0x0000_0001_ffff_fffe,
        0x5884_b7fa_0003_4802,
        0x998c_4fef_ecbc_4ff5,
        0x1824_b159_acc5_056f,
    ],
});

/// `R^2 = 2^512 mod q` in little-endian Montgomery form.
///
/// sage> mod(2^512, 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001)
/// sage> 0x748d9d99f59ff1105d314967254398f2b6cedcb87925c23c999e990f3f29c6d
#[allow(dead_code)]
const R2: Scalar = Scalar(blst_fr {
    l: [
        0xc999_e990_f3f2_9c6d,
        0x2b6c_edcb_8792_5c23,
        0x05d3_1496_7254_398f,
        0x0748_d9d9_9f59_ff11,
    ],
});

pub const S: u32 = 32;

impl fmt::Debug for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let be_bytes = self.to_bytes_be();
        write!(f, "Scalar(0x")?;
        for &b in be_bytes.iter() {
            write!(f, "{:02x}", b)?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Ord for Scalar {
    #[allow(clippy::comparison_chain)]
    fn cmp(&self, other: &Scalar) -> cmp::Ordering {
        for (a, b) in self.to_bytes_be().iter().zip(other.to_bytes_be().iter()) {
            if a > b {
                return cmp::Ordering::Greater;
            } else if a < b {
                return cmp::Ordering::Less;
            }
        }
        cmp::Ordering::Equal
    }
}

impl PartialOrd for Scalar {
    #[inline(always)]
    fn partial_cmp(&self, other: &Scalar) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l == other.0.l
    }
}

impl Eq for Scalar {}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0.l[0].ct_eq(&other.0.l[0])
            & self.0.l[1].ct_eq(&other.0.l[1])
            & self.0.l[2].ct_eq(&other.0.l[2])
            & self.0.l[3].ct_eq(&other.0.l[3])
    }
}

impl ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Scalar(blst_fr {
            l: [
                u64::conditional_select(&a.0.l[0], &b.0.l[0], choice),
                u64::conditional_select(&a.0.l[1], &b.0.l[1], choice),
                u64::conditional_select(&a.0.l[2], &b.0.l[2], choice),
                u64::conditional_select(&a.0.l[3], &b.0.l[3], choice),
            ],
        })
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
        let mut repr = [0u8; 32];
        repr[..8].copy_from_slice(&val.to_le_bytes());
        Scalar::from_bytes_le(&repr).unwrap()
    }
}

#[allow(clippy::from_over_into)]
impl Into<blst_scalar> for Scalar {
    fn into(self) -> blst_scalar {
        let mut out = blst_scalar::default();
        unsafe {
            blst_scalar_from_fr(&mut out, &self.0);
        }

        out
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

        let mut out = blst_fr::default();

        unsafe { blst_fr_from_scalar(&mut out, &self) };

        Ok(Scalar(out))
    }
}

impl Neg for &Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        let mut neg = *self;
        unsafe { blst_fr_cneg(&mut neg.0, &self.0, true) };
        neg
    }
}

impl Neg for Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        -&self
    }
}

impl Add<&Scalar> for &Scalar {
    type Output = Scalar;

    #[inline]
    fn add(self, rhs: &Scalar) -> Scalar {
        let mut out = *self;
        out += rhs;
        out
    }
}

impl Sub<&Scalar> for &Scalar {
    type Output = Scalar;

    #[inline]
    fn sub(self, rhs: &Scalar) -> Scalar {
        let mut out = *self;
        out -= rhs;
        out
    }
}

impl Mul<&Scalar> for &Scalar {
    type Output = Scalar;

    #[inline]
    fn mul(self, rhs: &Scalar) -> Scalar {
        let mut out = *self;
        out *= rhs;
        out
    }
}

impl AddAssign<&Scalar> for Scalar {
    #[inline]
    fn add_assign(&mut self, rhs: &Scalar) {
        unsafe { blst_fr_add(&mut self.0, &self.0, &rhs.0) };
    }
}

impl SubAssign<&Scalar> for Scalar {
    #[inline]
    fn sub_assign(&mut self, rhs: &Scalar) {
        unsafe { blst_fr_sub(&mut self.0, &self.0, &rhs.0) };
    }
}

impl MulAssign<&Scalar> for Scalar {
    #[inline]
    fn mul_assign(&mut self, rhs: &Scalar) {
        unsafe { blst_fr_mul(&mut self.0, &self.0, &rhs.0) };
    }
}

impl<T> Sum<T> for Scalar
where
    T: Borrow<Scalar>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Scalar::zero(), |sum, x| sum + x.borrow())
    }
}

impl<T> Product<T> for Scalar
where
    T: Borrow<Scalar>,
{
    fn product<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Scalar::one(), |product, x| product * x.borrow())
    }
}

impl_add_sub!(Scalar);
impl_add_sub_assign!(Scalar);
impl_mul!(Scalar);
impl_mul_assign!(Scalar);

/// The number of bits we should "shave" from a randomly sampled reputation.
const REPR_SHAVE_BITS: usize = 256 - Scalar::NUM_BITS as usize;

impl Field for Scalar {
    fn random(mut rng: impl RngCore) -> Self {
        loop {
            let mut raw = [0u64; 4];
            for int in raw.iter_mut() {
                *int = rng.next_u64();
            }

            // Mask away the unused most-significant bits.
            raw[3] &= 0xffffffffffffffff >> REPR_SHAVE_BITS;

            if let Some(scalar) = Scalar::from_u64s_le(&raw).into() {
                return scalar;
            }
        }
    }

    fn zero() -> Self {
        ZERO
    }

    fn one() -> Self {
        R
    }

    fn is_zero(&self) -> Choice {
        self.ct_eq(&ZERO)
    }

    fn square(&self) -> Self {
        let mut out = *self;
        out.square_assign();
        out
    }

    fn double(&self) -> Self {
        let mut out = *self;
        out += self;
        out
    }

    fn invert(&self) -> CtOption<Self> {
        let mut inv = blst_fr::default();
        unsafe { blst_fr_eucl_inverse(&mut inv, &self.0) };
        let is_invertible = !self.ct_eq(&Scalar::zero());
        CtOption::new(Scalar(inv), is_invertible)
    }

    fn sqrt(&self) -> CtOption<Self> {
        // Tonelli-Shank's algorithm for q mod 16 = 1
        // https://eprint.iacr.org/2012/685.pdf (page 12, algorithm 5)

        // w = self^((t - 1) // 2)
        //   = self^6104339283789297388802252303364915521546564123189034618274734669823
        let w = self.pow_vartime(&[
            0x7fff_2dff_7fff_ffff,
            0x04d0_ec02_a9de_d201,
            0x94ce_bea4_199c_ec04,
            0x0000_0000_39f6_d3a9,
        ]);

        let mut v = S;
        let mut x = self * w;
        let mut b = x * w;

        // Initialize z as the 2^S root of unity.
        let mut z = ROOT_OF_UNITY;

        for max_v in (1..=S).rev() {
            let mut k = 1;
            let mut tmp = b.square();
            let mut j_less_than_v: Choice = 1.into();

            for j in 2..max_v {
                let tmp_is_one = tmp.ct_eq(&Scalar::one());
                let squared = Scalar::conditional_select(&tmp, &z, tmp_is_one).square();
                tmp = Scalar::conditional_select(&squared, &tmp, tmp_is_one);
                let new_z = Scalar::conditional_select(&z, &squared, tmp_is_one);
                j_less_than_v &= !j.ct_eq(&v);
                k = u32::conditional_select(&j, &k, tmp_is_one);
                z = Scalar::conditional_select(&z, &new_z, j_less_than_v);
            }

            let result = x * z;
            x = Scalar::conditional_select(&result, &x, b.ct_eq(&Scalar::one()));
            z = z.square();
            b *= z;
            v = k;
        }

        CtOption::new(
            x,
            (x * x).ct_eq(self), // Only return Some if it's the square root.
        )
    }
}

/// Checks if the passed in bytes are less than the MODULUS. (both in non-Montgomery form and little endian).
/// Assumes that `a` is exactly 4 elements long.
#[allow(clippy::comparison_chain)]
fn is_valid(a: &[u64]) -> bool {
    debug_assert_eq!(a.len(), 4);
    for (a, b) in a.iter().zip(MODULUS.iter()).rev() {
        if a > b {
            return false;
        } else if a < b {
            return true;
        }
    }

    false
}

#[inline]
fn u64s_from_bytes(bytes: &[u8; 32]) -> [u64; 4] {
    [
        u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
        u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
        u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
        u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
    ]
}

impl PrimeField for Scalar {
    // Little-endian non-Montgomery form bigint mod p.
    type Repr = [u8; 32];

    const NUM_BITS: u32 = 255;
    const CAPACITY: u32 = Self::NUM_BITS - 1;
    const S: u32 = S;

    /// Converts a little-endian non-Montgomery form `repr` into a Montgomery form `Scalar`.
    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        Self::from_bytes_le(&repr)
    }

    fn from_repr_vartime(repr: Self::Repr) -> Option<Self> {
        let bytes_u64 = u64s_from_bytes(&repr);

        if !is_valid(&bytes_u64) {
            return None;
        }
        let mut out = blst_fr::default();
        unsafe { blst_fr_from_uint64(&mut out, bytes_u64.as_ptr()) };
        Some(Scalar(out))
    }

    /// Converts a Montgomery form `Scalar` into little-endian non-Montgomery from.
    fn to_repr(&self) -> Self::Repr {
        self.to_bytes_le()
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr()[0] & 1)
    }

    fn multiplicative_generator() -> Self {
        GENERATOR
    }

    fn root_of_unity() -> Self {
        ROOT_OF_UNITY
    }
}

#[cfg(not(target_pointer_width = "64"))]
type ReprBits = [u32; 8];

#[cfg(target_pointer_width = "64")]
type ReprBits = [u64; 4];

impl PrimeFieldBits for Scalar {
    // Representation in non-Montgomery form.
    type ReprBits = ReprBits;

    #[cfg(target_pointer_width = "64")]
    fn to_le_bits(&self) -> FieldBits<Self::ReprBits> {
        let mut limbs = [0u64; 4];
        unsafe { blst_uint64_from_fr(limbs.as_mut_ptr(), &self.0) };

        FieldBits::new(limbs)
    }

    #[cfg(not(target_pointer_width = "64"))]
    fn to_le_bits(&self) -> FieldBits<Self::ReprBits> {
        let bytes = self.to_bytes_le();
        let limbs = [
            u32::from_le_bytes(bytes[0..4].try_into().unwrap()),
            u32::from_le_bytes(bytes[4..8].try_into().unwrap()),
            u32::from_le_bytes(bytes[8..12].try_into().unwrap()),
            u32::from_le_bytes(bytes[12..16].try_into().unwrap()),
            u32::from_le_bytes(bytes[16..20].try_into().unwrap()),
            u32::from_le_bytes(bytes[20..24].try_into().unwrap()),
            u32::from_le_bytes(bytes[24..28].try_into().unwrap()),
            u32::from_le_bytes(bytes[28..32].try_into().unwrap()),
        ];
        FieldBits::new(limbs)
    }

    fn char_le_bits() -> FieldBits<Self::ReprBits> {
        #[cfg(not(target_pointer_width = "64"))]
        {
            FieldBits::new(MODULUS_LIMBS_32)
        }

        #[cfg(target_pointer_width = "64")]
        FieldBits::new(MODULUS)
    }
}

impl Scalar {
    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 32]) -> CtOption<Scalar> {
        let is_some =
            Choice::from(unsafe { blst_scalar_fr_check(&blst_scalar { b: *bytes }) as u8 });

        let mut out = blst_fr::default();
        let bytes_u64 = u64s_from_bytes(bytes);

        unsafe { blst_fr_from_uint64(&mut out, bytes_u64.as_ptr()) };

        CtOption::new(Scalar(out), is_some)
    }

    /// Attempts to convert a big-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_be(be_bytes: &[u8; 32]) -> CtOption<Scalar> {
        let mut le_bytes = *be_bytes;
        le_bytes.reverse();
        Self::from_bytes_le(&le_bytes)
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    #[inline]
    pub fn to_bytes_le(&self) -> [u8; 32] {
        let mut out = [0u64; 4];
        unsafe { blst_uint64_from_fr(out.as_mut_ptr(), &self.0) };
        out.as_byte_slice().try_into().unwrap()
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes_be(&self) -> [u8; 32] {
        let mut bytes = self.to_bytes_le();
        bytes.reverse();
        bytes
    }

    // `u64s` represent a little-endian non-Montgomery form integer mod p.
    pub fn from_u64s_le(bytes: &[u64; 4]) -> CtOption<Self> {
        let mut raw = blst_scalar::default();
        let mut out = blst_fr::default();

        unsafe { blst_scalar_from_uint64(&mut raw, bytes.as_ptr()) };
        let is_some = Choice::from(unsafe { blst_scalar_fr_check(&raw) as u8 });
        unsafe { blst_fr_from_scalar(&mut out, &raw) };

        CtOption::new(Scalar(out), is_some)
    }

    #[allow(clippy::match_like_matches_macro)]
    pub fn is_quad_res(&self) -> Choice {
        match self.legendre() {
            0 | 1 => Choice::from(1u8),
            _ => Choice::from(0u8),
        }
    }

    pub fn legendre(&self) -> i8 {
        const MOD_MINUS_1_OVER_2: [u64; 4] = [
            0x7fffffff80000000,
            0xa9ded2017fff2dff,
            0x199cec0404d0ec02,
            0x39f6d3a994cebea4,
        ];
        // s = self^((modulus - 1) // 2)
        let s = self.pow_vartime(MOD_MINUS_1_OVER_2);
        if s == Self::zero() {
            0
        } else if s == Self::one() {
            1
        } else {
            -1
        }
    }

    pub fn char() -> <Self as PrimeField>::Repr {
        MODULUS_REPR
    }

    pub fn num_bits(&self) -> u32 {
        let mut ret = 256;
        for i in self.to_bytes_be().iter() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != 8 {
                break;
            }
        }

        ret
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

    /// Calculates the `square` of this element.
    #[inline]
    pub fn square_assign(&mut self) {
        unsafe { blst_fr_sqr(&mut self.0, &self.0) };
    }
}

#[cfg(feature = "gpu")]
impl ec_gpu::GpuName for Scalar {
    fn name() -> String {
        ec_gpu::name!()
    }
}

#[cfg(feature = "gpu")]
impl ec_gpu::GpuField for Scalar {
    fn one() -> Vec<u32> {
        crate::u64_to_u32(&R.0.l[..])
    }

    fn r2() -> Vec<u32> {
        crate::u64_to_u32(&R2.0.l[..])
    }

    fn modulus() -> Vec<u32> {
        crate::u64_to_u32(&MODULUS[..])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    /// INV = -(q^{-1} mod 2^64) mod 2^64
    const INV: u64 = 0xfffffffeffffffff;

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
            inv = inv.wrapping_mul(MODULUS[0]);
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

        assert_ne!(Scalar::zero(), Scalar::one());
        assert_ne!(Scalar::one(), R2);
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
            R2,
        );

        // -1 should work
        assert!(bool::from(
            Scalar::from_bytes_le(&[
                0, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_some()
        ));

        // modulus is invalid
        assert!(bool::from(Scalar::from_bytes_le(&MODULUS_REPR).is_none()));

        // Anything larger than the modulus is invalid
        assert!(bool::from(
            Scalar::from_bytes_le(&[
                2, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes_le(&[
                1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 58, 51, 72, 125, 157, 41, 83, 167, 237, 115
            ])
            .is_none()
        ));
        assert!(bool::from(
            Scalar::from_bytes_le(&[
                1, 0, 0, 0, 255, 255, 255, 255, 254, 91, 254, 255, 2, 164, 189, 83, 5, 216, 161, 9,
                8, 216, 57, 51, 72, 125, 157, 41, 83, 167, 237, 116
            ])
            .is_none()
        ));
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
            Scalar(blst::blst_fr {
                l: [
                    0xfffffffeffffffff,
                    0x53bda402fffe5bfe,
                    0x3339d80809a1d805,
                    0x73eda753299d7d48
                ]
            })
        );

        let mut tmp = LARGEST;
        tmp += &Scalar(blst::blst_fr { l: [1, 0, 0, 0] });

        assert_eq!(tmp, Scalar::zero());
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar(blst::blst_fr { l: [1, 0, 0, 0] }));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar(blst::blst_fr { l: [1, 0, 0, 0] });
        assert_eq!(tmp, LARGEST);

        {
            let mut a = Scalar::zero();
            a = -a;

            assert!(bool::from(a.is_zero()));
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

            assert!(bool::from(a.is_zero()));
        }
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Scalar::zero());

        let mut tmp = Scalar::zero();
        tmp -= &LARGEST;

        let mut tmp2 = Scalar(blst::blst_fr { l: MODULUS });
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_multiplication() {
        let mut tmp = Scalar(blst::blst_fr {
            l: [
                0x6b7e9b8faeefc81a,
                0xe30a8463f348ba42,
                0xeff3cb67a8279c9c,
                0x3d303651bd7c774d,
            ],
        });
        tmp *= &Scalar(blst::blst_fr {
            l: [
                0x13ae28e3bc35ebeb,
                0xa10f4488075cae2c,
                0x8160e95a853c3b5d,
                0x5ae3f03b561a841d,
            ],
        });
        assert!(
            tmp == Scalar(blst::blst_fr {
                l: [
                    0x23717213ce710f71,
                    0xdbee1fe53a16e1af,
                    0xf565d3e1c2a48000,
                    0x4426507ee75df9d7
                ]
            })
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
        let mut r2 = r1;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r2.pow_vartime(&q_minus_2);

            assert_eq!(r1, r2);
            // Add R so we check something different next time around
            r1 += R;
            r2 = r1;
        }
    }

    #[test]
    fn test_sqrt() {
        {
            assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        }

        let mut square = Scalar(blst::blst_fr {
            l: [
                0x46cd85a5f273077e,
                0x1d30c47dd68fc735,
                0x77f656f60beca0eb,
                0x494aa01bdf32468d,
            ],
        });

        let mut none_count = 0;

        for _ in 0..100 {
            let square_root = square.sqrt();
            if square_root.is_none().into() {
                none_count += 1;
            } else {
                assert_eq!(square_root.unwrap() * square_root.unwrap(), square);
            }
            square -= Scalar::one();
        }

        assert_eq!(49, none_count);
    }

    #[test]
    fn test_double() {
        let a = Scalar::from_u64s_le(&[
            0x1fff3231233ffffd,
            0x4884b7fa00034802,
            0x998c4fefecbc4ff3,
            0x1824b159acc50562,
        ])
        .unwrap();

        assert_eq!(a.double(), a + a);
    }

    #[test]
    fn test_scalar_ordering() {
        fn assert_equality(a: Scalar, b: Scalar) {
            assert_eq!(a, b);
            assert!(a.cmp(&b) == core::cmp::Ordering::Equal);
        }

        fn assert_lt(a: Scalar, b: Scalar) {
            assert!(a < b);
            assert!(b > a);
        }

        assert_equality(
            Scalar::from_u64s_le(&[9999, 9999, 9999, 9999]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9999, 9999, 9999]).unwrap(),
        );
        assert_equality(
            Scalar::from_u64s_le(&[9999, 9998, 9999, 9999]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9998, 9999, 9999]).unwrap(),
        );
        assert_equality(
            Scalar::from_u64s_le(&[9999, 9999, 9999, 9997]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9999, 9999, 9997]).unwrap(),
        );
        assert_lt(
            Scalar::from_u64s_le(&[9999, 9997, 9999, 9998]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9997, 9999, 9999]).unwrap(),
        );
        assert_lt(
            Scalar::from_u64s_le(&[9999, 9997, 9998, 9999]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9997, 9999, 9999]).unwrap(),
        );
        assert_lt(
            Scalar::from_u64s_le(&[9, 9999, 9999, 9997]).unwrap(),
            Scalar::from_u64s_le(&[9999, 9999, 9999, 9997]).unwrap(),
        );
    }

    #[test]
    fn test_scalar_from_u64() {
        let a = Scalar::from(100);
        let mut expected_bytes = [0u8; 32];
        expected_bytes[0] = 100;
        assert_eq!(a.to_bytes_le(), expected_bytes);
    }

    #[test]
    fn test_scalar_is_odd() {
        assert!(bool::from(Scalar::from(0).is_even()));
        assert!(bool::from(Scalar::from(1).is_odd()));
        assert!(bool::from(Scalar::from(324834872).is_even()));
        assert!(bool::from(Scalar::from(324834873).is_odd()));
    }

    #[test]
    fn test_scalar_is_zero() {
        assert!(bool::from(Scalar::from(0).is_zero()));
        assert!(!bool::from(Scalar::from(1).is_zero()));
        assert!(!bool::from(
            Scalar::from_u64s_le(&[0, 0, 1, 0]).unwrap().is_zero()
        ));
    }

    #[test]
    fn test_scalar_num_bits() {
        assert_eq!(Scalar::NUM_BITS, 255);
        assert_eq!(Scalar::CAPACITY, 254);

        let mut a = Scalar::from(0);
        assert_eq!(0, a.num_bits());
        a = Scalar::from(1);
        assert_eq!(1, a.num_bits());
        for i in 2..Scalar::NUM_BITS {
            a = a.shl(1);
            assert_eq!(i, a.num_bits());
        }
    }

    #[test]
    fn test_scalar_legendre() {
        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        assert_eq!(Scalar::one().sqrt().unwrap(), Scalar::one());

        let e = Scalar::from_u64s_le(&[
            0x0dbc5349cd5664da,
            0x8ac5b6296e3ae29d,
            0x127cb819feceaa3b,
            0x3a6b21fb03867191,
        ])
        .unwrap();
        assert!(bool::from(e.is_quad_res()));

        let e = Scalar::from_u64s_le(&[
            0x96341aefd047c045,
            0x9b5f4254500a4d65,
            0x1ee08223b68ac240,
            0x31d9cd545c0ec7c6,
        ])
        .unwrap();
        assert!(!bool::from(e.is_quad_res()));
    }

    #[test]
    fn test_scalar_add_assign() {
        {
            // Random number
            let mut tmp = Scalar(blst::blst_fr {
                l: [
                    0x437ce7616d580765,
                    0xd42d1ccb29d1235b,
                    0xed8f753821bd1423,
                    0x4eede1c9c89528ca,
                ],
            });
            // assert!(tmp.is_valid());
            // Test that adding zero has no effect.
            tmp.add_assign(&Scalar(blst::blst_fr { l: [0, 0, 0, 0] }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ]
                })
            );
            // Add one and test for the result.
            tmp.add_assign(&Scalar(blst::blst_fr { l: [1, 0, 0, 0] }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0x437ce7616d580766,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ]
                })
            );
            // Add another random number that exercises the reduction.
            tmp.add_assign(&Scalar(blst::blst_fr {
                l: [
                    0x946f435944f7dc79,
                    0xb55e7ee6533a9b9b,
                    0x1e43b84c2f6194ca,
                    0x58717ab525463496,
                ],
            }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0xd7ec2abbb24fe3de,
                        0x35cdf7ae7d0d62f7,
                        0xd899557c477cd0e9,
                        0x3371b52bc43de018
                    ]
                })
            );
            // Add one to (r - 1) and test for the result.
            tmp = Scalar(blst::blst_fr {
                l: [
                    0xffffffff00000000,
                    0x53bda402fffe5bfe,
                    0x3339d80809a1d805,
                    0x73eda753299d7d48,
                ],
            });
            tmp.add_assign(&Scalar(blst::blst_fr { l: [1, 0, 0, 0] }));
            assert!(bool::from(tmp.is_zero()));
            // Add a random number to another one such that the result is r - 1
            tmp = Scalar(blst::blst_fr {
                l: [
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ],
            });
            tmp.add_assign(&Scalar(blst::blst_fr {
                l: [
                    0x521a525223349e70,
                    0xa99bb5f3d8231f31,
                    0xde8e397bebe477e,
                    0x1ad08e5041d7c321,
                ],
            }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0xffffffff00000000,
                        0x53bda402fffe5bfe,
                        0x3339d80809a1d805,
                        0x73eda753299d7d48
                    ]
                })
            );
            // Add one to the result and test for it.
            tmp.add_assign(&Scalar(blst::blst_fr { l: [1, 0, 0, 0] }));
            assert!(bool::from(tmp.is_zero()));
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
            let mut tmp = Scalar(blst::blst_fr {
                l: [
                    0x6a68c64b6f735a2b,
                    0xd5f4d143fe0a1972,
                    0x37c17f3829267c62,
                    0xa2f37391f30915c,
                ],
            });
            tmp.sub_assign(&Scalar(blst::blst_fr {
                l: [
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ],
            }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0xbc83189d92a7f89c,
                        0x7f908737d62d38a3,
                        0x45aa62cfe7e4c3e1,
                        0x24ffc5896108547d
                    ]
                })
            );

            // Test the opposite subtraction which doesn't test reduction.
            tmp = Scalar(blst::blst_fr {
                l: [
                    0xade5adacdccb6190,
                    0xaa21ee0f27db3ccd,
                    0x2550f4704ae39086,
                    0x591d1902e7c5ba27,
                ],
            });
            tmp.sub_assign(&Scalar(blst::blst_fr {
                l: [
                    0x6a68c64b6f735a2b,
                    0xd5f4d143fe0a1972,
                    0x37c17f3829267c62,
                    0xa2f37391f30915c,
                ],
            }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ]
                })
            );

            // Test for sensible results with zero
            tmp = Scalar(blst::blst_fr { l: [0, 0, 0, 0] });
            tmp.sub_assign(&Scalar(blst::blst_fr { l: [0, 0, 0, 0] }));
            assert!(bool::from(tmp.is_zero()));

            tmp = Scalar(blst::blst_fr {
                l: [
                    0x437ce7616d580765,
                    0xd42d1ccb29d1235b,
                    0xed8f753821bd1423,
                    0x4eede1c9c89528ca,
                ],
            });
            tmp.sub_assign(&Scalar(blst::blst_fr { l: [0, 0, 0, 0] }));
            assert_eq!(
                tmp,
                Scalar(blst::blst_fr {
                    l: [
                        0x437ce7616d580765,
                        0xd42d1ccb29d1235b,
                        0xed8f753821bd1423,
                        0x4eede1c9c89528ca
                    ]
                })
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
            assert!(bool::from(tmp1.is_zero()));
        }
    }

    #[test]
    fn test_scalar_mul_assign() {
        let mut tmp = Scalar(blst::blst_fr {
            l: [
                0x6b7e9b8faeefc81a,
                0xe30a8463f348ba42,
                0xeff3cb67a8279c9c,
                0x3d303651bd7c774d,
            ],
        });
        tmp.mul_assign(&Scalar(blst::blst_fr {
            l: [
                0x13ae28e3bc35ebeb,
                0xa10f4488075cae2c,
                0x8160e95a853c3b5d,
                0x5ae3f03b561a841d,
            ],
        }));
        assert!(
            tmp == Scalar(blst::blst_fr {
                l: [
                    0x23717213ce710f71,
                    0xdbee1fe53a16e1af,
                    0xf565d3e1c2a48000,
                    0x4426507ee75df9d7
                ]
            })
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
        let a = Scalar(blst::blst_fr {
            l: [
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x73eda753299d7d47,
            ],
        });
        // assert!(a.is_valid());
        assert_eq!(
            a.square(),
            Scalar::from_u64s_le(&[
                0xc0d698e7bde077b8,
                0xb79a310579e76ec2,
                0xac1da8d0a9af4e5f,
                0x13f629c49bf23e97
            ])
            .unwrap()
        );

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000000 {
            // Ensure that (a * a) = a^2
            let a = Scalar::random(&mut rng);

            let tmp = a.square();

            let mut tmp2 = a;
            tmp2.mul_assign(&a);

            assert_eq!(tmp, tmp2);
        }
    }

    #[test]
    fn test_scalar_inverse() {
        assert_eq!(Scalar::zero().invert().is_none().unwrap_u8(), 1);

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let one = Scalar::one();

        for i in 0..1000 {
            // Ensure that a * a^-1 = 1
            let mut a = Scalar::random(&mut rng);
            let ainv = a.invert().unwrap();
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
            let a = Scalar::random(&mut rng);
            let mut b = a;
            b.add_assign(&a);
            assert_eq!(a.double(), b);
        }
    }

    #[test]
    fn test_scalar_negate() {
        {
            let a = Scalar::zero();
            assert!(bool::from((-a).is_zero()));
        }

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..1000 {
            // Ensure (a - (-a)) = 0.
            let mut a = Scalar::random(&mut rng);
            a.add_assign(-a);
            assert!(bool::from(a.is_zero()));
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
            let target = a.pow_vartime(&[i]);
            let mut c = Scalar::one();
            for _ in 0..i {
                c.mul_assign(&a);
            }
            assert_eq!(c, target);
        }

        for _ in 0..1000 {
            // Exponentiating by the modulus should have no effect in a prime field.
            let a = Scalar::random(&mut rng);

            assert_eq!(a, a.pow_vartime(&MODULUS));
        }
    }

    #[test]
    fn test_scalar_sqrt() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        assert_eq!(Scalar::one().sqrt().unwrap(), Scalar::one());

        for _ in 0..1000 {
            // Ensure sqrt(a^2) = a or -a
            let a = Scalar::random(&mut rng);
            let a_new = a.square().sqrt().unwrap();
            assert!(a_new == a || a_new == -a);
        }

        for _ in 0..1000 {
            // Ensure sqrt(a)^2 = a for random a
            let a = Scalar::random(&mut rng);
            let sqrt = a.sqrt();
            if sqrt.is_some().into() {
                assert_eq!(sqrt.unwrap().square(), a);
            }
        }
    }

    #[test]
    fn test_scalar_from_into_repr() {
        // r + 1 should not be in the field
        assert!(bool::from(
            Scalar::from_u64s_le(&[
                0xffffffff00000002,
                0x53bda402fffe5bfe,
                0x3339d80809a1d805,
                0x73eda753299d7d48
            ])
            .is_none()
        ));

        // Modulus should not be in the field
        assert!(bool::from(Scalar::from_repr(Scalar::char()).is_none()));
        assert!(Scalar::from_repr_vartime(Scalar::char()).is_none());

        // Multiply some arbitrary representations to see if the result is as expected.
        let mut a = Scalar::from_u64s_le(&[
            0x25ebe3a3ad3c0c6a,
            0x6990e39d092e817c,
            0x941f900d42f5658e,
            0x44f8a103b38a71e0,
        ])
        .unwrap();
        let b = Scalar::from_u64s_le(&[
            0x264e9454885e2475,
            0x46f7746bb0308370,
            0x4683ef5347411f9,
            0x58838d7f208d4492,
        ])
        .unwrap();
        let c = Scalar::from_u64s_le(&[
            0x48a09ab93cfc740d,
            0x3a6600fbfc7a671,
            0x838567017501d767,
            0x7161d6da77745512,
        ])
        .unwrap();
        a.mul_assign(&b);
        assert_eq!(a, c);

        // Zero should be in the field.
        assert!(bool::from(Scalar::from_repr([0u8; 32]).unwrap().is_zero()));
        assert!(bool::from(
            Scalar::from_repr_vartime([0u8; 32]).unwrap().is_zero()
        ));

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for i in 0..1000 {
            // Try to turn Scalar elements into representations and back again, and compare.
            let a = Scalar::random(&mut rng);
            let a_again = Scalar::from_repr(a.to_repr()).unwrap();
            assert_eq!(a, a_again, "{}", i);
            let a_yet_again = Scalar::from_repr_vartime(a.to_repr()).unwrap();
            assert_eq!(a, a_yet_again);
        }
    }

    #[test]
    fn test_scalar_display() {
        assert_eq!(
            format!(
                "{}",
                Scalar::from_u64s_le(&[
                    0xc3cae746a3b5ecc7,
                    0x185ec8eb3f5b5aee,
                    0x684499ffe4b9dd99,
                    0x7c9bba7afb68faa
                ])
                .unwrap()
            ),
            "Scalar(0x07c9bba7afb68faa684499ffe4b9dd99185ec8eb3f5b5aeec3cae746a3b5ecc7)"
                .to_string()
        );
        assert_eq!(
            format!(
                "{}",
                Scalar::from_u64s_le(&[
                    0x44c71298ff198106,
                    0xb0ad10817df79b6a,
                    0xd034a80a2b74132b,
                    0x41cf9a1336f50719
                ])
                .unwrap()
            ),
            "Scalar(0x41cf9a1336f50719d034a80a2b74132bb0ad10817df79b6a44c71298ff198106)"
                .to_string()
        );
    }

    #[test]
    fn test_scalar_root_of_unity() {
        assert_eq!(Scalar::S, 32);
        assert_eq!(Scalar::multiplicative_generator(), Scalar::from(7));
        assert_eq!(
            Scalar::multiplicative_generator().pow_vartime([
                0xfffe5bfeffffffff,
                0x9a1d80553bda402,
                0x299d7d483339d808,
                0x73eda753
            ]),
            Scalar::root_of_unity()
        );
        assert_eq!(
            Scalar::root_of_unity().pow_vartime([1 << Scalar::S]),
            Scalar::one()
        );
        assert!(!bool::from(
            Scalar::multiplicative_generator().is_quad_res()
        ));
    }

    #[test]
    fn scalar_field_tests() {
        crate::tests::field::random_field_tests::<Scalar>();
        crate::tests::field::random_sqrt_tests::<Scalar>();
        crate::tests::field::from_str_tests::<Scalar>();
    }

    #[test]
    fn test_scalar_repr_conversion() {
        let a = Scalar::from(1);
        let mut expected_bytes = [0u8; 32];
        expected_bytes[0] = 1;
        assert_eq!(a, Scalar::from_repr(a.to_repr()).unwrap());
        assert_eq!(a.to_repr(), expected_bytes);
        assert_eq!(a, Scalar::from_repr(expected_bytes).unwrap());

        let a = Scalar::from(12);
        let mut expected_bytes = [0u8; 32];
        expected_bytes[0] = 12;
        assert_eq!(a, Scalar::from_repr(a.to_repr()).unwrap());
        assert_eq!(a.to_repr(), expected_bytes);
        assert_eq!(a, Scalar::from_repr(expected_bytes).unwrap());
    }

    #[test]
    fn test_scalar_repr_vartime_conversion() {
        let a = Scalar::from(1);
        let mut expected_bytes = [0u8; 32];
        expected_bytes[0] = 1;
        assert_eq!(a, Scalar::from_repr_vartime(a.to_repr()).unwrap());
        assert_eq!(a.to_repr(), expected_bytes);
        assert_eq!(a, Scalar::from_repr_vartime(expected_bytes).unwrap());

        let a = Scalar::from(12);
        let mut expected_bytes = [0u8; 32];
        expected_bytes[0] = 12;
        assert_eq!(a, Scalar::from_repr_vartime(a.to_repr()).unwrap());
        assert_eq!(a.to_repr(), expected_bytes);
        assert_eq!(a, Scalar::from_repr_vartime(expected_bytes).unwrap());
    }

    #[test]
    fn test_scalar_to_le_bits() {
        let mut bits = Scalar::one().to_le_bits().into_iter();
        assert!(bits.next().unwrap());
        for bit in bits {
            assert!(!bit);
        }

        let mut bits = Scalar::from(u64::MAX).to_le_bits().into_iter();
        for _ in 0..64 {
            assert!(bits.next().unwrap());
        }
        for _ in 64..Scalar::NUM_BITS {
            assert!(!bits.next().unwrap());
        }
        // Check that the final bit in the backing representation, i.e. the 256-th bit, is false.
        // This bit should always be `false` because it exceeds the field size modulus.
        assert!(!bits.next().unwrap());
        // Check that the bitvec's size does not exceed the size of the backing representation
        // `[u8; 32]`, i.e. 256-bits.
        assert!(bits.next().is_none());

        let mut neg1_bits = (-Scalar::one()).to_le_bits().into_iter();
        let mut modulus_bits = Scalar::char_le_bits().into_iter();
        assert_ne!(neg1_bits.next().unwrap(), modulus_bits.next().unwrap());
        for (b1, b2) in neg1_bits.zip(modulus_bits) {
            assert_eq!(b1, b2);
        }
    }

    #[test]
    fn m1_inv_bug() {
        // This fails on aarch64-darwin.
        let bad = Scalar::zero() - Scalar::from(7);

        let inv = bad.invert().unwrap();
        let check = inv * bad;
        assert_eq!(Scalar::one(), check);
    }
    #[test]
    fn m1_inv_bug_more() {
        let mut bad = Vec::new();
        for i in 1..1000000 {
            // Ensure that a * a^-1 = 1
            let a = Scalar::zero() - Scalar::from(i);
            let ainv = a.invert().unwrap();
            let check = a * ainv;
            let one = Scalar::one();

            if check != one {
                bad.push((i, a));
            }
        }
        assert_eq!(0, bad.len());
    }

    fn scalar_from_u64s(parts: [u64; 4]) -> Scalar {
        let mut le_bytes = [0u8; 32];
        le_bytes[0..8].copy_from_slice(&parts[0].to_le_bytes());
        le_bytes[8..16].copy_from_slice(&parts[1].to_le_bytes());
        le_bytes[16..24].copy_from_slice(&parts[2].to_le_bytes());
        le_bytes[24..32].copy_from_slice(&parts[3].to_le_bytes());
        let mut repr = <Scalar as PrimeField>::Repr::default();
        repr.as_mut().copy_from_slice(&le_bytes[..]);
        Scalar::from_repr_vartime(repr).expect("u64s exceed BLS12-381 scalar field modulus")
    }

    #[test]
    fn m1_inv_bug_special() {
        let maybe_bad = [scalar_from_u64s([
            0xb3fb72ea181b4e82,
            0x9435fcaf3a85c901,
            0x9eaf4fa6b9635037,
            0x2164d020b3bd14cc,
        ])];

        let mut yep_bad = Vec::new();

        for a in maybe_bad.iter() {
            let ainv = a.invert().unwrap();
            let check = a * ainv;
            let one = Scalar::one();

            if check != one {
                yep_bad.push(a);
            }
        }
        assert_eq!(0, yep_bad.len());
    }
}
