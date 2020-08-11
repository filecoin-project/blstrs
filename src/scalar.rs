//! An implementation of the BLS12-381 scalar field $\mathbb{F}_q$
//! where `q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001`

use core::{
    convert::TryInto,
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;

/// Represents an element of the scalar field $\mathbb{F}_q$ of the BLS12-381 elliptic
/// curve construction.
///
/// The inner representation is stored in Montgomery form.
#[derive(Default, Clone, Copy)]
pub struct Scalar(pub(crate) blst_fr);

/// Constant representing the modulus
/// q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
pub const MODULUS: Scalar = Scalar(blst_fr {
    l: [
        0xffffffff00000001,
        0x53bda402fffe5bfe,
        0x3339d80809a1d805,
        0x73eda753299d7d48,
    ],
});

impl fmt::Debug for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes_le();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl PartialEq for Scalar {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.0.l == other.0.l
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
        if !unsafe { blst_scalar_fr_check(&self as _) } {
            return Err(NotInFieldError);
        }

        // Safe because valid fr check was just made above.
        let fr: blst_fr = unsafe { std::mem::transmute(self) };

        Ok(Scalar(fr))
    }
}

impl From<u32> for Scalar {
    fn from(val: u32) -> Scalar {
        let mut raw = blst_scalar::default();

        unsafe { blst_scalar_from_uint32(&mut raw as *mut _, val as *const _) };

        raw.try_into().expect("u32 is always inside the field")
    }
}

impl From<u64> for Scalar {
    fn from(val: u64) -> Scalar {
        let mut raw = blst_scalar::default();

        unsafe { blst_scalar_from_uint64(&mut raw as *mut _, val as *const _) };

        raw.try_into().expect("u64 is always inside the field")
    }
}

impl Eq for Scalar {}

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0xfffffffeffffffff;

/// R = 2^256 mod q
const R: Scalar = Scalar(blst_fr {
    l: [
        0x00000001fffffffe,
        0x5884b7fa00034802,
        0x998c4fefecbc4ff5,
        0x1824b159acc5056f,
    ],
});

/// R^2 = 2^512 mod q
const R2: Scalar = Scalar(blst_fr {
    l: [
        0xc999e990f3f29c6d,
        0x2b6cedcb87925c23,
        0x05d314967254398f,
        0x0748d9d99f59ff11,
    ],
});

/// R^3 = 2^768 mod q
const R3: Scalar = Scalar(blst_fr {
    l: [
        0xc62c1807439b73af,
        0x1b3e0d188cf06990,
        0x73d13c71c7b5f418,
        0x6e2a5bb9c8db33e9,
    ],
});

const S: u32 = 32;

/// GENERATOR^t where t * 2^s + 1 = q
/// with t odd. In other words, this
/// is a 2^s root of unity.
///
/// `GENERATOR = 7 mod q` is a generator
/// of the q - 1 order multiplicative
/// subgroup.
const ROOT_OF_UNITY: Scalar = Scalar(blst_fr {
    l: [
        0xb9b58d8c5f0e466a,
        0x5b1b4c801819d7ec,
        0x0af53ae352a31e64,
        0x5bf3adda19e9b27b,
    ],
});

impl<'a> Neg for &'a Scalar {
    type Output = Scalar;

    #[inline]
    fn neg(self) -> Scalar {
        self.neg()
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

impl Scalar {
    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Scalar {
        Scalar(blst_fr { l: [0, 0, 0, 0] })
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Scalar {
        R
    }

    /// Doubles this field element.
    #[inline]
    pub fn double(&self) -> Scalar {
        todo!()
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Scalar`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 32]) -> Option<Scalar> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_scalar::default();

        unsafe {
            blst_scalar_from_lendian(&mut raw as _, in_v.as_mut_ptr());
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
            blst_scalar_from_bendian(&mut raw as _, in_v.as_mut_ptr());
        }

        raw.try_into().ok()
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Scalar` representation.
    pub fn from_raw(val: [u64; 4]) -> Self {
        let mut original = blst_fr::default();
        original.l.copy_from_slice(&val);

        let mut raw = blst_fr::default();
        unsafe { blst_fr_from(&mut raw as _, &original as _) }

        Scalar(raw)
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 32] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = Vec::with_capacity(32);
        // Safe because any valid blst_fr is also a valid blst_scalar.
        let scalar: blst_scalar = unsafe { std::mem::transmute(self.0) };

        unsafe {
            blst_lendian_from_scalar(out_v.as_mut_ptr(), &scalar);
        }
        let mut out = [0u8; 32];
        out.copy_from_slice(&out_v[..32]);

        out
    }

    /// Converts an element of `Scalar` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes_be(&self) -> [u8; 32] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = Vec::with_capacity(32);
        // Safe because any valid blst_fr is also a valid blst_scalar.
        let scalar: blst_scalar = unsafe { std::mem::transmute(self.0) };
        unsafe {
            blst_bendian_from_scalar(out_v.as_mut_ptr(), &scalar);
        }
        let mut out = [0u8; 32];
        out.copy_from_slice(&out_v[..32]);

        out
    }

    /// Squares this element.
    #[inline]
    pub fn square(&self) -> Scalar {
        let mut raw = blst_fr::default();
        unsafe { blst_fr_sqr(&mut raw as _, &self.0 as _) }

        Scalar(raw)
    }

    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> Option<Self> {
        todo!()
    }

    /// Computes the multiplicative inverse of this element, failing if the element is zero.
    pub fn invert(&self) -> Option<Self> {
        todo!()
    }

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_mul(&mut out as _, &self.0 as _, &rhs.0 as _) };

        Scalar(out)
    }

    /// Exponentiates `self` by `by`, where `by` is a little-endian order integer exponent.
    pub fn pow(&self, _by: &[u64; 4]) -> Self {
        todo!()
    }

    /// Subtracts `rhs` from `self`, returning the result.
    #[inline]
    pub fn sub(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_sub(&mut out as _, &self.0 as _, &rhs.0 as _) };

        Scalar(out)
    }

    /// Adds `rhs` to `self`, returning the result.
    #[inline]
    pub fn add(&self, rhs: &Self) -> Self {
        let mut out = blst_fr::default();

        unsafe { blst_fr_add(&mut out as _, &self.0 as _, &rhs.0 as _) };

        Scalar(out)
    }

    /// Negates `self`.
    #[inline]
    pub fn neg(&self) -> Self {
        let mut out = blst_fr::default();

        let flag = 0x1;

        unsafe { blst_fr_cneg(&mut out as _, &self.0 as _, flag) };

        Scalar(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const LARGEST: Scalar = Scalar(blst_fr {
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

    #[cfg(feature = "std")]
    #[test]
    fn test_debug() {
        assert_eq!(
            format!("{:?}", Scalar::zero()),
            "0x0000000000000000000000000000000000000000000000000000000000000000"
        );
        assert_eq!(
            format!("{:?}", Scalar::one()),
            "0x0000000000000000000000000000000000000000000000000000000000000001"
        );
        assert_eq!(
            format!("{:?}", R2),
            "0x1824b159acc5056f998c4fefecbc4ff55884b7fa0003480200000001fffffffe"
        );
    }

    #[test]
    fn test_equality() {
        assert_eq!(Scalar::zero(), Scalar::zero());
        assert_eq!(Scalar::one(), Scalar::one());
        assert_eq!(R2, R2);

        assert!(Scalar::zero() != Scalar::one());
        assert!(Scalar::one() != R2);
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
            Scalar(blst_fr {
                l: [
                    0xfffffffeffffffff,
                    0x53bda402fffe5bfe,
                    0x3339d80809a1d805,
                    0x73eda753299d7d48
                ]
            })
        );

        let mut tmp = LARGEST;
        tmp += &Scalar(blst_fr { l: [1, 0, 0, 0] });

        assert_eq!(tmp, Scalar::zero());
    }

    #[test]
    fn test_negation() {
        let tmp = -&LARGEST;

        assert_eq!(tmp, Scalar(blst_fr { l: [1, 0, 0, 0] }));

        let tmp = -&Scalar::zero();
        assert_eq!(tmp, Scalar::zero());
        let tmp = -&Scalar(blst_fr { l: [1, 0, 0, 0] });
        assert_eq!(tmp, LARGEST);
    }

    #[test]
    fn test_subtraction() {
        let mut tmp = LARGEST;
        tmp -= &LARGEST;

        assert_eq!(tmp, Scalar::zero());

        let mut tmp = Scalar::zero();
        tmp -= &LARGEST;

        let mut tmp2 = MODULUS;
        tmp2 -= &LARGEST;

        assert_eq!(tmp, tmp2);
    }

    #[test]
    fn test_multiplication() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp *= &cur;

            let mut tmp2 = Scalar::zero();
            for b in cur
                .to_bytes_le()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_squaring() {
        let mut cur = LARGEST;

        for _ in 0..100 {
            let mut tmp = cur;
            tmp = tmp.square();

            let mut tmp2 = Scalar::zero();
            for b in cur
                .to_bytes_le()
                .iter()
                .rev()
                .flat_map(|byte| (0..8).rev().map(move |i| ((byte >> i) & 1u8) == 1u8))
            {
                let tmp3 = tmp2;
                tmp2.add_assign(&tmp3);

                if b {
                    tmp2.add_assign(&cur);
                }
            }

            assert_eq!(tmp, tmp2);

            cur.add_assign(&LARGEST);
        }
    }

    #[test]
    fn test_inversion() {
        assert!(Scalar::zero().invert().is_none());
        assert_eq!(Scalar::one().invert().unwrap(), Scalar::one());
        assert_eq!((-&Scalar::one()).invert().unwrap(), -&Scalar::one());

        let mut tmp = R2;

        for _ in 0..100 {
            let mut tmp2 = tmp.invert().unwrap();
            tmp2.mul_assign(&tmp);

            assert_eq!(tmp2, Scalar::one());

            tmp.add_assign(&R2);
        }
    }

    #[test]
    fn test_invert_is_pow() {
        let q_minus_2 = [
            0xfffffffeffffffff,
            0x53bda402fffe5bfe,
            0x3339d80809a1d805,
            0x73eda753299d7d48,
        ];

        let mut r1 = R;
        let mut r2 = R;
        let mut r3 = R;

        for _ in 0..100 {
            r1 = r1.invert().unwrap();
            r2 = r3.pow(&q_minus_2);

            assert_eq!(r1, r2);
            // Add R so we check something different next time around
            r1.add_assign(&R);
            r2 = r1;
        }
    }

    #[test]
    fn test_sqrt() {
        {
            assert_eq!(Scalar::zero().sqrt().unwrap(), Scalar::zero());
        }

        let mut square = Scalar(blst_fr {
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

        assert_eq!(a.double(), a + a);
    }
}
