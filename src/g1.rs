//! An implementation of the $\mathbb{G}_1$ group of BLS12-381.

use core::{
    borrow::Borrow,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use rand_core::RngCore;

use crate::{Fp, Scalar};

/// This is an element of $\mathbb{G}_1$ represented in the affine coordinate space.
/// It is ideal to keep elements in this representation to reduce memory usage and
/// improve performance through the use of mixed curve model arithmetic.
#[derive(Copy, Clone, Debug)]
pub struct G1Affine(pub(crate) blst_p1_affine);

impl Default for G1Affine {
    fn default() -> G1Affine {
        G1Affine::zero()
    }
}

impl From<&G1Projective> for G1Affine {
    fn from(p: &G1Projective) -> G1Affine {
        let mut out = blst_p1_affine::default();

        unsafe { blst_p1_to_affine(&mut out, &p.0) };

        G1Affine(out)
    }
}

impl From<G1Projective> for G1Affine {
    fn from(p: G1Projective) -> G1Affine {
        G1Affine::from(&p)
    }
}

impl Eq for G1Affine {}
impl PartialEq for G1Affine {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p1_affine_is_equal(&self.0, &other.0) }
    }
}

impl Neg for &G1Affine {
    type Output = G1Affine;

    #[inline]
    fn neg(self) -> G1Affine {
        // Missing for affine in blst
        todo!()
    }
}

impl Neg for G1Affine {
    type Output = G1Affine;

    #[inline]
    fn neg(self) -> G1Affine {
        -&self
    }
}

impl<'a, 'b> Add<&'b G1Projective> for &'a G1Affine {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Projective) -> G1Projective {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b G1Affine> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Affine) -> G1Projective {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Projective> for &'a G1Affine {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Projective) -> G1Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Affine> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Affine) -> G1Projective {
        self + (-rhs)
    }
}

impl<T> Sum<T> for G1Projective
where
    T: Borrow<G1Projective>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl_binops_additive!(G1Projective, G1Affine);
impl_binops_additive_specify_output!(G1Affine, G1Projective, G1Projective);

impl G1Affine {
    /// Returns the additive identity.
    pub fn zero() -> Self {
        G1Affine(blst_p1_affine::default())
    }

    /// Returns a fixed generator of unknown exponent.
    pub fn one() -> Self {
        G1Affine(unsafe { BLS12_381_G1 })
    }

    /// Determines if this point represents the point at infinity; the additive identity.
    pub fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    /// Serializes this element into compressed form.
    pub fn to_compressed(&self) -> [u8; 48] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 48];

        unsafe {
            blst_p1_affine_compress(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 48];
        out.copy_from_slice(&out_v);

        out
    }

    /// Serializes this element into uncompressed form.
    pub fn to_uncompressed(&self) -> [u8; 96] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 96];

        unsafe {
            blst_p1_affine_serialize(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 96];
        out.copy_from_slice(&out_v);

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> Option<Self> {
        G1Affine::from_uncompressed_unchecked(bytes).and_then(|el| {
            if el.is_torsion_free() && el.is_on_curve() {
                Some(el)
            } else {
                None
            }
        })
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_p1_affine::default();

        if unsafe {
            blst_p1_deserialize(&mut raw as _, in_v.as_mut_ptr()) != BLST_ERROR::BLST_SUCCESS
        } {
            return None;
        }

        Some(G1Affine(raw))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> Option<Self> {
        G1Affine::from_compressed_unchecked(bytes).and_then(|el| {
            if el.is_torsion_free() && el.is_on_curve() {
                Some(el)
            } else {
                None
            }
        })
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> Option<Self> {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut in_v = bytes.to_vec();
        let mut raw = blst_p1_affine::default();

        if unsafe {
            blst_p1_uncompress(&mut raw as _, in_v.as_mut_ptr()) != BLST_ERROR::BLST_SUCCESS
        } {
            return None;
        }

        Some(G1Affine(raw))
    }

    /// Returns true if this point is free of an $h$-torsion component, and so it
    /// exists within the $q$-order subgroup $\mathbb{G}_1$. This should always return true
    /// unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> bool {
        unsafe { blst_p1_affine_in_g1(&self.0) }
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> bool {
        let on_curve = unsafe { blst_p1_affine_on_curve(&self.0) };
        // FIXME: is_zero check should happen in blst
        on_curve || self.is_zero()
    }

    /// Attempts to construct an affine point given an x-coordinate. The
    /// point is not guaranteed to be in the prime order subgroup.
    ///
    /// If and only if `greatest` is set will the lexicographically
    /// largest y-coordinate be selected.
    fn get_point_from_x(x: Fp, greatest: bool) -> Option<Self> {
        todo!()
    }

    fn scale_by_cofactor(&self) -> G1Projective {
        todo!()
    }

    pub fn from_raw_unchecked(x: Fp, y: Fp, infinity: bool) -> Self {
        let mut raw = blst_p1_affine::default();
        raw.x = x.0;
        raw.y = y.0;
        // FIXME: what about infinity?

        G1Affine(raw)
    }

    /// Returns the x coordinate.
    pub fn x(&self) -> Fp {
        Fp(self.0.x)
    }

    /// Returns the y coordinate.
    pub fn y(&self) -> Fp {
        Fp(self.0.y)
    }

    pub const fn uncompressed_size() -> usize {
        96
    }

    pub const fn compressed_size() -> usize {
        48
    }
}

/// This is an element of $\mathbb{G}_1$ represented in the projective coordinate space.
#[derive(Copy, Clone, Debug)]
pub struct G1Projective(pub(crate) blst_p1);

impl From<&G1Affine> for G1Projective {
    fn from(p: &G1Affine) -> G1Projective {
        let mut out = blst_p1::default();

        unsafe { blst_p1_from_affine(&mut out, &p.0) };

        G1Projective(out)
    }
}

impl From<G1Affine> for G1Projective {
    fn from(p: G1Affine) -> G1Projective {
        G1Projective::from(&p)
    }
}

impl Eq for G1Projective {}
impl PartialEq for G1Projective {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        // TODO: more efficiente method
        G1Affine::from(self) == G1Affine::from(other)
    }
}

impl<'a> Neg for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn neg(self) -> G1Projective {
        let mut out = *self;
        const FLAG: usize = 0x1;

        unsafe { blst_p1_cneg(&mut out.0, FLAG) }

        out
    }
}

impl Neg for G1Projective {
    type Output = G1Projective;

    #[inline]
    fn neg(self) -> G1Projective {
        -&self
    }
}

impl<'a, 'b> Add<&'b G1Projective> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn add(self, rhs: &'b G1Projective) -> G1Projective {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b G1Projective> for &'a G1Projective {
    type Output = G1Projective;

    #[inline]
    fn sub(self, rhs: &'b G1Projective) -> G1Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G1Projective {
    type Output = G1Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G1Affine {
    type Output = G1Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        G1Projective::from(self).multiply(&other)
    }
}

impl_binops_additive!(G1Projective, G1Projective);
impl_binops_multiplicative!(G1Projective, Scalar);
impl_binops_multiplicative_mixed!(G1Affine, Scalar, G1Projective);

impl G1Projective {
    /// Returns the additive identity.
    pub fn zero() -> Self {
        G1Projective(blst_p1::default())
    }

    /// Returns a fixed generator of unknown exponent.
    pub fn one() -> Self {
        G1Affine::one().into()
    }

    /// Determines if this point represents the point at infinity; the additive identity.
    pub fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    /// Serializes this element into compressed form.
    pub fn to_compressed(&self) -> [u8; 48] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 48];

        unsafe {
            blst_p1_compress(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 48];
        out.copy_from_slice(&out_v);

        out
    }

    /// Serializes this element into uncompressed form.
    pub fn to_uncompressed(&self) -> [u8; 96] {
        // TODO: figure out if there is a way to avoid this heap allocation
        let mut out_v = vec![0u8; 96];

        unsafe {
            blst_p1_serialize(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 96];
        out.copy_from_slice(&out_v);

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> Option<Self> {
        G1Affine::from_uncompressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        G1Affine::from_uncompressed_unchecked(bytes).map(Into::into)
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> Option<Self> {
        G1Affine::from_compressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> Option<Self> {
        G1Affine::from_compressed_unchecked(bytes).map(Into::into)
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> G1Projective {
        let mut out = blst_p1::default();

        unsafe { blst_p1_add_or_double(&mut out, &self.0, &self.0) };

        G1Projective(out)
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &G1Projective) -> G1Projective {
        let mut out = blst_p1::default();

        unsafe { blst_p1_add_or_double(&mut out, &self.0, &rhs.0) };

        G1Projective(out)
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &G1Affine) -> G1Projective {
        let mut out = blst_p1::default();

        unsafe { blst_p1_add_or_double_affine(&mut out, &self.0, &rhs.0) };

        G1Projective(out)
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> bool {
        unsafe { blst_p1_on_curve(&self.0) }
    }

    fn multiply(&self, by: &Scalar) -> G1Projective {
        let mut out = blst_p1::default();

        // Sclar is 255 bits wide.
        const NBITS: usize = 255;

        // Safe, because all bslt_fr are valid blst_scalar.
        let scalar: blst_scalar = by.into();

        unsafe { blst_p1_mult(&mut out, &self.0, &scalar, NBITS) };

        G1Projective(out)
    }

    pub fn random<R: RngCore>(rng: &mut R) -> Self {
        loop {
            let x = Fp::random(rng);
            let greatest = rng.next_u32() % 2 != 0;

            if let Some(p) = G1Affine::get_point_from_x(x, greatest) {
                let p = p.scale_by_cofactor();

                if !p.is_zero() {
                    return p;
                }
            }
        }
    }

    pub fn from_raw_unchecked(x: Fp, y: Fp, z: Fp) -> Self {
        let mut raw = blst_p1::default();
        raw.x = x.0;
        raw.y = y.0;
        raw.z = z.0;

        G1Projective(raw)
    }

    /// Returns the x coordinate.
    pub fn x(&self) -> Fp {
        Fp(self.0.x)
    }

    /// Returns the y coordinate.
    pub fn y(&self) -> Fp {
        Fp(self.0.y)
    }

    /// Returns the z coordinate.
    pub fn z(&self) -> Fp {
        Fp(self.0.z)
    }

    pub fn batch_normalization<S: std::borrow::BorrowMut<Self>>(v: &mut [S]) {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn curve_tests() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        // Negation edge case with zero.
        {
            let mut z = G1Projective::zero();
            z = z.neg();
            assert!(z.is_zero());
        }

        // Doubling edge case with zero.
        {
            let mut z = G1Projective::zero();
            z = z.double();
            assert!(z.is_zero());
        }

        // Addition edge cases with zero
        {
            let mut r = G1Projective::random(&mut rng);
            let rcopy = r;
            r += &G1Projective::zero();
            assert_eq!(r, rcopy);
            r += &G1Affine::zero();
            assert_eq!(r, rcopy);

            let mut z = G1Projective::zero();
            z += &G1Projective::zero();
            assert!(z.is_zero());
            z += &G1Affine::zero();
            assert!(z.is_zero());

            let mut z2 = z;
            z2.add_assign(&r);

            z += &G1Affine::from(r);

            assert_eq!(z, z2);
            assert_eq!(z, r);
        }

        // Transformations
        {
            let a = G1Projective::random(&mut rng);
            let b: G1Projective = G1Affine::from(a).into();
            let c = G1Projective::from(G1Affine::from(G1Projective::from(G1Affine::from(a))));

            assert_eq!(a, b);
            assert_eq!(b, c);
        }
    }

    #[test]
    fn test_is_on_curve() {
        assert!(G1Projective::zero().is_on_curve());
        assert!(G1Projective::one().is_on_curve());

        assert!(G1Affine::zero().is_on_curve());
        assert!(G1Affine::one().is_on_curve());

        let z = Fp::from_raw_unchecked([
            0xba7afa1f9a6fe250,
            0xfa0f5b595eafe731,
            0x3bdc477694c306e7,
            0x2149be4b3949fa24,
            0x64aa6e0649b2078c,
            0x12b108ac33643c3e,
        ]);

        let gen = G1Affine::one();
        let mut test =
            G1Projective::from_raw_unchecked(gen.x() * (z.square()), gen.y() * (z.square() * z), z);

        assert!(test.is_on_curve());

        test.0.x = z.0;
        assert!(!test.is_on_curve());
    }

    #[test]
    fn test_affine_point_equality() {
        let a = G1Affine::one();
        let b = G1Affine::zero();

        assert!(a == a);
        assert!(b == b);
        assert!(a != b);
        assert!(b != a);
    }

    #[test]
    fn test_projective_point_equality() {
        let a = G1Projective::one();
        let b = G1Projective::zero();

        assert!(a == a);
        assert!(b == b);
        assert!(a != b);
        assert!(b != a);

        let z = Fp::from_raw_unchecked([
            0xba7afa1f9a6fe250,
            0xfa0f5b595eafe731,
            0x3bdc477694c306e7,
            0x2149be4b3949fa24,
            0x64aa6e0649b2078c,
            0x12b108ac33643c3e,
        ]);

        let mut c =
            G1Projective::from_raw_unchecked(a.x() * (z.square()), a.y() * (z.square() * z), z);
        assert!(c.is_on_curve());

        assert!(a == c);
        assert!(b != c);
        assert!(c == a);
        assert!(c != b);

        c.0.y = (-c.y()).0;
        assert!(c.is_on_curve());

        assert!(a != c);
        assert!(b != c);
        assert!(c != a);
        assert!(c != b);

        c.0.y = (-c.y()).0;
        c.0.x = z.0;
        assert!(!c.is_on_curve());
        assert!(a != b);
        assert!(a != c);
        assert!(b != c);
    }

    #[test]
    fn test_projective_to_affine() {
        let a = G1Projective::one();
        let b = G1Projective::zero();

        assert!(G1Affine::from(a).is_on_curve());
        assert!(!G1Affine::from(a).is_zero());
        assert!(G1Affine::from(b).is_on_curve());
        assert!(G1Affine::from(b).is_zero());

        let z = Fp::from_raw_unchecked([
            0xba7afa1f9a6fe250,
            0xfa0f5b595eafe731,
            0x3bdc477694c306e7,
            0x2149be4b3949fa24,
            0x64aa6e0649b2078c,
            0x12b108ac33643c3e,
        ]);

        let c = G1Projective::from_raw_unchecked(a.x() * (z.square()), a.y() * (z.square() * z), z);

        assert_eq!(G1Affine::from(c), G1Affine::one());
    }

    #[test]
    fn test_affine_to_projective() {
        let a = G1Affine::one();
        let b = G1Affine::zero();

        assert!(G1Projective::from(a).is_on_curve());
        assert!(!G1Projective::from(a).is_zero());
        assert!(G1Projective::from(b).is_on_curve());
        assert!(G1Projective::from(b).is_zero());
    }

    #[test]
    fn test_doubling() {
        {
            let tmp = G1Projective::zero().double();
            assert!(tmp.is_zero());
            assert!(tmp.is_on_curve());
        }
        {
            let tmp = G1Projective::one().double();
            assert!(!tmp.is_zero());
            assert!(tmp.is_on_curve());

            assert_eq!(
                G1Affine::from(tmp),
                G1Affine::from_raw_unchecked(
                    Fp::from_raw_unchecked([
                        0x53e978ce58a9ba3c,
                        0x3ea0583c4f3d65f9,
                        0x4d20bb47f0012960,
                        0xa54c664ae5b2b5d9,
                        0x26b552a39d7eb21f,
                        0x8895d26e68785
                    ]),
                    Fp::from_raw_unchecked([
                        0x70110b3298293940,
                        0xda33c5393f1f6afc,
                        0xb86edfd16a5aa785,
                        0xaec6d1c9e7b1c895,
                        0x25cfc2b522d11720,
                        0x6361c83f8d09b15
                    ]),
                    false
                )
            );
        }
    }

    #[test]
    fn test_projective_addition() {
        {
            let a = G1Projective::zero();
            let b = G1Projective::zero();
            let c = a + b;
            assert!(c.is_zero());
            assert!(c.is_on_curve());
        }
        {
            let a = G1Projective::zero();
            let mut b = G1Projective::one();
            {
                let z = Fp::from_raw_unchecked([
                    0xba7afa1f9a6fe250,
                    0xfa0f5b595eafe731,
                    0x3bdc477694c306e7,
                    0x2149be4b3949fa24,
                    0x64aa6e0649b2078c,
                    0x12b108ac33643c3e,
                ]);

                b = G1Projective::from_raw_unchecked(
                    b.x() * (z.square()),
                    b.y() * (z.square() * z),
                    z,
                );
            }
            let c = a + b;
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(c == G1Projective::one());
        }
        {
            let a = G1Projective::zero();
            let mut b = G1Projective::one();
            {
                let z = Fp::from_raw_unchecked([
                    0xba7afa1f9a6fe250,
                    0xfa0f5b595eafe731,
                    0x3bdc477694c306e7,
                    0x2149be4b3949fa24,
                    0x64aa6e0649b2078c,
                    0x12b108ac33643c3e,
                ]);

                b = G1Projective::from_raw_unchecked(
                    b.x() * (z.square()),
                    b.y() * (z.square() * z),
                    z,
                );
            }
            let c = b + a;
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(c == G1Projective::one());
        }
        {
            let a = G1Projective::one().double().double(); // 4P
            let b = G1Projective::one().double(); // 2P
            let c = a + b;

            let mut d = G1Projective::one();
            for _ in 0..5 {
                d = d + G1Projective::one();
            }
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(!d.is_zero());
            assert!(d.is_on_curve());
            assert_eq!(c, d);
        }

        // Degenerate case
        {
            let beta = Fp::from_raw_unchecked([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ]);
            let beta = beta.square();
            let a = G1Projective::one().double().double();
            let b = G1Projective::from_raw_unchecked(a.x() * beta, -a.y(), a.z());
            assert!(a.is_on_curve());
            assert!(b.is_on_curve());

            let c = a + b;
            assert_eq!(
                G1Affine::from(c),
                G1Affine::from(G1Projective::from_raw_unchecked(
                    Fp::from_raw_unchecked([
                        0x29e1e987ef68f2d0,
                        0xc5f3ec531db03233,
                        0xacd6c4b6ca19730f,
                        0x18ad9e827bc2bab7,
                        0x46e3b2c5785cc7a9,
                        0x7e571d42d22ddd6
                    ]),
                    Fp::from_raw_unchecked([
                        0x94d117a7e5a539e7,
                        0x8e17ef673d4b5d22,
                        0x9d746aaf508a33ea,
                        0x8c6d883d2516c9a2,
                        0xbc3b8d5fb0447f7,
                        0x7bfa4c7210f4f44
                    ]),
                    Fp::one(),
                ))
            );
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
        }
    }

    #[test]
    fn test_mixed_addition() {
        {
            let a = G1Affine::zero();
            let b = G1Projective::zero();
            let c = a + b;
            assert!(c.is_zero());
            assert!(c.is_on_curve());
        }
        {
            let a = G1Affine::zero();
            let mut b = G1Projective::one();
            {
                let z = Fp::from_raw_unchecked([
                    0xba7afa1f9a6fe250,
                    0xfa0f5b595eafe731,
                    0x3bdc477694c306e7,
                    0x2149be4b3949fa24,
                    0x64aa6e0649b2078c,
                    0x12b108ac33643c3e,
                ]);

                b = G1Projective::from_raw_unchecked(
                    b.x() * (z.square()),
                    b.y() * (z.square() * z),
                    z,
                );
            }
            let c = a + b;
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(c == G1Projective::one());
        }
        {
            let a = G1Affine::zero();
            let mut b = G1Projective::one();
            {
                let z = Fp::from_raw_unchecked([
                    0xba7afa1f9a6fe250,
                    0xfa0f5b595eafe731,
                    0x3bdc477694c306e7,
                    0x2149be4b3949fa24,
                    0x64aa6e0649b2078c,
                    0x12b108ac33643c3e,
                ]);

                b = G1Projective::from_raw_unchecked(
                    b.x() * (z.square()),
                    b.y() * (z.square() * z),
                    z,
                );
            }
            let c = b + a;
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(c == G1Projective::one());
        }
        {
            let a = G1Projective::one().double().double(); // 4P
            let b = G1Projective::one().double(); // 2P
            let c = a + b;

            let mut d = G1Projective::one();
            for _ in 0..5 {
                d = d + G1Affine::one();
            }
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
            assert!(!d.is_zero());
            assert!(d.is_on_curve());
            assert_eq!(c, d);
        }

        // Degenerate case
        {
            let beta = Fp::from_raw_unchecked([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x3f97d6e83d050d2,
                0x18f0206554638741,
            ]);
            let beta = beta.square();
            let a = G1Projective::one().double().double();
            let b = G1Projective::from_raw_unchecked(a.x() * beta, -a.y(), a.z());
            let a = G1Affine::from(a);
            assert!(a.is_on_curve());
            assert!(b.is_on_curve());

            let c = a + b;
            assert_eq!(
                G1Affine::from(c),
                G1Affine::from(G1Projective::from_raw_unchecked(
                    Fp::from_raw_unchecked([
                        0x29e1e987ef68f2d0,
                        0xc5f3ec531db03233,
                        0xacd6c4b6ca19730f,
                        0x18ad9e827bc2bab7,
                        0x46e3b2c5785cc7a9,
                        0x7e571d42d22ddd6
                    ]),
                    Fp::from_raw_unchecked([
                        0x94d117a7e5a539e7,
                        0x8e17ef673d4b5d22,
                        0x9d746aaf508a33ea,
                        0x8c6d883d2516c9a2,
                        0xbc3b8d5fb0447f7,
                        0x7bfa4c7210f4f44
                    ]),
                    Fp::one()
                ))
            );
            assert!(!c.is_zero());
            assert!(c.is_on_curve());
        }
    }

    #[test]
    fn test_projective_negation_and_subtraction() {
        let a = G1Projective::one().double();
        assert_eq!(a + (-a), G1Projective::zero());
        assert_eq!(a + (-a), a - a);
    }

    #[test]
    fn test_affine_negation_and_subtraction() {
        let a = G1Affine::one();
        assert_eq!(G1Projective::from(a) + (-a), G1Projective::zero());
        assert_eq!(G1Projective::from(a) + (-a), G1Projective::from(a) - a);
    }

    #[test]
    fn test_projective_scalar_multiplication() {
        let g = G1Projective::one();
        let a = Scalar::from_raw([
            0x2b568297a56da71c,
            0xd8c39ecb0ef375d1,
            0x435c38da67bfbf96,
            0x8088a05026b659b2,
        ]);
        let b = Scalar::from_raw([
            0x785fdd9b26ef8b85,
            0xc997f25837695c18,
            0x4c8dbc39e7b756c1,
            0x70d9b6cc6d87df20,
        ]);
        let c = a * b;

        assert_eq!((g * a) * b, g * c);
    }

    #[test]
    fn test_affine_scalar_multiplication() {
        let g = G1Affine::one();
        let a = Scalar::from_raw([
            0x2b568297a56da71c,
            0xd8c39ecb0ef375d1,
            0x435c38da67bfbf96,
            0x8088a05026b659b2,
        ]);
        let b = Scalar::from_raw([
            0x785fdd9b26ef8b85,
            0xc997f25837695c18,
            0x4c8dbc39e7b756c1,
            0x70d9b6cc6d87df20,
        ]);
        let c = a * b;

        assert_eq!(G1Affine::from(g * a) * b, g * c);
    }
}
