//! An implementation of the $\mathbb{G}_2$ group of BLS12-381.

use core::{
    borrow::Borrow,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use rand_core::RngCore;

use crate::{Fp2, Scalar};

/// This is an element of $\mathbb{G}_2$ represented in the affine coordinate space.
/// It is ideal to keep elements in this representation to reduce memory usage and
/// improve performance through the use of mixed curve model arithmetic.
#[derive(Copy, Clone, Debug)]
pub struct G2Affine(pub(crate) blst_p2_affine);

impl Default for G2Affine {
    fn default() -> G2Affine {
        G2Affine::zero()
    }
}

impl From<&G2Projective> for G2Affine {
    fn from(p: &G2Projective) -> G2Affine {
        let mut out = blst_p2_affine::default();

        unsafe { blst_p2_to_affine(&mut out, &p.0) };

        G2Affine(out)
    }
}

impl From<G2Projective> for G2Affine {
    fn from(p: G2Projective) -> G2Affine {
        G2Affine::from(&p)
    }
}

impl Eq for G2Affine {}
impl PartialEq for G2Affine {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        unsafe { blst_p2_affine_is_equal(&self.0, &other.0) }
    }
}

impl Neg for &G2Affine {
    type Output = G2Affine;

    #[inline]
    fn neg(self) -> G2Affine {
        // Missing for affine in blst
        todo!()
    }
}

impl Neg for G2Affine {
    type Output = G2Affine;

    #[inline]
    fn neg(self) -> G2Affine {
        -&self
    }
}

impl<'a, 'b> Add<&'b G2Projective> for &'a G2Affine {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Projective) -> G2Projective {
        rhs.add_mixed(self)
    }
}

impl<'a, 'b> Add<&'b G2Affine> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Affine) -> G2Projective {
        self.add_mixed(rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Projective> for &'a G2Affine {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Projective) -> G2Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Affine> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Affine) -> G2Projective {
        self + (-rhs)
    }
}

impl<T> Sum<T> for G2Projective
where
    T: Borrow<G2Projective>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::zero(), |acc, item| acc + item.borrow())
    }
}

impl_binops_additive!(G2Projective, G2Affine);
impl_binops_additive_specify_output!(G2Affine, G2Projective, G2Projective);

impl G2Affine {
    /// Returns the additive identity.
    pub fn zero() -> Self {
        G2Affine(blst_p2_affine::default())
    }

    /// Returns a fixed generator of unknown exponent.
    pub fn one() -> Self {
        G2Affine(unsafe { BLS12_381_G2 })
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
            blst_p2_affine_compress(out_v.as_mut_ptr(), &self.0);
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
            blst_p2_affine_serialize(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 96];
        out.copy_from_slice(&out_v);

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_uncompressed_unchecked(bytes).and_then(|el| {
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
        let mut raw = blst_p2_affine::default();

        if unsafe { blst_p2_deserialize(&mut raw, in_v.as_mut_ptr()) != BLST_ERROR::BLST_SUCCESS } {
            return None;
        }

        Some(G2Affine(raw))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> Option<Self> {
        G2Affine::from_compressed_unchecked(bytes).and_then(|el| {
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
        let mut raw = blst_p2_affine::default();

        if unsafe { blst_p2_uncompress(&mut raw, in_v.as_mut_ptr()) != BLST_ERROR::BLST_SUCCESS } {
            return None;
        }

        Some(G2Affine(raw))
    }

    /// Returns true if this point is free of an $h$-torsion component, and so it
    /// exists within the $q$-order subgroup $\mathbb{G}_2$. This should always return true
    /// unless an "unchecked" API was used.
    pub fn is_torsion_free(&self) -> bool {
        unsafe { blst_p2_affine_in_g2(&self.0) }
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> bool {
        let on_curve = unsafe { blst_p2_affine_on_curve(&self.0) };
        // FIXME: is_zero check should happen in blst
        on_curve || self.is_zero()
    }

    /// Attempts to construct an affine point given an x-coordinate. The
    /// point is not guaranteed to be in the prime order subgroup.
    ///
    /// If and only if `greatest` is set will the lexicographically
    /// largest y-coordinate be selected.
    fn get_point_from_x(x: Fp2, greatest: bool) -> Option<Self> {
        todo!()
    }

    fn scale_by_cofactor(&self) -> G2Projective {
        todo!()
    }

    pub fn from_raw_unchecked(x: Fp2, y: Fp2, infinity: bool) -> Self {
        let mut raw = blst_p2_affine::default();
        raw.x = x.0;
        raw.y = y.0;
        // FIXME: what about infinity?

        G2Affine(raw)
    }

    /// Returns the x coordinate.
    pub fn x(&self) -> Fp2 {
        Fp2(self.0.x)
    }

    /// Returns the y coordinate.
    pub fn y(&self) -> Fp2 {
        Fp2(self.0.y)
    }

    pub const fn uncompressed_size() -> usize {
        96
    }

    pub const fn compressed_size() -> usize {
        48
    }

    /// Prepare this element for use in a pairing.
    pub fn prepare(&self) -> G2Prepared {
        let mut lines = vec![blst_fp6::default(); 68];
        unsafe { blst_precompute_lines(lines.as_mut_ptr(), &self.0) }
        G2Prepared(lines)
    }
}

/// This is an element of $\mathbb{G}_2$ represented in the projective coordinate space.
#[derive(Copy, Clone, Debug)]
pub struct G2Projective(pub(crate) blst_p2);

impl From<&G2Affine> for G2Projective {
    fn from(p: &G2Affine) -> G2Projective {
        let mut out = blst_p2::default();

        unsafe { blst_p2_from_affine(&mut out, &p.0) };

        G2Projective(out)
    }
}

impl From<G2Affine> for G2Projective {
    fn from(p: G2Affine) -> G2Projective {
        G2Projective::from(&p)
    }
}

impl Eq for G2Projective {}
impl PartialEq for G2Projective {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        // TODO: more efficiente method
        G2Affine::from(self) == G2Affine::from(other)
    }
}

impl<'a> Neg for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn neg(self) -> G2Projective {
        let mut out = *self;
        const FLAG: usize = 0x1;

        unsafe { blst_p2_cneg(&mut out.0, FLAG) }

        out
    }
}

impl Neg for G2Projective {
    type Output = G2Projective;

    #[inline]
    fn neg(self) -> G2Projective {
        -&self
    }
}

impl<'a, 'b> Add<&'b G2Projective> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn add(self, rhs: &'b G2Projective) -> G2Projective {
        self.add(rhs)
    }
}

impl<'a, 'b> Sub<&'b G2Projective> for &'a G2Projective {
    type Output = G2Projective;

    #[inline]
    fn sub(self, rhs: &'b G2Projective) -> G2Projective {
        self + (-rhs)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G2Projective {
    type Output = G2Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        self.multiply(&other)
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G2Affine {
    type Output = G2Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        G2Projective::from(self).multiply(&other)
    }
}

impl_binops_additive!(G2Projective, G2Projective);
impl_binops_multiplicative!(G2Projective, Scalar);
impl_binops_multiplicative_mixed!(G2Affine, Scalar, G2Projective);

impl G2Projective {
    /// Returns the additive identity.
    pub fn zero() -> Self {
        G2Projective(blst_p2::default())
    }

    /// Returns a fixed generator of unknown exponent.
    pub fn one() -> Self {
        G2Affine::one().into()
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
            blst_p2_compress(out_v.as_mut_ptr(), &self.0);
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
            blst_p2_serialize(out_v.as_mut_ptr(), &self.0);
        }

        let mut out = [0u8; 96];
        out.copy_from_slice(&out_v);

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_uncompressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_uncompressed_unchecked(bytes).map(Into::into)
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 48]) -> Option<Self> {
        G2Affine::from_compressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 48]) -> Option<Self> {
        G2Affine::from_compressed_unchecked(bytes).map(Into::into)
    }

    /// Computes the doubling of this point.
    pub fn double(&self) -> G2Projective {
        let mut out = blst_p2::default();

        unsafe { blst_p2_add_or_double(&mut out, &self.0, &self.0) };

        G2Projective(out)
    }

    /// Adds this point to another point.
    pub fn add(&self, rhs: &G2Projective) -> G2Projective {
        let mut out = blst_p2::default();

        unsafe { blst_p2_add_or_double(&mut out, &self.0, &rhs.0) };

        G2Projective(out)
    }

    /// Adds this point to another point in the affine model.
    pub fn add_mixed(&self, rhs: &G2Affine) -> G2Projective {
        let mut out = blst_p2::default();

        unsafe { blst_p2_add_or_double_affine(&mut out, &self.0, &rhs.0) };

        G2Projective(out)
    }

    /// Returns true if this point is on the curve. This should always return
    /// true unless an "unchecked" API was used.
    pub fn is_on_curve(&self) -> bool {
        unsafe { blst_p2_on_curve(&self.0) }
    }

    fn multiply(&self, by: &Scalar) -> G2Projective {
        let mut out = blst_p2::default();

        // Sclar is 255 bits wide.
        const NBITS: usize = 255;

        // Safe, because all bslt_fr are valid blst_scalar.
        let scalar: blst_scalar = by.into();

        unsafe { blst_p2_mult(&mut out, &self.0, &scalar, NBITS) };

        G2Projective(out)
    }

    pub fn random<R: RngCore>(rng: &mut R) -> Self {
        loop {
            let x = Fp2::random(rng);
            let greatest = rng.next_u32() % 2 != 0;

            if let Some(p) = G2Affine::get_point_from_x(x, greatest) {
                let p = p.scale_by_cofactor();

                if !p.is_zero() {
                    return p;
                }
            }
        }
    }

    pub fn from_raw_unchecked(x: Fp2, y: Fp2, z: Fp2) -> Self {
        let mut raw = blst_p2::default();
        raw.x = x.0;
        raw.y = y.0;
        raw.z = z.0;

        G2Projective(raw)
    }

    /// Returns the x coordinate.
    pub fn x(&self) -> Fp2 {
        Fp2(self.0.x)
    }

    /// Returns the y coordinate.
    pub fn y(&self) -> Fp2 {
        Fp2(self.0.y)
    }

    /// Returns the z coordinate.
    pub fn z(&self) -> Fp2 {
        Fp2(self.0.z)
    }
}

#[derive(Clone, Debug)]
pub struct G2Prepared(Vec<blst_fp6>);
