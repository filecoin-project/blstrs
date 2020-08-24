//! An implementation of the $\mathbb{G}_2$ group of BLS12-381.

use core::{
    borrow::Borrow,
    fmt,
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use fff::{Field, PrimeField, PrimeFieldRepr};
use rand_core::RngCore;

use crate::{Fp12, Fp2, G1Affine, Scalar, ScalarRepr};

/// This is an element of $\mathbb{G}_2$ represented in the affine coordinate space.
/// It is ideal to keep elements in this representation to reduce memory usage and
/// improve performance through the use of mixed curve model arithmetic.
#[derive(Copy, Clone, Debug)]
pub struct G2Affine(pub(crate) blst_p2_affine);

impl fmt::Display for G2Affine {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            write!(f, "G2Affine(Infinity)")
        } else {
            write!(f, "G2Affine(x={}, y={})", self.x(), self.y())
        }
    }
}

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
        let mut res = *self;

        // Missing for affine in blst
        if !self.is_zero() {
            let mut y = res.y();
            y.negate();
            res.0.y = y.0;
        }

        res
    }
}

impl Neg for G2Affine {
    type Output = G2Affine;

    #[inline]
    fn neg(self) -> G2Affine {
        (&self).neg()
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

impl groupy::CurveAffine for G2Affine {
    type Engine = crate::Bls12;
    type Scalar = Scalar;
    type Base = Fp2;
    type Projective = G2Projective;
    type Uncompressed = G2Uncompressed;
    type Compressed = G2Compressed;

    fn zero() -> Self {
        G2Affine(blst_p2_affine::default())
    }

    fn one() -> Self {
        G2Affine(unsafe { BLS12_381_G2 })
    }

    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    fn mul<S: Into<<Self::Scalar as PrimeField>::Repr>>(&self, by: S) -> Self::Projective {
        G2Projective::from(self).multiply(&by.into())
    }

    fn negate(&mut self) {
        *self = self.neg();
    }

    fn into_projective(&self) -> Self::Projective {
        (*self).into()
    }
}

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
    pub fn to_compressed(&self) -> [u8; 96] {
        let mut out = [0u8; 96];

        unsafe {
            blst_p2_affine_compress(out.as_mut_ptr(), &self.0);
        }

        out
    }

    /// Serializes this element into uncompressed form.
    pub fn to_uncompressed(&self) -> [u8; 192] {
        let mut out = [0u8; 192];

        unsafe {
            blst_p2_affine_serialize(out.as_mut_ptr(), &self.0);
        }

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 192]) -> Option<Self> {
        G2Affine::from_uncompressed_unchecked(bytes).and_then(|el| {
            if el.is_zero() || (el.is_torsion_free() && el.is_on_curve()) {
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
    pub fn from_uncompressed_unchecked(bytes: &[u8; 192]) -> Option<Self> {
        if bytes.iter().all(|&b| b == 0) {
            return Some(Self::zero());
        }

        let mut raw = blst_p2_affine::default();

        if unsafe { blst_p2_deserialize(&mut raw, bytes.as_ptr()) != BLST_ERROR::BLST_SUCCESS } {
            return None;
        }

        Some(G2Affine(raw))
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_compressed_unchecked(bytes).and_then(|el| {
            if el.is_zero() || (el.is_torsion_free() && el.is_on_curve()) {
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
    pub fn from_compressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        if bytes.iter().all(|&b| b == 0) {
            return Some(Self::zero());
        }

        let mut raw = blst_p2_affine::default();

        if unsafe { blst_p2_uncompress(&mut raw, bytes.as_ptr()) != BLST_ERROR::BLST_SUCCESS } {
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

    pub fn from_raw_unchecked(x: Fp2, y: Fp2, _infinity: bool) -> Self {
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
        192
    }

    pub const fn compressed_size() -> usize {
        96
    }

    fn perform_pairing(&self, other: &G1Affine) -> Fp12 {
        use crate::Engine;

        crate::Bls12::pairing(*other, *self)
    }
}

/// This is an element of $\mathbb{G}_2$ represented in the projective coordinate space.
#[derive(Copy, Clone, Debug)]
pub struct G2Projective(pub(crate) blst_p2);

impl fmt::Display for G2Projective {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", G2Affine::from(self))
    }
}

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
        self.multiply(&other.into_repr())
    }
}

impl<'a, 'b> Mul<&'b Scalar> for &'a G2Affine {
    type Output = G2Projective;

    fn mul(self, other: &'b Scalar) -> Self::Output {
        G2Projective::from(self).multiply(&other.into_repr())
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
        let mut out = [0u8; 48];

        unsafe {
            blst_p2_compress(out.as_mut_ptr(), &self.0);
        }

        out
    }

    /// Serializes this element into uncompressed form.
    pub fn to_uncompressed(&self) -> [u8; 96] {
        let mut out = [0u8; 96];

        unsafe {
            blst_p2_serialize(out.as_mut_ptr(), &self.0);
        }

        out
    }

    /// Attempts to deserialize an uncompressed element.
    pub fn from_uncompressed(bytes: &[u8; 192]) -> Option<Self> {
        G2Affine::from_uncompressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is on the curve and not checking if it is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_uncompressed()` instead.
    pub fn from_uncompressed_unchecked(bytes: &[u8; 192]) -> Option<Self> {
        G2Affine::from_uncompressed_unchecked(bytes).map(Into::into)
    }

    /// Attempts to deserialize a compressed element.
    pub fn from_compressed(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_compressed(bytes).map(Into::into)
    }

    /// Attempts to deserialize an uncompressed element, not checking if the
    /// element is in the correct subgroup.
    ///
    /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
    /// API invariants may be broken.** Please consider using `from_compressed()` instead.
    pub fn from_compressed_unchecked(bytes: &[u8; 96]) -> Option<Self> {
        G2Affine::from_compressed_unchecked(bytes).map(Into::into)
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

    fn multiply(&self, by: &ScalarRepr) -> G2Projective {
        let mut out = blst_p2::default();

        // Sclar is 255 bits wide.
        const NBITS: usize = 255;

        // Safe, because all bslt_fr are valid blst_scalar.
        let scalar: blst_scalar = unsafe { std::mem::transmute(by.0) };

        unsafe { blst_p2_mult(&mut out, &self.0, &scalar, NBITS) };

        G2Projective(out)
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

impl groupy::CurveProjective for G2Projective {
    type Engine = crate::Bls12;
    type Scalar = Scalar;
    type Base = Fp2;
    type Affine = G2Affine;

    fn random<R: RngCore>(rng: &mut R) -> Self {
        let mut out = blst_p2::default();
        let mut msg = [0u8; 64];
        rng.fill_bytes(&mut msg);
        const DST: [u8; 16] = [0; 16];
        const AUG: [u8; 16] = [0; 16];

        unsafe {
            blst_encode_to_g2(
                &mut out,
                msg.as_ptr(),
                msg.len(),
                DST.as_ptr(),
                DST.len(),
                AUG.as_ptr(),
                AUG.len(),
            )
        };

        G2Projective(out)
    }

    fn zero() -> Self {
        // The point at infinity is always represented by Z = 0.
        G2Projective(blst_p2::default())
    }

    fn one() -> Self {
        G2Affine::one().into()
    }

    // The point at infinity is always represented by
    // Z = 0.
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }

    fn is_normalized(&self) -> bool {
        self.is_zero() || self.z() == Fp2::one()
    }

    fn batch_normalization<S: std::borrow::BorrowMut<Self>>(v: &mut [S]) {
        for el in v {
            let el = el.borrow_mut();
            let mut out = blst_p2_affine::default();

            unsafe { blst_p2_to_affine(&mut out, &el.0) };

            el.0.x = out.x;
            el.0.y = out.y;
            el.0.z = Fp2::one().0;
        }
    }

    fn double(&mut self) {
        let mut out = blst_p2::default();

        unsafe { blst_p2_double(&mut out, &self.0) };

        self.0 = out;
    }

    fn add_assign(&mut self, other: &Self) {
        *self += other;
    }

    fn add_assign_mixed(&mut self, other: &Self::Affine) {
        *self = self.add_mixed(other);
    }

    fn negate(&mut self) {
        *self = self.neg();
    }

    fn mul_assign<S: Into<<Self::Scalar as PrimeField>::Repr>>(&mut self, other: S) {
        *self = self.multiply(&other.into());
    }

    fn into_affine(&self) -> Self::Affine {
        (*self).into()
    }

    fn recommended_wnaf_for_scalar(scalar: <Self::Scalar as PrimeField>::Repr) -> usize {
        let num_bits = scalar.num_bits() as usize;

        if num_bits >= 130 {
            4
        } else if num_bits >= 34 {
            3
        } else {
            2
        }
    }

    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        const RECOMMENDATIONS: [usize; 12] =
            [1, 3, 7, 20, 43, 120, 273, 563, 1630, 3128, 7933, 62569];

        let mut ret = 4;
        for r in &RECOMMENDATIONS {
            if num_scalars > *r {
                ret += 1;
            } else {
                break;
            }
        }

        ret
    }

    fn hash(_msg: &[u8]) -> Self {
        unimplemented!("not supported");
    }
}

#[derive(Clone, Debug)]
pub struct G2Prepared(pub(crate) Vec<blst_fp6>);

impl crate::PairingCurveAffine for G2Affine {
    type Prepared = G2Prepared;
    type Pair = G1Affine;
    type PairingResult = Fp12;

    fn prepare(&self) -> Self::Prepared {
        let mut lines = vec![blst_fp6::default(); 68];
        unsafe { blst_precompute_lines(lines.as_mut_ptr(), &self.0) }
        G2Prepared(lines)
    }

    fn pairing_with(&self, other: &Self::Pair) -> Self::PairingResult {
        self.perform_pairing(other)
    }
}

#[derive(Copy, Clone)]
pub struct G2Uncompressed([u8; 192]);

encoded_point_delegations!(G2Uncompressed);

impl fmt::Debug for G2Uncompressed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        self.0[..].fmt(formatter)
    }
}

impl groupy::EncodedPoint for G2Uncompressed {
    type Affine = G2Affine;

    fn empty() -> Self {
        G2Uncompressed([0; 192])
    }
    fn size() -> usize {
        192
    }
    fn into_affine(&self) -> Result<G2Affine, groupy::GroupDecodingError> {
        G2Affine::from_uncompressed(&self.0).ok_or(groupy::GroupDecodingError::NotInSubgroup)
    }

    fn into_affine_unchecked(&self) -> Result<G2Affine, groupy::GroupDecodingError> {
        G2Affine::from_uncompressed_unchecked(&self.0)
            .ok_or(groupy::GroupDecodingError::UnexpectedInformation)
    }

    fn from_affine(affine: G2Affine) -> Self {
        Self(affine.to_uncompressed())
    }
}

#[derive(Copy, Clone)]
pub struct G2Compressed([u8; 96]);

encoded_point_delegations!(G2Compressed);

impl groupy::EncodedPoint for G2Compressed {
    type Affine = G2Affine;

    fn empty() -> Self {
        G2Compressed([0; 96])
    }
    fn size() -> usize {
        96
    }
    fn into_affine(&self) -> Result<G2Affine, groupy::GroupDecodingError> {
        G2Affine::from_compressed(&self.0).ok_or(groupy::GroupDecodingError::NotInSubgroup)
    }

    fn into_affine_unchecked(&self) -> Result<G2Affine, groupy::GroupDecodingError> {
        G2Affine::from_compressed_unchecked(&self.0)
            .ok_or(groupy::GroupDecodingError::UnexpectedInformation)
    }

    fn from_affine(affine: G2Affine) -> Self {
        G2Compressed(affine.to_compressed())
    }
}

impl fmt::Debug for G2Compressed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        self.0[..].fmt(formatter)
    }
}

#[cfg(test)]
mod tests {
    use crate::{Fp, Fp2, FpRepr, G2Affine, G2Projective};
    use fff::{Field, PrimeField};
    use groupy::CurveProjective;

    #[test]
    fn g2_test_is_valid() {
        // Reject point on isomorphic twist (b = 3 * (u + 1))
        {
            let p = G2Affine::from_raw_unchecked(
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0xa757072d9fa35ba9,
                        0xae3fb2fb418f6e8a,
                        0xc1598ec46faa0c7c,
                        0x7a17a004747e3dbe,
                        0xcc65406a7c2e5a73,
                        0x10b8c03d64db4d0c,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0xd30e70fe2f029778,
                        0xda30772df0f5212e,
                        0x5b47a9ff9a233a50,
                        0xfb777e5b9b568608,
                        0x789bac1fec71a2b9,
                        0x1342f02e2da54405,
                    ]))
                    .unwrap(),
                ),
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0xfe0812043de54dca,
                        0xe455171a3d47a646,
                        0xa493f36bc20be98a,
                        0x663015d9410eb608,
                        0x78e82a79d829a544,
                        0x40a00545bb3c1e,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x4709802348e79377,
                        0xb5ac4dc9204bcfbd,
                        0xda361c97d02f42b2,
                        0x15008b1dc399e8df,
                        0x68128fd0548a3829,
                        0x16a613db5c873aaa,
                    ]))
                    .unwrap(),
                ),
                false,
            );
            assert!(!p.is_on_curve());
        }

        // Reject point on a twist (b = 2 * (u + 1))
        {
            let p = G2Affine::from_raw_unchecked(
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0xf4fdfe95a705f917,
                        0xc2914df688233238,
                        0x37c6b12cca35a34b,
                        0x41abba710d6c692c,
                        0xffcc4b2b62ce8484,
                        0x6993ec01b8934ed,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0xb94e92d5f874e26,
                        0x44516408bc115d95,
                        0xe93946b290caa591,
                        0xa5a0c2b7131f3555,
                        0x83800965822367e7,
                        0x10cf1d3ad8d90bfa,
                    ]))
                    .unwrap(),
                ),
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0xbf00334c79701d97,
                        0x4fe714f9ff204f9a,
                        0xab70b28002f3d825,
                        0x5a9171720e73eb51,
                        0x38eb4fd8d658adb7,
                        0xb649051bbc1164d,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x9225814253d7df75,
                        0xc196c2513477f887,
                        0xe05e2fbd15a804e0,
                        0x55f2b8efad953e04,
                        0x7379345eda55265e,
                        0x377f2e6208fd4cb,
                    ]))
                    .unwrap(),
                ),
                false,
            );
            assert!(!p.is_on_curve());
            assert!(!p.is_torsion_free());
        }

        // Reject point in an invalid subgroup
        // There is only one r-order subgroup, as r does not divide the cofactor.
        {
            let p = G2Affine::from_raw_unchecked(
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0x262cea73ea1906c,
                        0x2f08540770fabd6,
                        0x4ceb92d0a76057be,
                        0x2199bc19c48c393d,
                        0x4a151b732a6075bf,
                        0x17762a3b9108c4a7,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x26f461e944bbd3d1,
                        0x298f3189a9cf6ed6,
                        0x74328ad8bc2aa150,
                        0x7e147f3f9e6e241,
                        0x72a9b63583963fff,
                        0x158b0083c000462,
                    ]))
                    .unwrap(),
                ),
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0x91fb0b225ecf103b,
                        0x55d42edc1dc46ba0,
                        0x43939b11997b1943,
                        0x68cad19430706b4d,
                        0x3ccfb97b924dcea8,
                        0x1660f93434588f8d,
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0xaaed3985b6dcb9c7,
                        0xc1e985d6d898d9f4,
                        0x618bd2ac3271ac42,
                        0x3940a2dbb914b529,
                        0xbeb88137cf34f3e7,
                        0x1699ee577c61b694,
                    ]))
                    .unwrap(),
                ),
                false,
            );
            assert!(p.is_on_curve());
            assert!(!p.is_torsion_free());
        }
    }

    #[test]
    fn test_g2_addition_correctness() {
        let mut p = G2Projective::from_raw_unchecked(
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0x6c994cc1e303094e,
                    0xf034642d2c9e85bd,
                    0x275094f1352123a9,
                    0x72556c999f3707ac,
                    0x4617f2e6774e9711,
                    0x100b2fe5bffe030b,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0x7a33555977ec608,
                    0xe23039d1fe9c0881,
                    0x19ce4678aed4fcb5,
                    0x4637c4f417667e2e,
                    0x93ebe7c3e41f6acc,
                    0xde884f89a9a371b,
                ]))
                .unwrap(),
            ),
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0xe073119472e1eb62,
                    0x44fb3391fe3c9c30,
                    0xaa9b066d74694006,
                    0x25fd427b4122f231,
                    0xd83112aace35cae,
                    0x191b2432407cbb7f,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0xf68ae82fe97662f5,
                    0xe986057068b50b7d,
                    0x96c30f0411590b48,
                    0x9eaa6d19de569196,
                    0xf6a03d31e2ec2183,
                    0x3bdafaf7ca9b39b,
                ]))
                .unwrap(),
            ),
            Fp2::one(),
        );

        p.add_assign(&G2Projective::from_raw_unchecked(
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0xa8c763d25910bdd3,
                    0x408777b30ca3add4,
                    0x6115fcc12e2769e,
                    0x8e73a96b329ad190,
                    0x27c546f75ee1f3ab,
                    0xa33d27add5e7e82,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0x93b1ebcd54870dfe,
                    0xf1578300e1342e11,
                    0x8270dca3a912407b,
                    0x2089faf462438296,
                    0x828e5848cd48ea66,
                    0x141ecbac1deb038b,
                ]))
                .unwrap(),
            ),
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0xf5d2c28857229c3f,
                    0x8c1574228757ca23,
                    0xe8d8102175f5dc19,
                    0x2767032fc37cc31d,
                    0xd5ee2aba84fd10fe,
                    0x16576ccd3dd0a4e8,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0x4da9b6f6a96d1dd2,
                    0x9657f7da77f1650e,
                    0xbc150712f9ffe6da,
                    0x31898db63f87363a,
                    0xabab040ddbd097cc,
                    0x11ad236b9ba02990,
                ]))
                .unwrap(),
            ),
            Fp2::one(),
        ));

        let p = G2Affine::from(p);

        assert_eq!(
            p,
            G2Affine::from_raw_unchecked(
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0xcde7ee8a3f2ac8af,
                        0xfc642eb35975b069,
                        0xa7de72b7dd0e64b7,
                        0xf1273e6406eef9cc,
                        0xababd760ff05cb92,
                        0xd7c20456617e89
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0xd1a50b8572cbd2b8,
                        0x238f0ac6119d07df,
                        0x4dbe924fe5fd6ac2,
                        0x8b203284c51edf6b,
                        0xc8a0b730bbb21f5e,
                        0x1a3b59d29a31274
                    ]))
                    .unwrap(),
                ),
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0x9e709e78a8eaa4c9,
                        0xd30921c93ec342f4,
                        0x6d1ef332486f5e34,
                        0x64528ab3863633dc,
                        0x159384333d7cba97,
                        0x4cb84741f3cafe8
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x242af0dc3640e1a4,
                        0xe90a73ad65c66919,
                        0x2bd7ca7f4346f9ec,
                        0x38528f92b689644d,
                        0xb6884deec59fb21f,
                        0x3c075d3ec52ba90
                    ]))
                    .unwrap(),
                ),
                false,
            )
        );
    }

    #[test]
    fn test_g2_doubling_correctness() {
        let mut p = G2Projective::from_raw_unchecked(
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0x6c994cc1e303094e,
                    0xf034642d2c9e85bd,
                    0x275094f1352123a9,
                    0x72556c999f3707ac,
                    0x4617f2e6774e9711,
                    0x100b2fe5bffe030b,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0x7a33555977ec608,
                    0xe23039d1fe9c0881,
                    0x19ce4678aed4fcb5,
                    0x4637c4f417667e2e,
                    0x93ebe7c3e41f6acc,
                    0xde884f89a9a371b,
                ]))
                .unwrap(),
            ),
            Fp2::new(
                Fp::from_repr(FpRepr::new([
                    0xe073119472e1eb62,
                    0x44fb3391fe3c9c30,
                    0xaa9b066d74694006,
                    0x25fd427b4122f231,
                    0xd83112aace35cae,
                    0x191b2432407cbb7f,
                ]))
                .unwrap(),
                Fp::from_repr(FpRepr::new([
                    0xf68ae82fe97662f5,
                    0xe986057068b50b7d,
                    0x96c30f0411590b48,
                    0x9eaa6d19de569196,
                    0xf6a03d31e2ec2183,
                    0x3bdafaf7ca9b39b,
                ]))
                .unwrap(),
            ),
            Fp2::one(),
        );

        p.double();

        let p = G2Affine::from(p);

        assert_eq!(
            p,
            G2Affine::from_raw_unchecked(
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0x91ccb1292727c404,
                        0x91a6cb182438fad7,
                        0x116aee59434de902,
                        0xbcedcfce1e52d986,
                        0x9755d4a3926e9862,
                        0x18bab73760fd8024
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x4e7c5e0a2ae5b99e,
                        0x96e582a27f028961,
                        0xc74d1cf4ef2d5926,
                        0xeb0cf5e610ef4fe7,
                        0x7b4c2bae8db6e70b,
                        0xf136e43909fca0
                    ]))
                    .unwrap(),
                ),
                Fp2::new(
                    Fp::from_repr(FpRepr::new([
                        0x954d4466ab13e58,
                        0x3ee42eec614cf890,
                        0x853bb1d28877577e,
                        0xa5a2a51f7fde787b,
                        0x8b92866bc6384188,
                        0x81a53fe531d64ef
                    ]))
                    .unwrap(),
                    Fp::from_repr(FpRepr::new([
                        0x4c5d607666239b34,
                        0xeddb5f48304d14b3,
                        0x337167ee6e8e3cb6,
                        0xb271f52f12ead742,
                        0x244e6c2015c83348,
                        0x19e2deae6eb9b441
                    ]))
                    .unwrap(),
                ),
                false,
            )
        );
    }

    #[test]
    fn groupy_g2_curve_tests() {
        use groupy::tests::curve_tests;
        curve_tests::<G2Projective>();
    }
}
