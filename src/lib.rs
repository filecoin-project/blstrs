//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

#![deny(clippy::all, clippy::perf, clippy::correctness)]
#![allow(clippy::many_single_char_names)]

#[macro_use]
mod macros;

mod fp;
mod fp12;
mod fp2;
mod fp6;
mod g1;
mod g2;
mod pairing;
mod scalar;
mod traits;

pub use fff::*;
pub use fp::{Fp, FpRepr};
pub use fp12::Fp12;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use g1::{G1Affine, G1Projective};
pub use g2::{G2Affine, G2Prepared, G2Projective};
pub use pairing::*;
pub use scalar::{Scalar, ScalarRepr, S as SCALAR_S};
pub use traits::*;

#[cfg(test)]
mod tests;

/// Bls12-381 engine
#[derive(Debug, Copy, Clone)]
pub struct Bls12;

impl fff::ScalarEngine for Bls12 {
    type Fr = Scalar;
}

impl Engine for Bls12 {
    type G1 = G1Projective;
    type G1Affine = G1Affine;
    type G2 = G2Projective;
    type G2Affine = G2Affine;
    type Fq = Fp;
    type Fqe = Fp2;
    type Fqk = Fp12;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<
            Item = &'a (
                &'a <Self::G1Affine as PairingCurveAffine>::Prepared,
                &'a <Self::G2Affine as PairingCurveAffine>::Prepared,
            ),
        >,
    {
        let mut ret = blst::blst_fp12::default();
        for (p, q) in i.into_iter() {
            let mut tmp = blst::blst_fp12::default();
            unsafe {
                blst::blst_miller_loop_lines(&mut tmp, q.0.as_ptr(), &p.0);
                blst::blst_fp12_mul(&mut ret, &ret, &tmp);
            }
        }
        ret.into()
    }

    fn final_exponentiation(r: &Fp12) -> Option<Fp12> {
        let mut out = blst::blst_fp12::default();
        unsafe { blst::blst_final_exp(&mut out, &r.0) };

        // TODO: What about the None case?
        Some(out.into())
    }
}

#[test]
fn bls12_engine_tests() {
    crate::tests::engine::engine_tests::<Bls12>();
}
