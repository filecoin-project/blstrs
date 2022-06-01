//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

#![deny(clippy::all, clippy::perf, clippy::correctness)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::wrong_self_convention)]

#[cfg(not(target_endian = "little"))]
compile_error!("blstrs is only supported on little endian architectures");

#[macro_use]
mod macros;

mod fp;
mod fp12;
mod fp2;
mod fp6;
mod g1;
mod g2;
mod gt;
mod pairing;
mod scalar;
mod traits;

pub use g1::{G1Affine, G1Compressed, G1Projective, G1Uncompressed};
pub use g2::{G2Affine, G2Compressed, G2Prepared, G2Projective, G2Uncompressed};
pub use gt::Gt;
pub use pairing::*;
pub use scalar::Scalar;
pub use traits::Compress;

#[cfg(feature = "serde")]
mod serde_impl;

#[cfg(test)]
mod tests;

// export for benchmarking only
#[cfg(feature = "__private_bench")]
pub use crate::{fp::Fp, fp12::Fp12, fp2::Fp2};

use ff::Field;
use group::prime::PrimeCurveAffine;
use pairing_lib::{Engine, MultiMillerLoop, PairingCurveAffine};

/// Bls12-381 engine
#[derive(Debug, Copy, Clone)]
pub struct Bls12;

impl Engine for Bls12 {
    type Fr = Scalar;
    type G1 = G1Projective;
    type G1Affine = G1Affine;
    type G2 = G2Projective;
    type G2Affine = G2Affine;
    type Gt = Gt;

    fn pairing(p: &Self::G1Affine, q: &Self::G2Affine) -> Self::Gt {
        pairing(p, q)
    }
}

impl MultiMillerLoop for Bls12 {
    type G2Prepared = G2Prepared;
    type Result = MillerLoopResult;

    /// Computes $$\sum_{i=1}^n \textbf{ML}(a_i, b_i)$$ given a series of terms
    /// $$(a_1, b_1), (a_2, b_2), ..., (a_n, b_n).$$
    fn multi_miller_loop(terms: &[(&Self::G1Affine, &Self::G2Prepared)]) -> Self::Result {
        let mut res = blst::blst_fp12::default();

        for (i, (p, q)) in terms.iter().enumerate() {
            let mut tmp = blst::blst_fp12::default();
            if (p.is_identity() | q.is_identity()).into() {
                // Define pairing with zero as one, matching what `pairing` does.
                tmp = crate::fp12::Fp12::one().0;
            } else {
                unsafe {
                    blst::blst_miller_loop_lines(&mut tmp, q.lines.as_ptr(), &p.0);
                }
            }
            if i == 0 {
                res = tmp;
            } else {
                unsafe {
                    blst::blst_fp12_mul(&mut res, &res, &tmp);
                }
            }
        }

        MillerLoopResult(crate::fp12::Fp12(res))
    }
}

#[cfg(feature = "gpu")]
fn u64_to_u32(limbs: &[u64]) -> Vec<u32> {
    limbs
        .iter()
        .flat_map(|limb| vec![(limb & 0xFFFF_FFFF) as u32, (limb >> 32) as u32])
        .collect()
}

#[test]
fn bls12_engine_tests() {
    crate::tests::engine::engine_tests::<Bls12>();
}
