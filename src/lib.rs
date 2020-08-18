//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

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
