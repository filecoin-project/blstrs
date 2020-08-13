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

pub use fp::Fp;
pub use fp12::Fp12;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use g1::{G1Affine, G1Projective};
pub use g2::{G2Affine, G2Prepared, G2Projective};
pub use pairing::*;
pub use scalar::{Scalar, S as SCALAR_S};
