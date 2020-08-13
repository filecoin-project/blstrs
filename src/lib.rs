//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

#[macro_use]
mod macros;

mod scalar;
mod g1;
mod fp;

pub use scalar::Scalar;
pub use g1::{G1Affine, G1Projective};
pub use fp::Fp;
