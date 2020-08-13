//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

#[macro_use]
mod macros;

mod scalar;
mod g1;
mod g2;
mod fp;
mod fp2;
mod fp6;
mod fp12;

pub use scalar::Scalar;
pub use g1::{G1Affine, G1Projective};
pub use g2::{G2Affine, G2Projective};
pub use fp::Fp;
pub use fp2::Fp2;
pub use fp6::Fp6;
pub use fp12::Fp12;


