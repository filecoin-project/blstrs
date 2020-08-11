//! # `blstrs`
//!
//! An implementation of the BLS12-381 pairing-friendly elliptic curve construction.

#[macro_use]
mod macros;

pub mod scalar;

pub use scalar::Scalar;
