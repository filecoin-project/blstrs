use crate::{Fp12, G1Affine, G2Affine, G2Prepared};

use blst::*;

/// Execute a complete pairing operation `(p, q)`.
pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    let mut tmp = blst_fp12::default();
    unsafe { blst_miller_loop(&mut tmp, &q.0, &p.0) };

    let mut out = blst_fp12::default();
    unsafe { blst_final_exp(&mut out, &tmp) };

    out.into()
}

/// Execute a complete pairing operation on many elements.
pub fn pairing_many(pairs: &[(G1Affine, G2Prepared)]) -> Fp12 {
    // FIXME: use miller_loop_n when available
    let mut tmp = blst_fp12::default();

    todo!();
    // unsafe { blst_miller_loop(&mut tmp, &q.0, &p.0) };

    let mut out = blst_fp12::default();
    unsafe { blst_final_exp(&mut out, &tmp) };

    out.into()
}
