use crate::{Fp12, G1Affine, G2Affine};

use blst::*;

/// Execute a complete pairing operation `(p, q)`.
pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    let mut tmp = blst_fp12::default();
    unsafe { blst_miller_loop(&mut tmp, &q.0, &p.0) };

    let mut out = blst_fp12::default();
    unsafe { blst_final_exp(&mut out, &tmp) };

    out.into()
}
