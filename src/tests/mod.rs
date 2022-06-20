pub mod engine;
pub mod field;

#[cfg(feature = "gpu")]
#[test]
fn u64_to_u32_test() {
    use rand_core::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let seed = [0; 16];
    let mut rng = XorShiftRng::from_seed(seed);

    let u64_limbs: Vec<u64> = (0..6).map(|_| rng.next_u64()).collect();

    let u32_limbs = crate::u64_to_u32(&u64_limbs);

    let u64_le_bytes: Vec<u8> = u64_limbs
        .iter()
        .flat_map(|limb| limb.to_le_bytes())
        .collect();
    let u32_le_bytes: Vec<u8> = u32_limbs
        .iter()
        .flat_map(|limb| limb.to_le_bytes())
        .collect();

    assert_eq!(u64_le_bytes, u32_le_bytes);
}
