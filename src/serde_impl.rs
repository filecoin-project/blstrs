use core::convert::TryFrom;
use core::fmt;
use core::marker::PhantomData;
use std::fmt::Formatter;

use group::{prime::PrimeCurveAffine, Curve, GroupEncoding};
use serde::{
    de::{Error as DeserializeError, SeqAccess, Visitor},
    ser::SerializeTuple,
    Deserialize, Deserializer, Serialize, Serializer,
};

use crate::gt::GtRepr;
use crate::util::decode_hex_into_slice;
use crate::{
    fp::Fp, fp12::Fp12, fp2::Fp2, fp6::Fp6, util, G1Affine, G1Projective, G2Affine, G2Projective,
    Gt, MillerLoopResult, Scalar,
};

const ERR_CODE: &str = "deserialized bytes don't encode a group element";

impl Serialize for G1Projective {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.to_affine().serialize(s)
    }
}

impl<'de> Deserialize<'de> for G1Projective {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(G1Affine::deserialize(d)?.to_curve())
    }
}

impl Serialize for G1Affine {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        serialize_affine(self, s)
    }
}

impl<'de> Deserialize<'de> for G1Affine {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        deserialize_affine(d)
    }
}

impl Serialize for G2Projective {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.to_affine().serialize(s)
    }
}

impl<'de> Deserialize<'de> for G2Projective {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(G2Affine::deserialize(d)?.to_curve())
    }
}

impl Serialize for G2Affine {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        serialize_affine(self, s)
    }
}

impl<'de> Deserialize<'de> for G2Affine {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        deserialize_affine(d)
    }
}

/// Serialize a group element using its compressed representation.
fn serialize_affine<S: Serializer, C: PrimeCurveAffine + fmt::LowerHex>(
    c: &C,
    s: S,
) -> Result<S::Ok, S::Error> {
    if s.is_human_readable() {
        s.serialize_str(&format!("{:x}", c))
    } else {
        let compressed_bytes = c.to_bytes();
        let len = compressed_bytes.as_ref().len();
        let mut tup = s.serialize_tuple(len)?;
        for byte in compressed_bytes.as_ref() {
            tup.serialize_element(byte)?;
        }
        tup.end()
    }
}

/// Deserializes the compressed representation of a group element.
fn deserialize_affine<'de, D: Deserializer<'de>, C: PrimeCurveAffine>(d: D) -> Result<C, D::Error> {
    struct TupleVisitor<C> {
        _ph: PhantomData<C>,
    }

    impl<'de, C: PrimeCurveAffine> Visitor<'de> for TupleVisitor<C> {
        type Value = C;

        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let len = C::Repr::default().as_ref().len();
            write!(f, "a tuple of {} bytes", len)
        }

        #[inline]
        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<C, A::Error> {
            let mut compressed = C::Repr::default();
            for (i, byte) in compressed.as_mut().iter_mut().enumerate() {
                let len_err = || DeserializeError::invalid_length(i, &self);
                *byte = seq.next_element()?.ok_or_else(len_err)?;
            }
            let opt = C::from_bytes(&compressed);
            if opt.is_some().into() {
                Ok(opt.unwrap())
            } else {
                Err(DeserializeError::custom(ERR_CODE))
            }
        }
    }

    if d.is_human_readable() {
        let hex_str = <&str>::deserialize(d)?;
        let mut compressed = C::Repr::default();
        let writer = compressed.as_mut();
        if hex_str.len() / 2 != writer.len() {
            return Err(DeserializeError::custom(format!(
                "invalid length, expected {}, received {}",
                writer.len(),
                hex_str.len() / 2
            )));
        }
        util::decode_hex_into_slice(writer, hex_str.as_bytes());
        let opt = C::from_bytes(&compressed);
        if opt.is_some().into() {
            Ok(opt.unwrap())
        } else {
            Err(DeserializeError::custom(ERR_CODE))
        }
    } else {
        let len = C::Repr::default().as_ref().len();
        d.deserialize_tuple(len, TupleVisitor { _ph: PhantomData })
    }
}

impl Serialize for Scalar {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        if s.is_human_readable() {
            s.serialize_str(&format!("{:x}", self))
        } else {
            let bytes = self.to_bytes_be();
            let mut tup = s.serialize_tuple(bytes.len())?;
            for byte in bytes.as_ref() {
                tup.serialize_element(byte)?;
            }
            tup.end()
        }
    }
}

impl<'de> Deserialize<'de> for Scalar {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let buf = if d.is_human_readable() {
            let hex_str = <&str>::deserialize(d)?;
            let mut bytes = [0u8; 32];
            util::decode_hex_into_slice(&mut bytes, hex_str.as_bytes());
            bytes
        } else {
            <[u8; 32]>::deserialize(d)?
        };
        Option::<Scalar>::from(Scalar::from_bytes_be(&buf))
            .ok_or_else(|| D::Error::custom(ERR_CODE))
    }
}

impl Serialize for Fp {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        let bytes = self.to_bytes_le();
        let u64s = [
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[0..8]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[8..16]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[16..24]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[24..32]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[32..40]).unwrap()),
            u64::from_le_bytes(<[u8; 8]>::try_from(&bytes[40..48]).unwrap()),
        ];
        u64s.serialize(s)
    }
}

impl<'de> Deserialize<'de> for MillerLoopResult {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let fp12 = Fp12::deserialize(d)?;
        Ok(MillerLoopResult(fp12))
    }
}

impl Serialize for MillerLoopResult {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(s)
    }
}

impl Serialize for Gt {
    fn serialize<S>(&self, s: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        if s.is_human_readable() {
            format!("{:x}", self).serialize(s)
        } else {
            let repr = self.to_bytes();
            let mut tupler = s.serialize_tuple(repr.0.len())?;
            for byte in repr.as_ref() {
                tupler.serialize_element(byte)?;
            }
            tupler.end()
        }
    }
}

impl<'de> Deserialize<'de> for Gt {
    fn deserialize<D>(d: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        if d.is_human_readable() {
            let repr = <String>::deserialize(d)?;
            let mut bytes = [0u8; 288];
            decode_hex_into_slice(&mut bytes, repr.as_bytes());
            Option::<Gt>::from(Gt::from_bytes(&GtRepr(bytes)))
                .ok_or_else(|| serde::de::Error::custom("invalid compressed value"))
        } else {
            struct ArrayVisitor;

            impl<'de> Visitor<'de> for ArrayVisitor {
                type Value = Gt;

                fn expecting(&self, f: &mut Formatter) -> fmt::Result {
                    write!(f, "an array of 288 bytes")
                }

                fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
                where
                    A: SeqAccess<'de>,
                {
                    let mut repr = GtRepr::default();
                    for i in repr.as_mut() {
                        *i = seq
                            .next_element()?
                            .ok_or_else(|| serde::de::Error::custom("invalid compressed value"))?;
                    }
                    Option::<Gt>::from(Gt::from_bytes(&repr))
                        .ok_or_else(|| serde::de::Error::custom("invalid compressed value"))
                }
            }

            d.deserialize_tuple(288, ArrayVisitor)
        }
    }
}

#[derive(Serialize, Deserialize)]
struct Fp2Ser {
    c0: Fp,
    c1: Fp,
}

#[derive(Serialize, Deserialize)]
struct Fp6Ser {
    c0: Fp2,
    c1: Fp2,
    c2: Fp2,
}

#[derive(Serialize, Deserialize)]
struct Fp12Ser {
    c0: Fp6,
    c1: Fp6,
}

impl Serialize for Fp2 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        Fp2Ser {
            c0: self.c0(),
            c1: self.c1(),
        }
        .serialize(s)
    }
}

impl Serialize for Fp12 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        Fp12Ser {
            c0: self.c0(),
            c1: self.c1(),
        }
        .serialize(s)
    }
}

impl Serialize for Fp6 {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        Fp6Ser {
            c0: self.c0(),
            c1: self.c1(),
            c2: self.c2(),
        }
        .serialize(s)
    }
}

impl<'de> Deserialize<'de> for Fp {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let deser = <[u64; 6]>::deserialize(d)?;
        match Fp::from_u64s_le(&deser).into() {
            Some(fp) => Ok(fp),
            None => Err(D::Error::custom(ERR_CODE)),
        }
    }
}

impl<'de> Deserialize<'de> for Fp2 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let Fp2Ser { c0, c1 } = Fp2Ser::deserialize(d)?;
        Ok(Fp2::new(c0, c1))
    }
}

impl<'de> Deserialize<'de> for Fp6 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let Fp6Ser { c0, c1, c2 } = Fp6Ser::deserialize(d)?;
        Ok(Fp6::new(c0, c1, c2))
    }
}

impl<'de> Deserialize<'de> for Fp12 {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let Fp12Ser { c0, c1 } = Fp12Ser::deserialize(d)?;
        Ok(Fp12::new(c0, c1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use core::fmt::Debug;

    use crate::pairing;
    use ff::Field;
    use group::{Curve, Group};
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    fn test_roundtrip<T: Serialize + for<'a> Deserialize<'a> + Debug + PartialEq>(t: &T) {
        // dbg!(t);
        let ser = serde_json::to_vec(t).unwrap();
        assert_eq!(*t, serde_json::from_slice(&ser).unwrap());
        let ser = serde_bare::to_vec(t).unwrap();
        assert_eq!(*t, serde_bare::from_slice(&ser).unwrap());
    }

    #[test]
    fn serde_g1() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..100 {
            let g = G1Projective::random(&mut rng);
            test_roundtrip(&g);
            test_roundtrip(&g.to_affine());
        }

        let g = G1Projective::identity();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_str::<G1Projective>(&"\"C00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\"").unwrap(),
            g
        );
        test_roundtrip(&g.to_affine());
        assert_eq!(
            serde_json::from_str::<G1Affine>(&"\"C00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\"").unwrap(),
            g.to_affine(),
        );

        let g = G1Projective::GENERATOR;
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_str::<G1Projective>(&"\"97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb\"").unwrap(),
            g
        );
        test_roundtrip(&g.to_affine());
        assert_eq!(
            serde_json::from_str::<G1Affine>(&"\"97f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb\"").unwrap(),
            g.to_affine(),
        );
    }

    #[test]
    fn serde_g2() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..100 {
            let g = G2Projective::random(&mut rng);
            test_roundtrip(&g);
            test_roundtrip(&g.to_affine());
        }

        let g = G2Projective::identity();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_str::<G2Projective>(&"\"c00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\"").unwrap(),
            g
        );
        test_roundtrip(&g.to_affine());
        assert_eq!(
            serde_json::from_str::<G2Affine>(&"\"c00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000\"").unwrap(),
            g.to_affine(),
        );

        let g = G2Projective::GENERATOR;
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_str::<G2Projective>(&"\"93e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8\"").unwrap(),
            g
        );
        test_roundtrip(&g.to_affine());
        assert_eq!(
            serde_json::from_str::<G2Affine>(&"\"93e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8\"").unwrap(),
            g.to_affine()
        );
    }

    #[test]
    fn serde_scalar() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        for _ in 0..100 {
            let f = Scalar::random(&mut rng);
            test_roundtrip(&f);
        }

        let f = Scalar::ZERO;
        test_roundtrip(&f);
        // The hex string "5b302c302c302c305d" encodes the unicode string "[0,0,0,0]" where each
        // byte in the hex string encodes a unicode character: 0x58 = "[", 0x30 = "0", 0x2c = ",",
        // and 0x5d = "]".
        assert_eq!(
            serde_json::from_str::<Scalar>(
                "\"0000000000000000000000000000000000000000000000000000000000000000\""
            )
            .unwrap(),
            f
        );

        let f = Scalar::ONE;
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_str::<Scalar>(
                "\"0000000000000000000000000000000000000000000000000000000000000001\""
            )
            .unwrap(),
            f
        );
    }

    #[test]
    fn serde_fp() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..100 {
            let f = Fp::random(&mut rng);
            test_roundtrip(&f);
        }

        let f = Fp::ZERO;
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp>(&hex::decode("5b302c302c302c302c302c305d").unwrap())
                .unwrap(),
            f
        );

        let f = Fp::ONE;
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp>(&hex::decode("5b312c302c302c302c302c305d").unwrap())
                .unwrap(),
            f
        );
    }

    #[test]
    fn serde_fp12() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        for _ in 0..100 {
            let f = Fp12::random(&mut rng);
            test_roundtrip(&f);
        }

        let f = Fp12::ZERO;
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp12>(&hex::decode("7b226330223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d2c226331223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d7d").unwrap()).unwrap(),
            f
        );

        let f = Fp12::ONE;
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp12>(&hex::decode("7b226330223a7b226330223a7b226330223a5b312c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d2c226331223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d7d").unwrap()).unwrap(),
            f
        );
    }

    #[test]
    fn serde_compatible_bls12_381() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        let s1 = Scalar::random(&mut rng);
        let bytes = s1.to_bytes_be();
        let s2 = bls12_381_plus::Scalar::from_be_bytes(&bytes).unwrap();
        assert_eq!(s1.to_bytes_le(), s2.to_le_bytes());

        assert_eq!(format!("{:x}", s1), format!("{:x}", s2));
        let ss1 = serde_json::to_string(&s1).unwrap();
        let s2 = serde_json::from_str::<bls12_381_plus::Scalar>(&ss1).unwrap();
        assert_eq!(s1.to_bytes_le(), s2.to_le_bytes());

        let p1 = G1Projective::random(&mut rng);
        let s1 = serde_json::to_string(&p1).unwrap();
        let p2 = serde_json::from_str::<bls12_381_plus::G1Projective>(&s1).unwrap();
        assert_eq!(
            p1.to_affine().to_compressed(),
            p2.to_affine().to_compressed()
        );

        let p1 = G2Projective::random(&mut rng);
        let s2 = serde_json::to_string(&p1).unwrap();
        let p2 = serde_json::from_str::<bls12_381_plus::G2Projective>(&s2).unwrap();
        assert_eq!(
            p1.to_affine().to_compressed(),
            p2.to_affine().to_compressed()
        );
    }

    #[test]
    fn serde_gt() {
        let gt = pairing(&G1Affine::generator(), &G2Affine::generator());
        let json = serde_json::to_string(&gt).unwrap();
        let gt2 = serde_json::from_str(&json);
        assert!(gt2.is_ok());
        assert_eq!(gt, gt2.unwrap());

        let bare = serde_bare::to_vec(&gt).unwrap();
        let gt2 = serde_bare::from_slice(&bare);
        assert!(gt2.is_ok());
        assert_eq!(gt, gt2.unwrap());
    }
}
