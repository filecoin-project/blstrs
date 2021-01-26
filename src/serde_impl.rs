use std::fmt;
use std::marker::PhantomData;

use fff::PrimeField;
use groupy::{CurveAffine, CurveProjective, EncodedPoint};

use serde::de::{Error as DeserializeError, SeqAccess, Visitor};
use serde::ser::SerializeTuple;
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::{
    Fp, Fp12, Fp2, Fp6, FpRepr, G1Affine, G1Projective, G2Affine, G2Projective, Scalar, ScalarRepr,
};

const ERR_CODE: &str = "deserialized bytes don't encode a group element";

impl Serialize for G1Projective {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.into_affine().serialize(s)
    }
}

impl<'de> Deserialize<'de> for G1Projective {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(G1Affine::deserialize(d)?.into_projective())
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
        self.into_affine().serialize(s)
    }
}

impl<'de> Deserialize<'de> for G2Projective {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(G2Affine::deserialize(d)?.into_projective())
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

/// Serializes a group element using its compressed representation.
fn serialize_affine<S: Serializer, C: CurveAffine>(c: &C, s: S) -> Result<S::Ok, S::Error> {
    let len = C::Compressed::size();
    let mut tup = s.serialize_tuple(len)?;
    for byte in c.into_compressed().as_ref() {
        tup.serialize_element(byte)?;
    }
    tup.end()
}

/// Deserializes the compressed representation of a group element.
fn deserialize_affine<'de, D: Deserializer<'de>, C: CurveAffine>(d: D) -> Result<C, D::Error> {
    struct TupleVisitor<C> {
        _ph: PhantomData<C>,
    }

    impl<'de, C: CurveAffine> Visitor<'de> for TupleVisitor<C> {
        type Value = C;

        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            let len = C::Compressed::size();
            write!(f, "a tuple of size {}", len)
        }

        #[inline]
        fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<C, A::Error> {
            let mut compressed = C::Compressed::empty();
            for (i, byte) in compressed.as_mut().iter_mut().enumerate() {
                let len_err = || DeserializeError::invalid_length(i, &self);
                *byte = seq.next_element()?.ok_or_else(len_err)?;
            }
            let to_err = |_| DeserializeError::custom(ERR_CODE);
            compressed.into_affine().map_err(to_err)
        }
    }

    let len = C::Compressed::size();
    d.deserialize_tuple(len, TupleVisitor { _ph: PhantomData })
}

impl Serialize for Scalar {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.into_repr().serialize(s)
    }
}

impl<'de> Deserialize<'de> for Scalar {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Scalar::from_repr(ScalarRepr::deserialize(d)?).map_err(|_| D::Error::custom(ERR_CODE))
    }
}

impl Serialize for ScalarRepr {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(s)
    }
}

impl<'de> Deserialize<'de> for ScalarRepr {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(ScalarRepr(<_>::deserialize(d)?))
    }
}

impl Serialize for Fp {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.into_repr().serialize(s)
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
        Fp::from_repr(FpRepr::deserialize(d)?).map_err(|_| D::Error::custom(ERR_CODE))
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

impl Serialize for FpRepr {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.0.serialize(s)
    }
}

impl<'de> Deserialize<'de> for FpRepr {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Ok(FpRepr(<_>::deserialize(d)?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Fp12;

    use std::fmt::Debug;

    use fff::Field;
    use rand_core::SeedableRng;
    use rand_xorshift::XorShiftRng;

    fn test_roundtrip<T: Serialize + for<'a> Deserialize<'a> + Debug + PartialEq>(t: &T) {
        dbg!(t);
        let ser = serde_json::to_vec(t).unwrap();
        assert_eq!(*t, serde_json::from_slice(&ser).unwrap());
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
            test_roundtrip(&g.into_affine());
        }

        let g = G1Projective::zero();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_slice::<G1Projective>(&hex::decode("5b3139322c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c305d").unwrap()).unwrap(),
            g
        );
        test_roundtrip(&g.into_affine());
        assert_eq!(
            serde_json::from_slice::<G1Affine>(&hex::decode("5b3139322c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c305d").unwrap()).unwrap(),
            g.into_affine(),
        );

        let g = G1Projective::one();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_slice::<G1Projective>(&hex::decode("5b3135312c3234312c3231312c3136372c34392c3135312c3231352c3134382c33382c3134392c39392c3134302c37392c3136392c3137322c31352c3139352c3130342c3134302c37392c3135312c3131362c3138352c352c3136312c37382c35382c36332c32332c32372c3137322c38382c3130382c38352c3233322c36332c3234392c3132322c32362c3233392c3235312c35382c3234302c31302c3231392c33342c3139382c3138375d").unwrap()).unwrap(),
            g
        );
        test_roundtrip(&g.into_affine());
        assert_eq!(
            serde_json::from_slice::<G1Affine>(&hex::decode("5b3135312c3234312c3231312c3136372c34392c3135312c3231352c3134382c33382c3134392c39392c3134302c37392c3136392c3137322c31352c3139352c3130342c3134302c37392c3135312c3131362c3138352c352c3136312c37382c35382c36332c32332c32372c3137322c38382c3130382c38352c3233322c36332c3234392c3132322c32362c3233392c3235312c35382c3234302c31302c3231392c33342c3139382c3138375d").unwrap()).unwrap(),
            g.into_affine(),
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
            test_roundtrip(&g.into_affine());
        }

        let g = G2Projective::zero();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_slice::<G2Projective>(&hex::decode("5b3139322c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c305d").unwrap()).unwrap(),
            g
        );
        test_roundtrip(&g.into_affine());
        assert_eq!(
            serde_json::from_slice::<G2Affine>(&hex::decode("5b3139322c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c302c305d").unwrap()).unwrap(),
            g.into_affine(),
        );

        let g = G2Projective::one();
        test_roundtrip(&g);
        assert_eq!(
            serde_json::from_slice::<G2Projective>(&hex::decode("5b3134372c3232342c34332c39362c38322c3131332c3135392c39362c3132352c3137322c3231312c3136302c3133362c33392c37392c3130312c38392c3130372c3230382c3230382c3135332c33322c3138322c32362c3138312c3231382c39372c3138372c3232302c3132372c38302c37332c35312c37362c3234312c31382c31392c3134382c39332c38372c3232392c3137322c3132352c352c39332c342c34332c3132362c322c37342c3136322c3137382c3234302c3134332c31302c3134352c33382c382c352c33392c34352c3139372c31362c38312c3139382c3232382c3132322c3231322c3235302c36342c35392c322c3138302c38312c31312c3130302c3132322c3232372c3230392c3131392c31312c3137322c332c33382c3136382c352c3138372c3233392c3231322c3132382c38362c3230302c3139332c33332c3138392c3138345d").unwrap()).unwrap(),
            g
        );
        test_roundtrip(&g.into_affine());
        assert_eq!(
            serde_json::from_slice::<G2Affine>(&hex::decode("5b3134372c3232342c34332c39362c38322c3131332c3135392c39362c3132352c3137322c3231312c3136302c3133362c33392c37392c3130312c38392c3130372c3230382c3230382c3135332c33322c3138322c32362c3138312c3231382c39372c3138372c3232302c3132372c38302c37332c35312c37362c3234312c31382c31392c3134382c39332c38372c3232392c3137322c3132352c352c39332c342c34332c3132362c322c37342c3136322c3137382c3234302c3134332c31302c3134352c33382c382c352c33392c34352c3139372c31362c38312c3139382c3232382c3132322c3231322c3235302c36342c35392c322c3138302c38312c31312c3130302c3132322c3232372c3230392c3131392c31312c3137322c332c33382c3136382c352c3138372c3233392c3231322c3132382c38362c3230302c3139332c33332c3138392c3138345d").unwrap()).unwrap(),
            g.into_affine()
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
            test_roundtrip(&f.into_repr());
        }

        let f = Scalar::zero();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Scalar>(&hex::decode("5b302c302c302c305d").unwrap()).unwrap(),
            f
        );
        test_roundtrip(&f.into_repr());
        assert_eq!(
            serde_json::from_slice::<ScalarRepr>(&hex::decode("5b302c302c302c305d").unwrap())
                .unwrap(),
            f.into_repr(),
        );

        let f = Scalar::one();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Scalar>(&hex::decode("5b312c302c302c305d").unwrap()).unwrap(),
            f
        );
        test_roundtrip(&f.into_repr());
        assert_eq!(
            serde_json::from_slice::<ScalarRepr>(&hex::decode("5b312c302c302c305d").unwrap())
                .unwrap(),
            f.into_repr(),
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
            test_roundtrip(&f.into_repr());
        }

        let f = Fp::zero();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp>(&hex::decode("5b302c302c302c302c302c305d").unwrap())
                .unwrap(),
            f
        );
        test_roundtrip(&f.into_repr());
        assert_eq!(
            serde_json::from_slice::<FpRepr>(&hex::decode("5b302c302c302c302c302c305d").unwrap())
                .unwrap(),
            f.into_repr(),
        );

        let f = Fp::one();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp>(&hex::decode("5b312c302c302c302c302c305d").unwrap())
                .unwrap(),
            f
        );
        test_roundtrip(&f.into_repr());
        assert_eq!(
            serde_json::from_slice::<FpRepr>(&hex::decode("5b312c302c302c302c302c305d").unwrap())
                .unwrap(),
            f.into_repr(),
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

        let f = Fp12::zero();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp12>(&hex::decode("7b226330223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d2c226331223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d7d").unwrap()).unwrap(),
            f
        );

        let f = Fp12::one();
        test_roundtrip(&f);
        assert_eq!(
            serde_json::from_slice::<Fp12>(&hex::decode("7b226330223a7b226330223a7b226330223a5b312c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d2c226331223a7b226330223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226331223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d2c226332223a7b226330223a5b302c302c302c302c302c305d2c226331223a5b302c302c302c302c302c305d7d7d7d").unwrap()).unwrap(),
            f
        );
    }
}
