pub fn decode_hex_into_slice(buffer: &mut [u8], bytes: &[u8]) {
    debug_assert_eq!(buffer.len(), bytes.len() / 2);
    let mut i = 0;
    while i < buffer.len() {
        buffer[i] = decode_hex_byte([bytes[2 * i], bytes[2 * i + 1]]);
        i += 1;
    }
}

/// Decode a single byte encoded as two hexadecimal characters.
pub fn decode_hex_byte(bytes: [u8; 2]) -> u8 {
    let mut i = 0;
    let mut result = 0u8;

    while i < 2 {
        result <<= 4;
        result |= match bytes[i] {
            b @ b'0'..=b'9' => b - b'0',
            b @ b'a'..=b'f' => 10 + b - b'a',
            b @ b'A'..=b'F' => 10 + b - b'A',
            b => {
                assert!(
                    matches!(b, b'0'..=b'9' | b'a' ..= b'f' | b'A'..=b'F'),
                    "invalid hex byte"
                );
                0
            }
        };

        i += 1;
    }

    result
}
