//! Encoding and decoding of the (23, 12, 7) standard Golay code.

use binfield_matrix::{matrix_mul, matrix_mul_systematic};

/// Encode the given 12 data bits into a 23-bit codeword.
pub fn encode(data: u16) -> u32 {
    assert_eq!(data >> 12, 0);

    // Compute wG.
    matrix_mul_systematic(data, CORE_XPOSE)
}

/// Try to decode the given 23-bit word to the nearest codeword, correcting up to 3
/// errors and detecting 4 errors.
///
/// If decoding was successful, return `Some((data, err))`, where `data` is the 12
/// data bits and `err` is the number of corrected bits. Otherwise, return `None` to
/// indicate an unrecoverable error.
pub fn decode(word: u32) -> Option<(u16, usize)> {
    assert_eq!(word >> 23, 0);

    // Strip off parity bits.
    let data = (word >> 11) as u16;

    // Check for 1 to 3 errors isolated in the parity bits.
    let s: u16 = matrix_mul(word, PAR);
    let n = s.count_ones() as usize;

    if n <= 3 {
        return Some((data, n));
    }

    // Check for cases with 1 error in the data bits and 0 to 2 errors in the parity bits.
    for (i, &syn) in SYN.iter().enumerate() {
        let n = (s ^ syn).count_ones() as usize;

        if n <= 2 {
            return Some((data ^ 1 << i, n + 1));
        }
    }

    // Check for 2 or 3 errors isolated to the data bits (except data MSB).
    let s: u16 = matrix_mul(rotate_11(word), PAR);
    let n = s.count_ones() as usize;

    if n <= 3 {
        return Some((data ^ s, n));
    }

    // Check for cases with 2 or 3 errors in the data bits (one being the data MSB) and
    // possibly 1 error in the parity bits or 2 errors in the data bits (exluding data
    // MSB) and 1 error in the parity bits.
    for (i, &syn) in SYN.iter().enumerate() {
        let r = s ^ syn;
        let n = r.count_ones() as usize;

        if n <= 2 {
            // The 0th syndrome corresponds to the data MSB (which now lies in the upper
            // 12 bits of the rotated word), so it must be flipped, but the following
            // syndromes correspond to parity bits, which don't need to be flipped.
            return if i == 0 {
                Some((data ^ r ^ 1 << 11, n + 1))
            } else {
                Some((data ^ r, 3))
            };
        }
    }

    None
}

/// Circularly shift the given 23-bit word right by 11 bits.
fn rotate_11(word: u32) -> u32 {
    let parity = word & 0x7FF;
    word >> 11 | parity << 12
}

/// Transpose of generator parity submatrix with extended code's LSB parity bit removed,
/// also known as **A**<sup>T</sup>.
const CORE_XPOSE: &[u16] = &[
    0b101001001111,
    0b111101101000,
    0b011110110100,
    0b001111011010,
    0b000111101101,
    0b101010111001,
    0b111100010011,
    0b110111000110,
    0b011011100011,
    0b100100111110,
    0b010010011111,
];

/// Parity-check matrix, also known as **H** = [ **A**<sup>T</sup> | I ].
///
/// This is derived from the standard-form generator matrix in the standard way.
///
/// Note that the top 12 rows of **H**<sup>T</sup> equal **A**.
const PAR: &[u32] = &[
    0b10100100111110000000000,
    0b11110110100001000000000,
    0b01111011010000100000000,
    0b00111101101000010000000,
    0b00011110110100001000000,
    0b10101011100100000100000,
    0b11110001001100000010000,
    0b11011100011000000001000,
    0b01101110001100000000100,
    0b10010011111000000000010,
    0b01001001111100000000001,
];

/// Syndromes for each single-bit error in the upper 12 bits.
///
/// The first syndrome corresponds to an LSB (bit 12) error, and the last corresponds to
/// an MSB (bit 23) error.
const SYN: &[u16] = &[
    0b10001110101,
    0b10010011111,
    0b10101001011,
    0b11011100011,
    0b00110110011,
    0b01101100110,
    0b11011001100,
    0b00111101101,
    0b01111011010,
    0b11110110100,
    0b01100011101,
    0b11000111010,
];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_rotate11() {
        assert_eq!(rotate_11(0b111111111111_00000000000), 0b000000000001_11111111111);
        assert_eq!(rotate_11(0b000000000000_11111111111), 0b111111111110_00000000000);
        assert_eq!(rotate_11(0b100000000000_00000000000), 0b000000000001_00000000000);
    }

    #[test]
    fn test_encode() {
        assert_eq!(encode(0), 0);
        assert_eq!(encode(0b111111111111), 0b111111111111_11111111111);
        assert_eq!(encode(0b111111000000), 0b111111000000_11001101000);
        assert_eq!(encode(0b000000111111), 0b000000111111_00110010111);
        assert_eq!(encode(0b100000000001), 0b100000000001_01001001111);
    }

    #[test]
    fn test_decode() {
        let w = 0b101010101010;
        let e = encode(w);
        assert_eq!(e, 0b1010101010_1000101111001);

        assert_eq!(decode(e^0b01000000000000000000010), Some((w, 2)));
        assert_eq!(decode(e^0b00100000000000000000100), Some((w, 2)));
        assert_eq!(decode(e^0b00010000000000000001000), Some((w, 2)));
        assert_eq!(decode(e^0b00001000000000000010000), Some((w, 2)));
        assert_eq!(decode(e^0b00000100000000000100000), Some((w, 2)));
        assert_eq!(decode(e^0b00000010000000001000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000001000000010000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000000100000100000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000000010001000000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000000001010000000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000000010000000000001), Some((w, 2)));
        assert_eq!(decode(e^0b00000000100000000000010), Some((w, 2)));
        assert_eq!(decode(e^0b00000001000000000000100), Some((w, 2)));
        assert_eq!(decode(e^0b00000010000000000001000), Some((w, 2)));
        assert_eq!(decode(e^0b00000100000000000010000), Some((w, 2)));
        assert_eq!(decode(e^0b00001000000000000100000), Some((w, 2)));
        assert_eq!(decode(e^0b00010000000000001000000), Some((w, 2)));
        assert_eq!(decode(e^0b00100000000000010000000), Some((w, 2)));
        assert_eq!(decode(e^0b01000000000000100000000), Some((w, 2)));
        assert_eq!(decode(e^0b10000000000001000000000), Some((w, 2)));
        assert_eq!(decode(e^0b10000000000010000000000), Some((w, 2)));
        assert_eq!(decode(e^0b11100000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b01110000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00111000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00011100000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00001110000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000111000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000011100000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000001110000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000111000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000011100000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000001110000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000111000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000011100000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000001110000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000111000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000011100000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000001110000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000000111000), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000000011100), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000000001110), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000000000111), Some((w, 3)));
        assert_eq!(decode(e^0b00000000000000000000000), Some((w, 0)));
        assert_eq!(decode(e^0b00000000000000000000001), Some((w, 1)));
        assert_eq!(decode(e^0b00000000000000000000011), Some((w, 2)));
        assert_eq!(decode(e^0b00000000000000000000111), Some((w, 3)));
        assert_eq!(decode(e^0b00000001000000000000000), Some((w, 1)));
        assert_eq!(decode(e^0b00000011000000000000000), Some((w, 2)));
        assert_eq!(decode(e^0b00000111000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b00000100000000000000001), Some((w, 2)));
        assert_eq!(decode(e^0b00000110000000000000001), Some((w, 3)));
        assert_eq!(decode(e^0b00000100000000000000011), Some((w, 3)));

        // Exhaustively test that all codewords are decoded correctly.
        for w in 0..1<<12 {
            assert_eq!(decode(encode(w)), Some((w, 0)));
        }

        let w = encode(0b110111101110);

        // Exhaustively test that all 0 through 3-bit errors are detected.
        for ((i, j), k) in (0..23).zip(0..23).zip(0..23) {
            let e: u32 = 1 << i | 1 << j | 1 << k;
            let n = e.count_ones() as usize;

            assert_eq!(decode(w ^ e), Some((0b110111101110, n)));
        }

        // Exhaustively test that all 4-bit errors are detected.
        for (((h, i), j), k) in (0..23).zip(0..23).zip(0..23).zip(0..23) {
            let e: u32 = 1 << h | 1 << i | 1 << j | 1 << k;
            let n = e.count_ones() as usize;

            if n >= 4 {
                assert_eq!(decode(w ^ e), None);
            }
        }
    }
}
