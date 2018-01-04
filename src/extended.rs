//! Encoding and decoding of the (24, 12, 8) extended Golay code.

use binfield_matrix::{matrix_mul, matrix_mul_systematic};

/// Encode the given 12 data bits into a 24-bit codeword.
pub fn encode(data: u16) -> u32 {
    assert_eq!(data >> 12, 0);

    // Compute wG = w[ I | A ].
    matrix_mul_systematic(data, CORE_XPOSE)
}

/// Try to decode the given 24-bit word to the nearest codeword, correcting up to 3
/// errors and detecting 4 errors.
///
/// If decoding was successful, return `Some((data, err))`, where `data` is the 12
/// data bits and `err` is the number of corrected bits. Otherwise, return `None` to
/// indicate an unrecoverable error.
pub fn decode(word: u32) -> Option<(u16, usize)> {
    assert_eq!(word >> 24, 0);

    // Strip off parity bits.
    let data = (word >> 12) as u16;

    // Compute wG<sup>T</sup> to check for errors isolated to upper 12 bits.
    let s: u16 = matrix_mul(word, PAR_ALT);
    let n = s.count_ones() as usize;

    if n <= 3 {
        return Some((data ^ s, n));
    }

    // Check for cases with one error in lower 12 bits and one or two errors in upper
    // 12 bits [2. p261].
    //
    // For each e<sub>i</sub> = 1 << i, compute (w + e<sub>i</sub>)G<sup>T</sup> =
    // wG<sup>T</sup> + e<sub>i</sub>G<sup>T</sup> = s + a<sub>i</sub>, where
    // a<sub>i</sub> is the i'th row from the bottom of A</sup>T</sup>.
    //
    // Since e<sub>i<sub> isn't used to repair the data bits, we instead just loop
    // over all the words in A<sup>T</sup>.
    for &q in CORE_XPOSE.iter() {
        let syn = s ^ q;
        let n = syn.count_ones() as usize;

        if n <= 2 {
            return Some((data ^ syn, n + 1));
        }
    }

    // Compute wH<sup>T</sup> to check for errors isolated to lower 12 bits.
    let s: u16 = matrix_mul(word, PAR);
    let n = s.count_ones() as usize;

    if n <= 3 {
        return Some((data, n));
    }

    // Check for cases with one error in upper 12 bits and 2 errors in lower 12 bits [2,
    // p261].
    //
    // For each e<sub>i</sub>, compute (w + e<sub>i+12</sub>)H<sup>T</sup> =
    // wH<sup>T</sup> + e<sub>i+12</sub>H<sup>T</sup> = s + b<sub>i</sub>, where
    // b<sub>i</sub> is the (i+12)'th row from the bottom of H<sup>T</sup>, which
    // equals the i'th row from the bottom of A.
    for (i, &q) in CORE.iter().enumerate() {
        let syn = s ^ q;

        if syn.count_ones() <= 2 {
            let err = 1 << 11 >> i;
            return Some((data ^ err, 3));
        }
    }

    None
}

/// Generator parity submatrix, also known as **A**.
const CORE: &[u16] = &[
    0b110001110101,
    0b011000111011,
    0b111101101000,
    0b011110110100,
    0b001111011010,
    0b110110011001,
    0b011011001101,
    0b001101100111,
    0b110111000110,
    0b101010010111,
    0b100100111110,
    0b100011101011,
];

/// Transpose of generator parity submatrix, also known as **A**<sup>T</sup>.
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
    0b110001110101,
];

/// Parity-check matrix, also known as **H** = [ **A**<sup>T</sup> | I ].
///
/// This is derived from the standard-form generator matrix in the standard way.
///
/// Note that the top 12 rows of **H**<sup>T</sup> equal **A**.
const PAR: &[u32] = &[
    0b101001001111100000000000,
    0b111101101000010000000000,
    0b011110110100001000000000,
    0b001111011010000100000000,
    0b000111101101000010000000,
    0b101010111001000001000000,
    0b111100010011000000100000,
    0b110111000110000000010000,
    0b011011100011000000001000,
    0b100100111110000000000100,
    0b010010011111000000000010,
    0b110001110101000000000001,
];

/// Alternative parity-check matrix, which equals the standard-form generator matrix **G**
/// = **H**<sup>⊥</sup> = [ **I** | **A** ].
///
/// Note that the bottom 12 rows of **G**<sup>T</sup> equal **A**<sup>T</sup>.
///
/// Since the Golay code is self-dual, this can also be used as a parity-check matrix [2,
/// p258].
const PAR_ALT: &[u32] = &[
    0b100000000000110001110101,
    0b010000000000011000111011,
    0b001000000000111101101000,
    0b000100000000011110110100,
    0b000010000000001111011010,
    0b000001000000110110011001,
    0b000000100000011011001101,
    0b000000010000001101100111,
    0b000000001000110111000110,
    0b000000000100101010010111,
    0b000000000010100100111110,
    0b000000000001100011101011,
];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_encode() {
        assert_eq!(encode(0), 0);
        assert_eq!(encode(0b111111111111), 0b111111111111_111111111111);
        assert_eq!(encode(0b111111000000), 0b111111000000_110011010001);
        assert_eq!(encode(0b000000111111), 0b000000111111_001100101110);
        assert_eq!(encode(0b100000000001), 0b100000000001_010010011110);
    }

    #[test]
    fn test_decode() {
        let w = 0b111111101010;
        let e = encode(w);
        assert_eq!(e, 0b111111101010_111011100100);

        assert_eq!(decode(e^0b100000000000000000000010), Some((w, 2)));
        assert_eq!(decode(e^0b010000000000000000000001), Some((w, 2)));
        assert_eq!(decode(e^0b001000000000000000000010), Some((w, 2)));
        assert_eq!(decode(e^0b000100000000000000000100), Some((w, 2)));
        assert_eq!(decode(e^0b000010000000000000001000), Some((w, 2)));
        assert_eq!(decode(e^0b000001000000000000010000), Some((w, 2)));
        assert_eq!(decode(e^0b000000100000000000100000), Some((w, 2)));
        assert_eq!(decode(e^0b000000010000000001000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000001000000010000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000000100000100000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000000010001000000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000000001010000000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000000010000000000001), Some((w, 2)));
        assert_eq!(decode(e^0b000000000100000000000010), Some((w, 2)));
        assert_eq!(decode(e^0b000000001000000000000100), Some((w, 2)));
        assert_eq!(decode(e^0b000000010000000000001000), Some((w, 2)));
        assert_eq!(decode(e^0b000000100000000000010000), Some((w, 2)));
        assert_eq!(decode(e^0b000001000000000000100000), Some((w, 2)));
        assert_eq!(decode(e^0b000010000000000001000000), Some((w, 2)));
        assert_eq!(decode(e^0b000100000000000010000000), Some((w, 2)));
        assert_eq!(decode(e^0b001000000000000100000000), Some((w, 2)));
        assert_eq!(decode(e^0b010000000000001000000000), Some((w, 2)));
        assert_eq!(decode(e^0b010000000000010000000000), Some((w, 2)));
        assert_eq!(decode(e^0b111000000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b011100000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b001110000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000111000000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000011100000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000001110000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000111000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000011100000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000001110000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000111000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000011100000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000001110000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000111000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000011100000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000001110000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000111000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000011100000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000001110000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000000111000), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000000011100), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000000001110), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000000000111), Some((w, 3)));
        assert_eq!(decode(e^0b000000000000000000000000), Some((w, 0)));
        assert_eq!(decode(e^0b000000000000000000000001), Some((w, 1)));
        assert_eq!(decode(e^0b000000000000000000000011), Some((w, 2)));
        assert_eq!(decode(e^0b000000000000000000000111), Some((w, 3)));
        assert_eq!(decode(e^0b000000001000000000000000), Some((w, 1)));
        assert_eq!(decode(e^0b000000011000000000000000), Some((w, 2)));
        assert_eq!(decode(e^0b000000111000000000000000), Some((w, 3)));
        assert_eq!(decode(e^0b000000100000000000000001), Some((w, 2)));
        assert_eq!(decode(e^0b000000110000000000000001), Some((w, 3)));
        assert_eq!(decode(e^0b000000100000000000000011), Some((w, 3)));

        // Exhaustively test that all codewords are decoded correctly.
        for w in 0..1<<12 {
            assert_eq!(decode(encode(w)), Some((w, 0)));
        }

        let w = encode(0b110110100110);

        // Exhaustively test that all 0 through 3-bit errors are detected.
        for ((i, j), k) in (0..24).zip(0..24).zip(0..24) {
            let e: u32 = 1 << i | 1 << j | 1 << k;
            let n = e.count_ones() as usize;

            assert_eq!(decode(w ^ e), Some((0b110110100110, n)));
        }

        // Exhaustively test that all 4-bit errors are detected.
        for (((h, i), j), k) in (0..24).zip(0..24).zip(0..24).zip(0..24) {
            let e: u32 = 1 << h | 1 << i | 1 << j | 1 << k;
            let n = e.count_ones() as usize;

            if n >= 4 {
                assert_eq!(decode(w ^ e), None);
            }
        }
    }
}
