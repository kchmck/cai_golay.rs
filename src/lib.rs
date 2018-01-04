//! Encoding and decoding of the (23, 12, 7) standard and (24, 12, 8) extended Golay
//! codes.
//!
//! The generator polynomial for these codes is given by
//!
//! > g(x) = x<sup>11</sup> + x<sup>10</sup> + x<sup>6</sup> + x<sup>5</sup> +
//! x<sup>4</sup> + x<sup>2</sup> + 1
//!
//! The decoding algorithm for the extended code is derived from [1, p79] and [2, p261],
//! and the decoding algorithm for the standard code is derived from [3].
//!
//! ## References
//!
//! 1. *Coding Theory and Cryptography: The Essentials*, Hankerson et al, 2000
//! 2. *Coding and Information Theory*, Roman, 1992
//! 3. "High-Speed Decoding of the Binary Golay Code", Lee et al, 2013

extern crate binfield_matrix;

pub mod extended;
pub mod standard;
