//! Validated scalar envelopes on (x, y, z, xi), where xi = ct - z.
//!
//! All coordinates are metres; values are dimensionless S = a_rms^2, with
//! polarization averaging already applied. This module does not implement Field.
//! Future sampling must return the raised laboratory gradient
//! (S_xi, -S_x, -S_y, -S_z + S_xi), in inverse metres.

// The grid is intentionally not connected to a simulation field yet.
#![allow(dead_code)]

mod error;
mod grid;

pub use error::GridError;
// Public entry points are reserved for the next interpolation/integration steps.
#[allow(unused_imports)]
pub use grid::{EnvelopeGrid, EnvelopeMetadata, UniformAxis};
