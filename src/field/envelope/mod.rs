//! Validated scalar envelopes on (x, y, z, xi), where xi = ct - z.
//!
//! All coordinates are metres; values are dimensionless S = a_rms^2, with
//! polarization averaging already applied. This module does not implement Field.
//! Sampling returns stored-coordinate derivatives (S_x, S_y, S_z, S_xi),
//! in inverse metres. Laboratory sampling returns the raised (+---) four-gradient.

// The grid is intentionally not connected to a simulation field yet.
#![allow(dead_code)]

mod error;
mod grid;
mod laboratory;
mod multilinear;

pub use error::GridError;
// Public entry points are reserved for the next interpolation/integration steps.
#[allow(unused_imports)]
pub use error::SampleError;
#[allow(unused_imports)]
pub use grid::{EnvelopeGrid, EnvelopeMetadata, UniformAxis};
#[allow(unused_imports)]
pub use laboratory::LabEnvelopeSample;
#[allow(unused_imports)]
pub use multilinear::StoredEnvelopeSample;
