//! Validated scalar envelopes on (x, y, z, xi), where xi = ct - z.
//!
//! All coordinates are metres; values are dimensionless S = a_rms^2, with
//! polarization averaging already applied. This module does not implement Field.
//! Sampling returns stored-coordinate derivatives (S_x, S_y, S_z, S_xi),
//! in inverse metres. Laboratory sampling returns the raised (+---) four-gradient.

// The grid is intentionally not connected to a simulation field yet.
#![allow(dead_code)]

mod cubic;
mod error;
mod grid;
mod laboratory;
mod multilinear;

#[cfg(test)]
mod validation;

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

/// Explicit interpolation choice. Existing sampling entry points use Multilinear.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InterpolationMethod {
    Multilinear,
    /// Shared-slope Hermite interpolation; signed overshoot is permitted.
    Cubic,
}

impl EnvelopeGrid {
    /// Sample S and stored derivatives with an explicit method.
    /// Cubic requires >=3 nodes per axis and is diagnostic until positivity is addressed.
    pub fn sample_stored_with_method(
        &self,
        q: [f64; 4],
        method: InterpolationMethod,
    ) -> Result<StoredEnvelopeSample, SampleError> {
        match method {
            InterpolationMethod::Multilinear => self.sample_stored(q),
            InterpolationMethod::Cubic => cubic::sample(self, q),
        }
    }
}
