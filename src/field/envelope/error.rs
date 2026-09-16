use std::fmt;

/// Invalid in-memory envelope data. Axis indices refer to (x, y, z, xi).
#[derive(Debug, PartialEq)]
pub enum GridError {
    Axis {
        index: Option<usize>,
        reason: &'static str,
    },
    Wavelength,
    SizeOverflow {
        shape: [usize; 4],
    },
    Shape {
        expected: usize,
        actual: usize,
    },
    Value {
        index: usize,
        value: f64,
    },
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Axis { index, reason } => match index {
                Some(i) => write!(
                    f,
                    "invalid {} axis: {}",
                    ["x", "y", "z", "xi"].get(*i).unwrap_or(&"unknown"),
                    reason
                ),
                None => write!(f, "invalid axis: {}", reason),
            },
            Self::Wavelength => write!(
                f,
                "reference wavelength must be finite and positive, in metres"
            ),
            Self::SizeOverflow { shape } => {
                write!(f, "grid shape {:?} exceeds addressable f64 storage", shape)
            }
            Self::Shape { expected, actual } => write!(
                f,
                "grid requires {} scalar values, got {}",
                expected, actual
            ),
            Self::Value { index, value } => write!(
                f,
                "scalar value at flat index {} must be finite and nonnegative, got {}",
                index, value
            ),
        }
    }
}

impl std::error::Error for GridError {}

/// Sampling failure. Stored axis indices refer to (x, y, z, xi);
/// laboratory component indices refer to (ct, x, y, z).
#[derive(Debug, PartialEq)]
pub enum SampleError {
    NonFiniteCoordinate {
        axis: usize,
        coordinate: f64,
    },
    OutOfDomain {
        axis: usize,
        coordinate: f64,
        min: f64,
        max: f64,
    },
    NonFiniteLabCoordinate {
        component: usize,
        coordinate: f64,
    },
    NonFiniteXi {
        ct: f64,
        z: f64,
    },
    TransformedOutOfDomain {
        position: [f64; 4],
        axis: usize,
        coordinate: f64,
        min: f64,
        max: f64,
    },
    NonFiniteResult {
        quantity: &'static str,
        value: f64,
    },
}

impl fmt::Display for SampleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let axis_name = |axis: usize| {
            ["x", "y", "z", "xi"]
                .get(axis)
                .copied()
                .unwrap_or("unknown")
        };
        match self {
            Self::NonFiniteLabCoordinate { component, coordinate } => write!(f,
                "nonfinite laboratory {} coordinate: {} (expected finite metres)",
                ["ct", "x", "y", "z"].get(*component).unwrap_or(&"unknown"), coordinate),
            Self::NonFiniteXi { ct, z } => write!(f,
                "nonfinite transformed xi = ct - z from finite laboratory ct={} m, z={} m", ct, z),
            Self::TransformedOutOfDomain { position, axis, coordinate, min, max } => write!(f,
                "laboratory (ct,x,y,z)={:?} m maps via xi=ct-z to stored {}={} m outside canonical domain [{}, {}] m",
                position, axis_name(*axis), coordinate, min, max),
            Self::NonFiniteResult { quantity, value } => write!(f,
                "cannot return laboratory envelope sample: nonfinite {} result ({})", quantity, value),
            Self::NonFiniteCoordinate { axis, coordinate } => write!(
                f,
                "nonfinite {} coordinate: {} (expected finite metres)",
                axis_name(*axis),
                coordinate
            ),
            Self::OutOfDomain {
                axis,
                coordinate,
                min,
                max,
            } => write!(
                f,
                "{} coordinate {} m is outside canonical domain [{}, {}] m",
                axis_name(*axis),
                coordinate,
                min,
                max
            ),
        }
    }
}

impl std::error::Error for SampleError {}
