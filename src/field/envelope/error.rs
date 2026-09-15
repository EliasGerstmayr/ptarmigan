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
