use super::GridError;
use crate::field::Polarization;

/// Uniform, increasing coordinates in metres. At least two nodes are required.
#[derive(Clone, Copy, Debug)]
pub struct UniformAxis {
    origin: f64,
    spacing: f64,
    len: usize,
}

impl UniformAxis {
    /// Validate supplied nodes, then store the uniform grid defined by its endpoints.
    /// Residual tolerance is 1e-10 of a cell plus eight coordinate-roundoff units.
    /// Axes with roundoff allowance above 1e-3 of a cell are rejected.
    /// Supplied nodes are not retained. Reconstruction can round even the last
    /// endpoint; sampling bounds are the canonical coordinates returned
    /// by `coordinate(0)` and `coordinate(len - 1)`, not the supplied endpoints.
    pub fn try_from_coordinates(nodes: &[f64]) -> Result<Self, GridError> {
        let invalid = |reason| GridError::Axis {
            index: None,
            reason,
        };
        if nodes.len() < 2 {
            return Err(invalid("at least two nodes are required"));
        }
        if nodes.iter().any(|x| !x.is_finite()) {
            return Err(invalid("coordinates must be finite"));
        }
        if nodes.windows(2).any(|w| w[1] <= w[0]) {
            return Err(invalid("coordinates must be strictly increasing"));
        }
        let origin = nodes[0];
        let end = nodes[nodes.len() - 1];
        let spacing = (end - origin) / (nodes.len() - 1) as f64;
        if !spacing.is_finite() || spacing <= 0.0 || !spacing.recip().is_finite() {
            return Err(invalid(
                "spacing and inverse spacing must be finite and positive",
            ));
        }
        let roundoff = 8.0 * f64::EPSILON * origin.abs().max(end.abs());
        if roundoff > 1.0e-3 * spacing {
            return Err(invalid("origin is too large to resolve spacing reliably"));
        }
        let tolerance = 1.0e-10 * spacing + roundoff;
        for (i, &node) in nodes.iter().enumerate() {
            let expected = (i as f64).mul_add(spacing, origin);
            if !expected.is_finite() || (node - expected).abs() > tolerance {
                return Err(invalid("coordinates must be uniformly spaced"));
            }
        }
        Ok(Self {
            origin,
            spacing,
            len: nodes.len(),
        })
    }

    pub fn origin(&self) -> f64 {
        self.origin
    }
    pub fn spacing(&self) -> f64 {
        self.spacing
    }
    pub fn len(&self) -> usize {
        self.len
    }

    /// Canonical node coordinate; invalid indices have no coordinate.
    pub fn coordinate(&self, index: usize) -> Option<f64> {
        (index < self.len).then(|| (index as f64).mul_add(self.spacing, self.origin))
    }
}

/// S already includes cycle averaging: LP S_peak = a0^2/2, CP S_peak = a0^2.
/// Polarization is descriptive here; no additional scaling is applied.
#[derive(Clone, Copy)]
pub struct EnvelopeMetadata {
    reference_wavelength_m: f64,
    polarization: Polarization,
}

impl EnvelopeMetadata {
    pub fn try_new(
        reference_wavelength_m: f64,
        polarization: Polarization,
    ) -> Result<Self, GridError> {
        if !reference_wavelength_m.is_finite() || reference_wavelength_m <= 0.0 {
            return Err(GridError::Wavelength);
        }
        Ok(Self {
            reference_wavelength_m,
            polarization,
        })
    }

    pub fn reference_wavelength_m(&self) -> f64 {
        self.reference_wavelength_m
    }
    pub fn polarization(&self) -> Polarization {
        self.polarization
    }
}

/// Owned, validated S = a_rms^2 in C order [nx, ny, nz, nxi], xi fastest.
pub struct EnvelopeGrid {
    axes: [UniformAxis; 4],
    shape: [usize; 4],
    values: Vec<f64>,
    metadata: EnvelopeMetadata,
}

fn element_count(shape: [usize; 4]) -> Result<usize, GridError> {
    let invalid = || GridError::SizeOverflow { shape };
    let count = shape
        .iter()
        .try_fold(1usize, |n, &len| n.checked_mul(len))
        .ok_or_else(invalid)?;
    let bytes = count
        .checked_mul(std::mem::size_of::<f64>())
        .ok_or_else(invalid)?;
    if bytes > isize::MAX as usize {
        return Err(invalid());
    }
    Ok(count)
}

impl EnvelopeGrid {
    /// Validate axes in (x,y,z,xi) order and consume scalar storage without copying.
    pub fn try_new(
        coordinates: [&[f64]; 4],
        values: Vec<f64>,
        metadata: EnvelopeMetadata,
    ) -> Result<Self, GridError> {
        let axis = |i: usize| {
            UniformAxis::try_from_coordinates(coordinates[i]).map_err(|error| match error {
                GridError::Axis { reason, .. } => GridError::Axis {
                    index: Some(i),
                    reason,
                },
                other => other,
            })
        };
        let axes = [axis(0)?, axis(1)?, axis(2)?, axis(3)?];
        let shape = axes.map(|axis| axis.len());
        let expected = element_count(shape)?;
        if values.len() != expected {
            return Err(GridError::Shape {
                expected,
                actual: values.len(),
            });
        }
        for (index, &value) in values.iter().enumerate() {
            if !value.is_finite() || value < 0.0 {
                return Err(GridError::Value { index, value });
            }
        }
        Ok(Self {
            axes,
            shape,
            values,
            metadata,
        })
    }

    pub fn axes(&self) -> &[UniformAxis; 4] {
        &self.axes
    }
    pub fn shape(&self) -> [usize; 4] {
        self.shape
    }
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    pub fn metadata(&self) -> &EnvelopeMetadata {
        &self.metadata
    }

    pub fn flat_index(&self, indices: [usize; 4]) -> Option<usize> {
        if indices.iter().zip(self.shape.iter()).any(|(i, n)| i >= n) {
            return None;
        }
        let [x, y, z, xi] = indices;
        Some(((x * self.shape[1] + y) * self.shape[2] + z) * self.shape[3] + xi)
    }

    pub fn node(&self, indices: [usize; 4]) -> Option<f64> {
        self.flat_index(indices).map(|index| self.values[index])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metadata(pol: Polarization) -> EnvelopeMetadata {
        EnvelopeMetadata::try_new(0.8e-6, pol).unwrap()
    }

    #[test]
    fn shifted_axes_and_indexing() {
        let x = [-3.0, -1.0];
        let y = [2.0, 2.5, 3.0];
        let z = [-6.0, -3.0, 0.0, 3.0];
        let xi = [7.0, 11.0, 15.0, 19.0, 23.0];
        let mut values = Vec::new();
        for i in 0..2 {
            for j in 0..3 {
                for k in 0..4 {
                    for l in 0..5 {
                        values.push((1000 * i + 100 * j + 10 * k + l) as f64);
                    }
                }
            }
        }
        let ptr = values.as_ptr();
        let grid = EnvelopeGrid::try_new([&x, &y, &z, &xi], values, metadata(Polarization::Linear))
            .unwrap();
        assert_eq!(grid.shape(), [2, 3, 4, 5]);
        assert_eq!(grid.values().as_ptr(), ptr);
        assert_eq!(grid.axes()[1].origin(), 2.0);
        assert_eq!(grid.axes()[1].spacing(), 0.5);
        assert_eq!(grid.axes()[3].coordinate(4), Some(23.0));
        assert_eq!(grid.axes()[3].coordinate(5), None);
        for i in 0..2 {
            for j in 0..3 {
                for k in 0..4 {
                    for l in 0..5 {
                        assert_eq!(
                            grid.node([i, j, k, l]),
                            Some((1000 * i + 100 * j + 10 * k + l) as f64)
                        );
                    }
                }
            }
        }
        for axis in 0..4 {
            let mut index = [0; 4];
            index[axis] = grid.shape()[axis];
            assert_eq!(grid.flat_index(index), None);
        }
        assert_eq!(grid.node([usize::MAX; 4]), None);
    }

    #[test]
    fn axis_validation() {
        for nodes in [
            vec![],
            vec![0.0],
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, f64::NAN],
            vec![0.0, f64::INFINITY],
            vec![0.0, 1.0, 2.01],
            vec![-f64::MAX, f64::MAX],
            vec![0.0, f64::from_bits(1)],
            vec![1e16, 1e16 + 2.0],
        ] {
            assert!(
                UniformAxis::try_from_coordinates(&nodes).is_err(),
                "{:?}",
                nodes
            );
        }
        let nodes: Vec<_> = (0..100).map(|i| -3.7e-6 + i as f64 * 0.13e-6).collect();
        assert!(UniformAxis::try_from_coordinates(&nodes).is_ok());
        assert!(UniformAxis::try_from_coordinates(&[0.0, 1.0 + 1e-12, 2.0]).is_ok());
        assert!(UniformAxis::try_from_coordinates(&[0.0, 1.0 + 1e-8, 2.0]).is_err());
    }

    #[test]
    fn rounded_input_reconstructs_canonical_coordinates() {
        // Decimal nodes as supplied by a rounded file export: 0.3 micrometre
        // spacing. Compare to these independent physical coordinates, rather
        // than recomputing the implementation's endpoint/FMA formula.
        let supplied = [
            -3.7e-6, -3.4e-6, -3.1e-6, -2.8e-6, -2.5e-6, -2.2e-6, -1.9e-6, -1.6e-6, -1.3e-6,
            -1.0e-6,
        ];
        let axis = UniformAxis::try_from_coordinates(&supplied).unwrap();
        for (i, &expected) in supplied.iter().enumerate() {
            let actual = axis.coordinate(i).unwrap();
            assert!((actual - expected).abs() < 2.0e-21);
            if i > 0 {
                let previous = axis.coordinate(i - 1).unwrap();
                assert!(actual > previous);
                assert!((actual - previous - 0.3e-6).abs() < 2.0e-21);
            }
        }
        assert_eq!(axis.coordinate(0), Some(supplied[0]));
        // This fixture specifically exercises endpoint reconstruction rounding:
        // the supplied upper endpoint lies just inside the canonical endpoint.
        let canonical_end = axis.coordinate(supplied.len() - 1).unwrap();
        assert!(canonical_end > supplied[supplied.len() - 1]);
        assert_eq!(axis.coordinate(supplied.len()), None);
    }

    #[test]
    fn invalid_grid_data() {
        let a = [0.0, 1.0];
        let m = metadata(Polarization::Linear);
        assert!(matches!(
            EnvelopeGrid::try_new([&a; 4], vec![0.0; 15], m),
            Err(GridError::Shape {
                expected: 16,
                actual: 15
            })
        ));
        assert!(matches!(
            EnvelopeGrid::try_new([&a, &a, &[1.0, 0.0], &a], vec![], m),
            Err(GridError::Axis { index: Some(2), .. })
        ));
        for bad in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut values = vec![0.0; 16];
            values[7] = bad;
            assert!(matches!(
                EnvelopeGrid::try_new([&a; 4], values, m),
                Err(GridError::Value { index: 7, .. })
            ));
        }
    }

    #[test]
    fn checked_storage_size() {
        assert_eq!(element_count([2, 3, 4, 5]), Ok(120));
        for shape in [
            [usize::MAX, 2, 2, 2],
            [usize::MAX / 8, 2, 2, 2],
            [isize::MAX as usize / 8 + 1, 1, 1, 1],
        ] {
            assert!(matches!(
                element_count(shape),
                Err(GridError::SizeOverflow { .. })
            ));
        }
    }

    #[test]
    fn metadata_and_normalization() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                EnvelopeMetadata::try_new(bad, Polarization::Linear),
                Err(GridError::Wavelength)
            ));
        }
        let a = [0.0, 1.0];
        for pol in [Polarization::Linear, Polarization::Circular] {
            let grid = EnvelopeGrid::try_new([&a; 4], vec![0.5; 16], metadata(pol)).unwrap();
            assert!(grid.metadata().polarization() == pol);
            assert_eq!(grid.metadata().reference_wavelength_m(), 0.8e-6);
            assert_eq!(grid.node([1; 4]), Some(0.5));
        }
        assert!(
            EnvelopeGrid::try_new([&a; 4], vec![0.0; 16], metadata(Polarization::Linear)).is_ok()
        );
    }

    #[test]
    fn informative_errors() {
        let error = GridError::Axis {
            index: Some(3),
            reason: "coordinates must be finite",
        };
        assert_eq!(
            error.to_string(),
            "invalid xi axis: coordinates must be finite"
        );
        assert!(GridError::Shape {
            expected: 16,
            actual: 15
        }
        .to_string()
        .contains("got 15"));
        assert!(GridError::Value {
            index: 7,
            value: -1.0
        }
        .to_string()
        .contains("index 7"));
    }
}
