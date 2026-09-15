use super::{EnvelopeGrid, SampleError, UniformAxis};

/// One sample of the multilinear interpolant. This is not a laboratory four-vector.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoredEnvelopeSample {
    /// Dimensionless, already cycle-averaged S = a_rms^2.
    pub a_sqd: f64,
    /// [S_x, S_y, S_z, S_xi], in inverse metres; other stored coordinates fixed.
    pub derivatives: [f64; 4],
}

/// Locate using comparisons to canonical nodes, never a coordinate tolerance.
fn locate(axis: &UniformAxis, q: f64, axis_index: usize) -> Result<(usize, f64), SampleError> {
    if !q.is_finite() {
        return Err(SampleError::NonFiniteCoordinate {
            axis: axis_index,
            coordinate: q,
        });
    }
    let min = axis.coordinate(0).unwrap();
    let max = axis.coordinate(axis.len() - 1).unwrap();
    if q < min || q > max {
        return Err(SampleError::OutOfDomain {
            axis: axis_index,
            coordinate: q,
            min,
            max,
        });
    }
    if q == max {
        return Ok((axis.len() - 2, 1.0));
    }

    let mut cell = (((q - min) / axis.spacing()) as usize).min(axis.len() - 2);
    // Division may put an exact knot on the wrong side of an integer. Also
    // handle neighbouring representable queries without snapping them to knots.
    while cell > 0 && q < axis.coordinate(cell).unwrap() {
        cell -= 1;
    }
    while cell < axis.len() - 2 && q >= axis.coordinate(cell + 1).unwrap() {
        cell += 1;
    }
    let lower = axis.coordinate(cell).unwrap();
    let fraction = (q - lower) / axis.spacing();
    // The query has already passed strict canonical bounds and cell comparisons.
    // Canonical subtraction/division can round a fraction just above one. Bound
    // only this arithmetic result; never expand the domain or change the cell.
    Ok((cell, fraction.min(1.0)))
}

impl EnvelopeGrid {
    /// Sample at stored (x,y,z,xi), all in metres, with xi = ct - z.
    ///
    /// The canonical domain is closed. Interior knots use the positive-side cell;
    /// the upper endpoint uses the last cell at fraction one. Values are continuous
    /// to floating-point accuracy; derivatives can jump at cell boundaries.
    /// No heap allocation, field integration or laboratory conversion is performed.
    pub fn sample_stored(&self, point: [f64; 4]) -> Result<StoredEnvelopeSample, SampleError> {
        let mut cells = [0; 4];
        let mut fractions = [0.0; 4];
        for axis in 0..4 {
            let (cell, fraction) = locate(&self.axes()[axis], point[axis], axis)?;
            cells[axis] = cell;
            fractions[axis] = fraction;
        }

        // Corner bit a selects the upper node along stored axis a. The backing
        // storage remains xi-fastest; only these 16 scalar values are on the stack.
        let mut corners = [0.0; 16];
        for (bits, value) in corners.iter_mut().enumerate() {
            let indices = std::array::from_fn(|axis| cells[axis] + ((bits >> axis) & 1));
            *value = self.node(indices).unwrap();
        }

        // Pair lower/upper derivative terms before weighting: this is exactly
        // the +/-1/spacing rule, but cancels a constant offset before division.
        let mut derivatives = [0.0; 4];
        for axis in 0..4 {
            let bit = 1 << axis;
            for lower in 0..16 {
                if lower & bit != 0 {
                    continue;
                }
                let mut weight = 1.0;
                for other in 0..4 {
                    if other != axis {
                        weight *= if lower & (1 << other) == 0 {
                            1.0 - fractions[other]
                        } else {
                            fractions[other]
                        };
                    }
                }
                // Skip zero-weight edges, avoiding 0 * overflow for extreme data.
                if weight != 0.0 {
                    derivatives[axis] += weight
                        * ((corners[lower | bit] - corners[lower]) / self.axes()[axis].spacing());
                }
            }
        }

        // Nested convex combinations are the 16-corner product-weight sum, with
        // fewer operations and without accumulating a large common offset 16 times.
        let mut count = 16;
        for u in fractions {
            for i in 0..count / 2 {
                corners[i] = (1.0 - u) * corners[2 * i] + u * corners[2 * i + 1];
            }
            count /= 2;
        }
        Ok(StoredEnvelopeSample {
            a_sqd: corners[0],
            derivatives,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{envelope::EnvelopeMetadata, Polarization};

    // Float neighbours without next_up/next_down, preserving the crate's MSRV.
    fn next_up(x: f64) -> f64 {
        if x == 0.0 {
            f64::from_bits(1)
        } else {
            f64::from_bits(if x > 0.0 {
                x.to_bits() + 1
            } else {
                x.to_bits() - 1
            })
        }
    }
    fn next_down(x: f64) -> f64 {
        -next_up(-x)
    }

    fn grid(nodes: [&[f64]; 4], f: impl Fn([f64; 4]) -> f64) -> EnvelopeGrid {
        let axes = nodes.map(|n| UniformAxis::try_from_coordinates(n).unwrap());
        let mut values = Vec::new();
        for x in 0..axes[0].len() {
            for y in 0..axes[1].len() {
                for z in 0..axes[2].len() {
                    for xi in 0..axes[3].len() {
                        let indices = [x, y, z, xi];
                        values.push(f(std::array::from_fn(|a| {
                            axes[a].coordinate(indices[a]).unwrap()
                        })));
                    }
                }
            }
        }
        EnvelopeGrid::try_new(
            nodes,
            values,
            EnvelopeMetadata::try_new(0.8e-6, Polarization::Linear).unwrap(),
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64, scale: f64) {
        // These fixtures have O(1) values and well-resolved cells. 128 eps allows
        // accumulated corner arithmetic/subtraction error, scaled by S or S/h.
        let tolerance = 128.0 * f64::EPSILON * scale;
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={:e}, expected={:e}, tolerance={:e}",
            actual,
            expected,
            tolerance
        );
    }

    #[test]
    fn constant_value_and_zero_derivatives() {
        let a = [-1.0, 0.0, 1.0];
        let g = grid([&a; 4], |_| 3.25);
        for q in [[-0.71, 0.32, -0.13, 0.91], [-1.0; 4], [1.0; 4], [0.0; 4]] {
            let s = g.sample_stored(q).unwrap();
            close(s.a_sqd, 3.25, 3.25);
            assert_eq!(s.derivatives, [0.0; 4]); // paired edges cancel exactly
        }
    }

    #[test]
    fn affine_shifted_unequal_axes() {
        let nodes: [&[f64]; 4] = [
            &[-2.0, 0.0],
            &[1.0, 1.5, 2.0],
            &[-3.0, -2.0, -1.0, 0.0],
            &[2.0, 2.25, 2.5, 2.75, 3.0],
        ];
        let coeff = [2.0, 3.0, 5.0, 7.0];
        let f = |q: [f64; 4]| 30.0 + (0..4).map(|a| coeff[a] * q[a]).sum::<f64>();
        let g = grid(nodes, f);
        for q in [[-1.3, 1.22, -1.8, 2.31], [-0.1, 1.83, -0.13, 2.91]] {
            let s = g.sample_stored(q).unwrap();
            close(s.a_sqd, f(q), 60.0);
            for a in 0..4 {
                close(s.derivatives[a], coeff[a], 60.0 / g.axes()[a].spacing());
            }
        }
    }

    #[test]
    fn positive_multiaffine_product() {
        let a = [-0.5, 0.0, 0.5, 1.0];
        let b = [0.2, 0.3, 0.4, 0.5];
        let f = |q: [f64; 4]| (0..4).map(|i| 2.0 + b[i] * q[i]).product::<f64>();
        let g = grid([&a; 4], f);
        for q in [[-0.31, 0.12, 0.38, 0.87], [0.51, -0.27, 0.93, 0.11]] {
            let s = g.sample_stored(q).unwrap();
            close(s.a_sqd, f(q), 40.0);
            for i in 0..4 {
                let expected = b[i]
                    * (0..4)
                        .filter(|&j| j != i)
                        .map(|j| 2.0 + b[j] * q[j])
                        .product::<f64>();
                close(s.derivatives[i], expected, 40.0 / 0.5);
            }
        }
    }

    #[test]
    fn realistic_micrometre_affine() {
        let a = [-3.7e-6, -3.4e-6, -3.1e-6, -2.8e-6];
        let coeff = [2e5, -3e5, 5e5, 7e5];
        let f = |q: [f64; 4]| 10.0 + (0..4).map(|i| coeff[i] * q[i]).sum::<f64>();
        let g = grid([&a; 4], f);
        let q = [-3.61e-6, -3.22e-6, -2.91e-6, -3.47e-6];
        let s = g.sample_stored(q).unwrap();
        close(s.a_sqd, f(q), 10.0);
        for i in 0..4 {
            close(s.derivatives[i], coeff[i], 10.0 / 0.3e-6);
        }
    }

    #[test]
    fn xi_and_stored_z_are_independent() {
        let a = [-1e-6, 0.0, 1e-6];
        for varying in [2, 3] {
            let g = grid([&a; 4], |q| 2.0 + 4e5 * q[varying]);
            let s = g
                .sample_stored([0.17e-6, -0.29e-6, 0.33e-6, -0.41e-6])
                .unwrap();
            for i in 0..4 {
                close(
                    s.derivatives[i],
                    if i == varying { 4e5 } else { 0.0 },
                    3.0 / 1e-6,
                );
            }
        }
    }

    const ROUNDED: [f64; 10] = [
        -3.7e-6, -3.4e-6, -3.1e-6, -2.8e-6, -2.5e-6, -2.2e-6, -1.9e-6, -1.6e-6, -1.3e-6, -1.0e-6,
    ];

    #[test]
    fn knot_sidedness_and_neighbours_on_every_axis() {
        // Nodal values i^2 give different slopes in every cell. Rounded canonical
        // knots exercise division rounding that a globally affine test hides.
        let a = UniformAxis::try_from_coordinates(&ROUNDED).unwrap();
        for varying in 0..4 {
            let g = grid([&ROUNDED; 4], |q| {
                let i = (0..a.len())
                    .find(|&i| a.coordinate(i).unwrap() == q[varying])
                    .unwrap();
                (i * i) as f64
            });
            for i in 1..a.len() - 1 {
                let knot = a.coordinate(i).unwrap();
                for (q, expected_cell) in [(next_down(knot), i - 1), (knot, i), (next_up(knot), i)]
                {
                    let (cell, u) = locate(&a, q, varying).unwrap();
                    assert_eq!(cell, expected_cell);
                    if q == knot {
                        assert_eq!(u, 0.0);
                    }
                    let mut point = [-2.91e-6; 4];
                    point[varying] = q;
                    let s = g.sample_stored(point).unwrap();
                    close(
                        s.derivatives[varying],
                        (2 * expected_cell + 1) as f64 / a.spacing(),
                        20.0 / a.spacing(),
                    );
                    if q == knot {
                        close(s.a_sqd, (i * i) as f64, 100.0);
                    }
                }
            }
            let mut point = [-2.91e-6; 4];
            point[varying] = a.coordinate(a.len() - 1).unwrap();
            let s = g.sample_stored(point).unwrap();
            close(s.a_sqd, 81.0, 100.0);
            close(
                s.derivatives[varying],
                17.0 / a.spacing(),
                20.0 / a.spacing(),
            );
            assert_eq!(locate(&a, point[varying], varying).unwrap(), (8, 1.0));
        }
    }

    #[test]
    fn canonical_faces_corners_and_adjacent_queries() {
        let g = grid([&ROUNDED; 4], |_| 2.0);
        let min = g.axes()[0].coordinate(0).unwrap();
        let max = g.axes()[0].coordinate(9).unwrap();
        for bits in 0..16 {
            let q = std::array::from_fn(|a| if bits & (1 << a) == 0 { min } else { max });
            close(g.sample_stored(q).unwrap().a_sqd, 2.0, 2.0);
        }
        for axis in 0..4 {
            for value in [min, next_up(min), next_down(max), max] {
                let mut q = [-2.91e-6; 4];
                q[axis] = value;
                close(g.sample_stored(q).unwrap().a_sqd, 2.0, 2.0);
            }
            for coordinate in [next_down(min), next_up(max)] {
                let mut q = [-2.91e-6; 4];
                q[axis] = coordinate;
                assert_eq!(
                    g.sample_stored(q),
                    Err(SampleError::OutOfDomain {
                        axis,
                        coordinate,
                        min,
                        max
                    })
                );
            }
        }
    }

    #[test]
    fn supplied_endpoints_follow_strict_canonical_bounds() {
        // Original regression: canonical upper endpoint moves outward.
        // This separate decimal fixture has an inward canonical upper endpoint.
        let inward = [-1.3e-6, -1.0e-6, -0.7e-6, -0.4e-6, -0.1e-6];
        let mut movements = [false; 2];
        for nodes in [&ROUNDED[..], &inward[..]] {
            let g = grid([nodes; 4], |_| 1.0);
            let a = &g.axes()[0];
            let min = a.coordinate(0).unwrap();
            let max = a.coordinate(a.len() - 1).unwrap();
            let supplied = nodes[nodes.len() - 1];
            if supplied < max {
                movements[0] = true;
            }
            if supplied > max {
                movements[1] = true;
            }
            for axis in 0..4 {
                let mut point = [min; 4];
                point[axis] = supplied;
                if supplied > max {
                    assert_eq!(
                        g.sample_stored(point),
                        Err(SampleError::OutOfDomain {
                            axis,
                            coordinate: supplied,
                            min,
                            max
                        })
                    );
                } else {
                    close(g.sample_stored(point).unwrap().a_sqd, 1.0, 1.0);
                }
                point[axis] = max;
                assert!(g.sample_stored(point).is_ok());
            }
        }
        assert_eq!(movements, [true, true]);
    }

    #[test]
    fn nonfinite_queries_identify_each_axis() {
        let a = [1.0, 2.0];
        let g = grid([&a; 4], |_| 1.0);
        for axis in 0..4 {
            for coordinate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut point = [1.5; 4];
                point[axis] = coordinate;
                let error = g.sample_stored(point).unwrap_err();
                assert!(matches!(error, SampleError::NonFiniteCoordinate {axis:a,..} if a==axis));
                assert!(error.to_string().contains(&format!(
                    "nonfinite {} coordinate",
                    ["x", "y", "z", "xi"][axis]
                )));
            }
            let mut point = [1.5; 4];
            point[axis] = 3.0;
            let message = g.sample_stored(point).unwrap_err().to_string();
            assert!(message.contains(&format!("{} coordinate", ["x", "y", "z", "xi"][axis])));
            assert!(message.contains("outside canonical domain [1, 2] m"));
        }
    }
}
