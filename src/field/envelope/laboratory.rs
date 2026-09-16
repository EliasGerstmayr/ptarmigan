//! Laboratory sampling; see [design](../../../docs/numerical-envelope.md#implemented-laboratory-sampling).

use super::{EnvelopeGrid, InterpolationMethod, SampleError, StoredEnvelopeSample};
use crate::geometry::FourVector;

/// Dimensionless S and its raised laboratory four-gradient, distinct from a stored sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LabEnvelopeSample {
    pub a_sqd: f64,
    /// (d/d(ct), -d/dx, -d/dy, -d/dz) S, all in inverse metres.
    pub grad_a_sqd: FourVector,
}

pub(super) fn finite(value: f64, quantity: &'static str) -> Result<(), SampleError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(SampleError::NonFiniteResult { quantity, value })
    }
}

fn checked_lab_sample(stored: StoredEnvelopeSample) -> Result<LabEnvelopeSample, SampleError> {
    finite(stored.a_sqd, "a_sqd")?;
    for (value, name) in
        stored
            .derivatives
            .iter()
            .zip(["stored S_x", "stored S_y", "stored S_z", "stored S_xi"])
    {
        finite(*value, name)?;
    }
    let [sx, sy, sz, sxi] = stored.derivatives;
    let gradient = [sxi, -sx, -sy, -sz + sxi];
    for (value, name) in gradient.iter().zip([
        "laboratory grad[ct]",
        "laboratory grad[x]",
        "laboratory grad[y]",
        "laboratory grad[z]",
    ]) {
        finite(*value, name)?;
    }
    Ok(LabEnvelopeSample {
        a_sqd: stored.a_sqd,
        grad_a_sqd: gradient.into(),
    })
}

impl EnvelopeGrid {
    /// Sample at r=(ct,x,y,z) in metres, reusing the stored interpolant.
    ///
    /// Since xi=ct-z, dS/d(ct)=S_xi and dS/dz at fixed ct is S_z-S_xi.
    /// Raising the derivative index with (+---) gives
    /// (S_xi,-S_x,-S_y,-S_z+S_xi), in inverse metres. S_z here holds stored xi
    /// fixed. The first component differentiates ct, not t: no extra c is needed.
    /// See the [design](../../../docs/numerical-envelope.md#implemented-laboratory-sampling).
    ///
    /// Nonfinite inputs, transformation overflow, transformed queries outside
    /// strict canonical bounds, and nonfinite arithmetic results are errors.
    /// No coordinate clamp, zero padding, or particle-termination policy is added.
    pub fn sample_lab(&self, r: FourVector) -> Result<LabEnvelopeSample, SampleError> {
        self.sample_lab_with_method(r, InterpolationMethod::Multilinear)
    }

    /// Explicit interpolant, with the same coordinate conversion and finite checks.
    /// Finite negative cubic samples are returned unchanged.
    pub fn sample_lab_with_method(
        &self,
        r: FourVector,
        method: InterpolationMethod,
    ) -> Result<LabEnvelopeSample, SampleError> {
        let position = [r[0], r[1], r[2], r[3]];
        for (component, &coordinate) in position.iter().enumerate() {
            if !coordinate.is_finite() {
                return Err(SampleError::NonFiniteLabCoordinate {
                    component,
                    coordinate,
                });
            }
        }
        let xi = r[0] - r[3];
        if !xi.is_finite() {
            return Err(SampleError::NonFiniteXi { ct: r[0], z: r[3] });
        }
        let stored = self
            .sample_stored_with_method([r[1], r[2], r[3], xi], method)
            .map_err(|error| match error {
                SampleError::OutOfDomain {
                    axis,
                    coordinate,
                    min,
                    max,
                } => SampleError::TransformedOutOfDomain {
                    position,
                    axis,
                    coordinate,
                    min,
                    max,
                },
                other => other,
            })?;
        checked_lab_sample(stored)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{
        envelope::{EnvelopeMetadata, UniformAxis},
        Polarization,
    };

    fn grid(nodes: [&[f64]; 4], f: impl Fn([f64; 4]) -> f64) -> EnvelopeGrid {
        let axes = nodes.map(|n| UniformAxis::try_from_coordinates(n).unwrap());
        let mut values = Vec::new();
        for x in 0..axes[0].len() {
            for y in 0..axes[1].len() {
                for z in 0..axes[2].len() {
                    for xi in 0..axes[3].len() {
                        let i = [x, y, z, xi];
                        values.push(f(std::array::from_fn(|a| {
                            axes[a].coordinate(i[a]).unwrap()
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
        assert!(
            (actual - expected).abs() <= 128.0 * f64::EPSILON * scale,
            "actual={:e}, expected={:e}, scale={:e}",
            actual,
            expected,
            scale
        );
    }

    fn check_gradient(actual: FourVector, expected: [f64; 4], scale: f64) {
        for i in 0..4 {
            close(actual[i as i32], expected[i], scale);
        }
    }

    #[test]
    fn constant_lab_sample() {
        let a = [-1.0, 0.0, 1.0];
        let g = grid([&a; 4], |_| 3.25);
        let s = g.sample_lab([0.7, 0.2, -0.3, 0.4].into()).unwrap();
        close(s.a_sqd, 3.25, 3.25);
        assert_eq!(s.grad_a_sqd, [0.0; 4].into());
    }

    #[test]
    fn affine_lab_gradient_shifted_unequal_axes() {
        let nodes: [&[f64]; 4] = [
            &[-2.0, 0.0],
            &[1.0, 1.5, 2.0],
            &[-3.0, -2.0, -1.0, 0.0],
            &[2.0, 2.25, 2.5, 2.75, 3.0],
        ];
        let g = grid(nodes, |q| {
            30.0 + 2.0 * q[0] + 3.0 * q[1] + 5.0 * q[2] + 7.0 * q[3]
        });
        // Independent mapping: ct=1, z=-1.5 gives xi=2.5.
        let s = g.sample_lab([1.0, -1.0, 1.25, -1.5].into()).unwrap();
        close(s.a_sqd, 41.75, 60.0);
        check_gradient(s.grad_a_sqd, [7.0, -2.0, -3.0, 2.0], 60.0 / 0.25);
    }

    #[test]
    fn xi_only_values_and_simultaneous_shift() {
        let a = [-1.0, 0.0, 1.0, 2.0];
        let g = grid([&a; 4], |q| 5.0 + 3.0 * q[3]);
        for r in [[0.75, 0.2, 0.3, 0.25], [1.25, 0.2, 0.3, 0.75]] {
            let s = g.sample_lab(r.into()).unwrap();
            close(s.a_sqd, 6.5, 10.0); // xi=0.5 for both independent positions
            check_gradient(s.grad_a_sqd, [3.0, 0.0, 0.0, 3.0], 10.0);
        }
        close(
            g.sample_lab([1.0, 0.2, 0.3, 0.25].into()).unwrap().a_sqd,
            7.25,
            10.0,
        );
    }

    #[test]
    fn stored_z_evolution_survives_lab_conversion() {
        let a = [-1.0, 0.0, 1.0];
        let g = grid([&a; 4], |q| 5.0 + 2.0 * q[2]);
        for ct in [0.0, 0.5] {
            let s = g.sample_lab([ct, 0.2, 0.3, 0.25].into()).unwrap();
            close(s.a_sqd, 5.5, 7.0);
            check_gradient(s.grad_a_sqd, [0.0, 0.0, 0.0, -2.0], 7.0);
        }
        close(
            g.sample_lab([0.5, 0.2, 0.3, 0.5].into()).unwrap().a_sqd,
            6.0,
            7.0,
        );
    }

    #[test]
    fn mixed_function_lab_finite_differences() {
        let a = [-1.0, 0.0, 1.0];
        let g = grid([&a; 4], |q| {
            (2.0 + 0.2 * q[0]) * (2.0 + 0.3 * q[1]) * (2.0 + 0.4 * q[2]) * (2.0 + 0.5 * q[3])
        });
        for r in [[0.65, 0.2, 0.3, 0.25], [0.9, 0.4, 0.6, 0.3]] {
            let s = g.sample_lab(r.into()).unwrap();
            // Independent central differences perturb one laboratory input at
            // a time. In particular ct is held fixed when laboratory z changes.
            // Along z this multi-affine function becomes quadratic, so central
            // differences have zero truncation error in exact arithmetic.
            // All steps stay inside one cell. Smaller h amplifies roundoff ~S/h.
            for h in [1e-2, 1e-3, 1e-4] {
                for i in 0..4 {
                    let mut plus = r;
                    plus[i] += h;
                    let mut minus = r;
                    minus[i] -= h;
                    let ds = (g.sample_lab(plus.into()).unwrap().a_sqd
                        - g.sample_lab(minus.into()).unwrap().a_sqd)
                        / (2.0 * h);
                    let raised = if i == 0 { ds } else { -ds };
                    let bound = 128.0 * f64::EPSILON * 40.0 / h;
                    assert!((s.grad_a_sqd[i as i32] - raised).abs() <= bound);
                }
            }
        }
    }

    #[test]
    fn micrometre_mapping_to_known_stored_point() {
        let x = [-2e-6, -1e-6, 0.0];
        let y = [1e-6, 2e-6, 3e-6];
        let z = [3e-6, 4e-6, 5e-6];
        let xi = [-3e-6, -2e-6, -1e-6];
        let g = grid([&x, &y, &z, &xi], |q| {
            10.0 + 2e5 * q[0] + 3e5 * q[1] + 5e5 * q[2] + 7e5 * q[3]
        });
        let s = g
            .sample_lab([2e-6, -1.5e-6, 1.5e-6, 4.5e-6].into())
            .unwrap();
        let stored = g.sample_stored([-1.5e-6, 1.5e-6, 4.5e-6, -2.5e-6]).unwrap();
        close(s.a_sqd, 10.65, 15.0);
        close(s.a_sqd, stored.a_sqd, 15.0);
        check_gradient(s.grad_a_sqd, [7e5, -2e5, -3e5, 2e5], 15.0 / 1e-6);
    }

    #[test]
    fn transformed_boundaries_and_outside_context() {
        // Binary-exact nodes and ct=z+xi ensure these intended corners really map
        // to canonical bounds, independently asserted before invoking sample_lab.
        let a = [0.5, 1.0, 1.5];
        let g = grid([&a; 4], |_| 2.0);
        for bits in 0..16 {
            let q: [f64; 4] = std::array::from_fn(|i| if bits & (1 << i) == 0 { 0.5 } else { 1.5 });
            let r = [q[2] + q[3], q[0], q[1], q[2]];
            assert_eq!(r[0] - r[3], q[3]);
            for i in 0..4 {
                assert_eq!(
                    g.axes()[i].coordinate(if q[i] == 0.5 { 0 } else { 2 }),
                    Some(q[i])
                );
            }
            close(g.sample_lab(r.into()).unwrap().a_sqd, 2.0, 2.0);
        }
        for (r, axis, coordinate) in [
            ([3.0, 1.0, 1.0, 1.0], 3, 2.0),
            ([2.0, 2.0, 1.0, 1.0], 0, 2.0),
            ([2.0, 1.0, 0.0, 1.0], 1, 0.0),
            ([2.0, 1.0, 1.0, 2.0], 2, 2.0),
        ] {
            let error = g.sample_lab(r.into()).unwrap_err();
            assert_eq!(
                error,
                SampleError::TransformedOutOfDomain {
                    position: r,
                    axis,
                    coordinate,
                    min: 0.5,
                    max: 1.5
                }
            );
            assert!(error.to_string().contains("maps via xi=ct-z to stored"));
        }
        // Immediately beyond the transformed xi upper endpoint is still outside.
        let ct = f64::from_bits(2.0f64.to_bits() + 1);
        assert!(ct - 0.5 > 1.5);
        assert!(matches!(
            g.sample_lab([ct, 1.0, 1.0, 0.5].into()),
            Err(SampleError::TransformedOutOfDomain { axis: 3, .. })
        ));
    }

    #[test]
    fn invalid_lab_inputs_and_xi_overflow() {
        let a = [0.0, 1.0];
        let g = grid([&a; 4], |_| 1.0);
        for component in 0..4 {
            for coordinate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut r = [0.5; 4];
                r[component] = coordinate;
                let error = g.sample_lab(r.into()).unwrap_err();
                assert!(
                    matches!(error,SampleError::NonFiniteLabCoordinate {component:c,..} if c==component)
                );
                assert!(error.to_string().contains(&format!(
                    "laboratory {} coordinate",
                    ["ct", "x", "y", "z"][component]
                )));
            }
        }
        for (ct, z) in [(f64::MAX, -f64::MAX), (-f64::MAX, f64::MAX)] {
            let error = g.sample_lab([ct, 0.5, 0.5, z].into()).unwrap_err();
            assert_eq!(error, SampleError::NonFiniteXi { ct, z });
            assert!(error.to_string().contains("transformed xi = ct - z"));
        }
    }

    #[test]
    fn derivative_and_longitudinal_overflow() {
        let a = [0.0, 0.5];
        let g = grid([&a; 4], |q| if q[0] == 0.0 { 0.0 } else { 1e308 });
        let error = g.sample_lab([0.5, 0.25, 0.25, 0.25].into()).unwrap_err();
        assert!(matches!(
            error,
            SampleError::NonFiniteResult {
                quantity: "stored S_x",
                ..
            }
        ));
        assert!(error.to_string().contains("nonfinite stored S_x"));

        let g = grid([&a; 4], |q| 5e307 + 1e308 * q[3] - 1e308 * q[2]);
        let stored = g.sample_stored([0.25; 4]).unwrap();
        assert!(stored.derivatives.iter().all(|v| v.is_finite()));
        let error = g.sample_lab([0.5, 0.25, 0.25, 0.25].into()).unwrap_err();
        assert!(matches!(
            error,
            SampleError::NonFiniteResult {
                quantity: "laboratory grad[z]",
                ..
            }
        ));
    }

    #[test]
    fn rejects_nonfinite_scalar_and_each_stored_derivative() {
        // Synthetic samples exercise the result guard directly, including NaN;
        // end-to-end finite-input overflow is covered separately above.
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                checked_lab_sample(StoredEnvelopeSample {
                    a_sqd: value,
                    derivatives: [0.0; 4]
                }),
                Err(SampleError::NonFiniteResult {
                    quantity: "a_sqd",
                    ..
                })
            ));
            for i in 0..4 {
                let mut derivatives = [0.0; 4];
                derivatives[i] = value;
                assert!(matches!(
                    checked_lab_sample(StoredEnvelopeSample {
                        a_sqd: 1.0,
                        derivatives
                    }),
                    Err(SampleError::NonFiniteResult { .. })
                ));
            }
        }
        // Large finite values are not rejected merely for exceeding a threshold.
        assert!(checked_lab_sample(StoredEnvelopeSample {
            a_sqd: f64::MAX,
            derivatives: [0.0; 4]
        })
        .is_ok());
    }
}
