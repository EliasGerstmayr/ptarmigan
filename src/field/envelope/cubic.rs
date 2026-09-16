//! Tensor-product Hermite interpolation with shared second-order nodal slopes.
//! Finite signed overshoot is intentional: this is not yet a particle-field policy.
use super::{
    laboratory::finite, multilinear::locate, EnvelopeGrid, SampleError, StoredEnvelopeSample,
    UniformAxis,
};

#[derive(Clone, Copy, Default)]
struct Stencil {
    start: usize,
    len: usize,
    value: [f64; 4],
    derivative: [f64; 4],
}

fn stencil(axis: &UniformAxis, cell: usize, u: f64) -> Stencil {
    let n = axis.len();
    let boundary = cell == 0 || cell == n - 2;
    let start = if cell == 0 { 0 } else { cell - 1 };
    let mut s = Stencil {
        start,
        len: if boundary { 3 } else { 4 },
        ..Stencil::default()
    };
    let u2 = u * u;
    let u3 = u2 * u;
    let h = [
        2.0 * u3 - 3.0 * u2 + 1.0,
        u3 - 2.0 * u2 + u,
        -2.0 * u3 + 3.0 * u2,
        u3 - u2,
    ];
    let dh = [
        6.0 * u2 - 6.0 * u,
        3.0 * u2 - 4.0 * u + 1.0,
        -6.0 * u2 + 6.0 * u,
        3.0 * u2 - 2.0 * u,
    ];
    for (weights, basis) in [(&mut s.value, h), (&mut s.derivative, dh)] {
        weights[cell - start] += basis[0];
        weights[cell + 1 - start] += basis[2];
        for (node, factor) in [(cell, basis[1]), (cell + 1, basis[3])] {
            // Coefficients of h*m: second-order one-sided at endpoints,
            // centred at all shared interior nodes. Applied to differences below.
            let (indices, slope) = if node == 0 {
                ([0, 1, 2], [-1.5, 2.0, -0.5])
            } else if node == n - 1 {
                ([n - 3, n - 2, n - 1], [0.5, -2.0, 1.5])
            } else {
                ([node - 1, node, node + 1], [-0.5, 0.0, 0.5])
            };
            for k in 0..3 {
                weights[indices[k] - start] += factor * slope[k];
            }
        }
    }
    for d in &mut s.derivative {
        *d /= axis.spacing();
    }
    s
}

fn evaluate(grid: &EnvelopeGrid, s: [Stencil; 4]) -> Result<StoredEnvelopeSample, SampleError> {
    let anchor = grid.node(s.map(|a| a.start)).unwrap();
    let mut result = StoredEnvelopeSample {
        a_sqd: anchor,
        derivatives: [0.0; 4],
    };
    for x in 0..s[0].len {
        for y in 0..s[1].len {
            for z in 0..s[2].len {
                for xi in 0..s[3].len {
                    let j = [x, y, z, xi];
                    if j == [0; 4] {
                        continue;
                    } // anchor already read
                      // Subtract a common value before applying the linear operators.
                      // This preserves constant fields exactly and reduces cancellation.
                    let delta = grid
                        .node(std::array::from_fn(|a| s[a].start + j[a]))
                        .unwrap()
                        - anchor;
                    let w: [f64; 4] = std::array::from_fn(|a| s[a].value[j[a]]);
                    result.a_sqd += w.iter().product::<f64>() * delta;
                    for a in 0..4 {
                        let weight = (0..4)
                            .map(|b| if a == b { s[b].derivative[j[b]] } else { w[b] })
                            .product::<f64>();
                        result.derivatives[a] += weight * delta;
                    }
                }
            }
        }
    }
    finite(result.a_sqd, "a_sqd")?;
    for (d, name) in
        result
            .derivatives
            .iter()
            .zip(["stored S_x", "stored S_y", "stored S_z", "stored S_xi"])
    {
        finite(*d, name)?;
    }
    Ok(result)
}

pub(super) fn sample(
    grid: &EnvelopeGrid,
    q: [f64; 4],
) -> Result<StoredEnvelopeSample, SampleError> {
    let mut stencils = [Stencil::default(); 4];
    for a in 0..4 {
        let axis = &grid.axes()[a];
        if axis.len() < 3 {
            return Err(SampleError::UnsupportedCubicAxis {
                axis: a,
                nodes: axis.len(),
            });
        }
        // Reuse the unchanged canonical locator and knot sidedness.
        let (cell, u) = locate(axis, q[a], a)?;
        stencils[a] = stencil(axis, cell, u);
    }
    evaluate(grid, stencils)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{
        envelope::{EnvelopeMetadata, InterpolationMethod},
        Polarization,
    };
    const CUBIC: InterpolationMethod = InterpolationMethod::Cubic;

    fn grid(nodes: [&[f64]; 4], f: impl Fn([f64; 4]) -> f64) -> EnvelopeGrid {
        let axes = nodes.map(|a| UniformAxis::try_from_coordinates(a).unwrap());
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
    fn close(a: f64, b: f64, scale: f64) {
        assert!(
            (a - b).abs() <= 512.0 * f64::EPSILON * scale,
            "actual={:e}, expected={:e}, scale={:e}",
            a,
            b,
            scale
        );
    }
    fn quadratic(q: [f64; 4]) -> (f64, [f64; 4]) {
        let f: [f64; 4] =
            std::array::from_fn(|a| 2.0 + (a + 1) as f64 * 0.1 * q[a] + 0.2 * q[a] * q[a]);
        (
            f.iter().product(),
            std::array::from_fn(|a| {
                ((a + 1) as f64 * 0.1 + 0.4 * q[a])
                    * (0..4).filter(|&b| b != a).map(|b| f[b]).product::<f64>()
            }),
        )
    }

    #[test]
    fn constants_and_shifted_unequal_affine() {
        let axes: [&[f64]; 4] = [
            &[-2.0, 0.0, 2.0],
            &[1.0, 1.5, 2.0, 2.5],
            &[-3.0, -2.0, -1.0, 0.0, 1.0],
            &[2.0, 2.25, 2.5, 2.75, 3.0, 3.25],
        ];
        let constant = grid(axes, |_| 1e20);
        let q = [-1.3, 1.8, -0.7, 2.9];
        assert_eq!(
            sample(&constant, q).unwrap(),
            StoredEnvelopeSample {
                a_sqd: 1e20,
                derivatives: [0.0; 4]
            }
        );
        let coeff = [2.0, 3.0, 5.0, 7.0];
        let f = |q: [f64; 4]| 40.0 + (0..4).map(|a| coeff[a] * q[a]).sum::<f64>();
        let g = grid(axes, f);
        let s = g.sample_stored_with_method(q, CUBIC).unwrap();
        close(s.a_sqd, f(q), 100.0);
        for a in 0..4 {
            close(s.derivatives[a], coeff[a], 100.0 / g.axes()[a].spacing());
        }
        let lab = g
            .sample_lab_with_method([q[2] + q[3], q[0], q[1], q[2]].into(), CUBIC)
            .unwrap();
        close(lab.a_sqd, f(q), 100.0);
        for (a, v) in [7.0, -2.0, -3.0, 2.0].iter().enumerate() {
            close(lab.grad_a_sqd[a as i32], *v, 400.0);
        }
    }

    #[test]
    fn tensor_quadratic_interior_boundaries_and_nodes() {
        let a = [-1.0, -0.5, 0.0, 0.5, 1.0];
        let g = grid([&a; 4], |q| quadratic(q).0);
        for q in [
            [-0.83, 0.73, -0.29, 0.19],
            [0.3, -0.2, 0.1, -0.4],
            [-1.0, 1.0, 1.0, -1.0],
        ] {
            let (value, d) = quadratic(q);
            let s = sample(&g, q).unwrap();
            close(s.a_sqd, value, 40.0);
            for i in 0..4 {
                close(s.derivatives[i], d[i], 80.0);
            }
        }
        for x in 0..5 {
            for y in 0..5 {
                for z in 0..5 {
                    for xi in 0..5 {
                        let i = [x, y, z, xi];
                        let q = std::array::from_fn(|a| g.axes()[a].coordinate(i[a]).unwrap());
                        close(sample(&g, q).unwrap().a_sqd, g.node(i).unwrap(), 40.0);
                    }
                }
            }
        }
        let three = grid([&a[..3]; 4], |q| quadratic(q).0);
        close(
            sample(&three, [-0.7; 4]).unwrap().a_sqd,
            quadratic([-0.7; 4]).0,
            40.0,
        );
    }

    #[test]
    fn independently_selected_adjacent_polynomials_are_c1() {
        // Rounded micrometre nodes, non-polynomial values, all faces including
        // transitions from one-sided boundary stencils to centred interior ones.
        let a = [-3.7e-6, -3.4e-6, -3.1e-6, -2.8e-6, -2.5e-6];
        let g = grid([&a; 4], |q| {
            (0..4)
                .map(|i| 2.0 + (q[i] * 1e6 * (i + 1) as f64).sin())
                .product()
        });
        for direction in 0..4 {
            for knot in 1..4 {
                let mut left = std::array::from_fn(|i| stencil(&g.axes()[i], 1, 0.37));
                let mut right = left;
                let axis = &g.axes()[direction];
                // Actual canonical separation, rather than assuming it equals h.
                let u = (axis.coordinate(knot).unwrap() - axis.coordinate(knot - 1).unwrap())
                    / axis.spacing();
                left[direction] = stencil(axis, knot - 1, u);
                right[direction] = stencil(axis, knot, 0.0);
                let l = evaluate(&g, left).unwrap();
                let r = evaluate(&g, right).unwrap();
                close(l.a_sqd, r.a_sqd, 100.0);
                for i in 0..4 {
                    close(l.derivatives[i], r.derivatives[i], 100.0 / axis.spacing());
                }
            }
        }
    }

    #[test]
    fn scalar_finite_differences_in_stored_and_lab_micrometres() {
        let a = [-2e-6, -1e-6, 0.0, 1e-6, 2e-6];
        let g = grid([&a; 4], |q| {
            quadratic(q.map(|x| x / 1e-6)).0 + (q[0] / 1e-6 + 2.0 * q[2] / 1e-6).sin()
        });
        let q = [-1.3e-6, 0.3e-6, -0.4e-6, 1.2e-6];
        let r = [q[2] + q[3], q[0], q[1], q[2]];
        let s = sample(&g, q).unwrap();
        let lab = g.sample_lab_with_method(r.into(), CUBIC).unwrap();
        for epsilon in [1e-3, 1e-4, 1e-5] {
            let h = epsilon * 1e-6;
            for i in 0..4 {
                let mut plus = q;
                let mut minus = q;
                plus[i] += h;
                minus[i] -= h;
                let fd = (sample(&g, plus).unwrap().a_sqd - sample(&g, minus).unwrap().a_sqd)
                    / (2.0 * h);
                let tol =
                    (20.0 * epsilon * epsilon + 512.0 * f64::EPSILON / epsilon) * 100.0 / 1e-6;
                assert!((fd - s.derivatives[i]).abs() < tol);
                let mut plus = r;
                let mut minus = r;
                plus[i] += h;
                minus[i] -= h;
                // Perturb z alone: ct stays fixed and xi changes oppositely.
                let fd = (g.sample_lab_with_method(plus.into(), CUBIC).unwrap().a_sqd
                    - g.sample_lab_with_method(minus.into(), CUBIC).unwrap().a_sqd)
                    / (2.0 * h)
                    * if i == 0 { 1.0 } else { -1.0 };
                assert!((fd - lab.grad_a_sqd[i as i32]).abs() < tol);
            }
        }
    }

    #[test]
    fn strict_bounds_finiteness_and_small_grid_errors() {
        let a = [0.5, 1.0, 1.5];
        let g = grid([&a; 4], |_| 2.0);
        for bits in 0..16 {
            let q = std::array::from_fn(|i| if bits & (1 << i) == 0 { 0.5 } else { 1.5 });
            assert_eq!(sample(&g, q).unwrap().a_sqd, 2.0);
        }
        for axis in 0..4 {
            for coordinate in [
                f64::from_bits(0.5f64.to_bits() - 1),
                f64::from_bits(1.5f64.to_bits() + 1),
            ] {
                let mut q = [1.0; 4];
                q[axis] = coordinate;
                assert!(matches!(sample(&g,q),Err(SampleError::OutOfDomain{axis:a,..}) if a==axis));
            }
            for coordinate in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut q = [1.0; 4];
                q[axis] = coordinate;
                assert!(
                    matches!(sample(&g,q),Err(SampleError::NonFiniteCoordinate{axis:a,..}) if a==axis)
                );
                assert!(
                    matches!(g.sample_lab_with_method(q.into(),CUBIC),Err(SampleError::NonFiniteLabCoordinate{component:a,..}) if a==axis)
                );
            }
            let mut axes = [&a[..]; 4];
            axes[axis] = &a[..2];
            let short = grid(axes, |_| 2.0);
            assert_eq!(
                sample(&short, [0.75; 4]),
                Err(SampleError::UnsupportedCubicAxis { axis, nodes: 2 })
            );
            assert!(short.sample_stored([0.75; 4]).is_ok());
            assert!(short.sample_lab([1.5, 0.75, 0.75, 0.75].into()).is_ok());
        }
        let huge = grid([&[0.0, 0.25, 0.5]; 4], |q| {
            if q[0] == 0.0 {
                0.0
            } else {
                1e308
            }
        });
        assert!(matches!(
            sample(&huge, [0.1; 4]),
            Err(SampleError::NonFiniteResult { .. })
        ));
    }

    #[test]
    fn signed_overshoot_and_unchanged_multilinear_defaults() {
        let a = [0.0, 1.0, 2.0, 3.0];
        let g = grid([&a; 4], |q| if q[0] < 2.0 { 0.0 } else { 1.0 });
        let q = [0.25, 0.5, 0.5, 0.5];
        let r = [1.0, 0.25, 0.5, 0.5];
        let s = g.sample_stored_with_method(q, CUBIC).unwrap();
        // On this boundary cell the polynomial is (x^2-x)/2.
        close(s.a_sqd, -3.0 / 32.0, 1.0);
        close(s.derivatives[0], -0.25, 1.0);
        close(sample(&g, [0.5; 4]).unwrap().a_sqd, -0.125, 1.0);
        let h = 1e-4;
        let mut plus = q;
        let mut minus = q;
        plus[0] += h;
        minus[0] -= h;
        close(
            (sample(&g, plus).unwrap().a_sqd - sample(&g, minus).unwrap().a_sqd) / (2.0 * h),
            -0.25,
            1.0 / h,
        );
        let lab = g.sample_lab_with_method(r.into(), CUBIC).unwrap();
        close(lab.a_sqd, s.a_sqd, 1.0);
        close(lab.grad_a_sqd[1], 0.25, 1.0);
        assert_eq!(g.sample_stored(q).unwrap().a_sqd, 0.0);
        assert_eq!(
            g.sample_stored(q),
            g.sample_stored_with_method(q, InterpolationMethod::Multilinear)
        );
        assert_eq!(
            g.sample_lab(r.into()),
            g.sample_lab_with_method(r.into(), InterpolationMethod::Multilinear)
        );
        assert_eq!(g.sample_lab(r.into()).unwrap().a_sqd, 0.0);
    }
}
