use super::{gaussian::Gaussian, Result};
use crate::field::envelope::{EnvelopeGrid, EnvelopeMetadata, UniformAxis};

pub fn allocation(cells: usize, limit: usize) -> Result<(usize, usize)> {
    if cells == 0 {
        return Err("cells must be positive".into());
    }
    let nodes = cells.checked_add(1).ok_or("node count overflow")?;
    let count = nodes.checked_pow(4).ok_or("scalar count overflow")?;
    let bytes = count.checked_mul(8).ok_or("scalar byte count overflow")?;
    if bytes > limit || bytes > isize::MAX as usize {
        return Err(format!(
            "estimated scalar allocation {} bytes exceeds limit {}",
            bytes, limit
        )
        .into());
    }
    Ok((count, bytes))
}

pub fn axes(g: &Gaussian, cells: usize, extent: f64) -> Result<[Vec<f64>; 4]> {
    if cells == 0 || !extent.is_finite() || extent <= 0.0 {
        return Err("invalid cells/domain extent".into());
    }
    let axes: [Vec<f64>; 4] = std::array::from_fn(|a| {
        (0..=cells)
            .map(|i| {
                g.params.offsets[a]
                    + extent * g.lengths()[a] * (2.0 * i as f64 / cells as f64 - 1.0)
            })
            .collect()
    });
    for axis in &axes {
        UniformAxis::try_from_coordinates(axis)?;
    }
    Ok(axes)
}

pub fn build(g: &Gaussian, cells: usize, extent: f64, limit: usize) -> Result<EnvelopeGrid> {
    let (count, bytes) = allocation(cells, limit)?;
    eprintln!(
        "cells={} nodes={} estimated scalar allocation={} bytes (not measured peak memory)",
        cells,
        cells + 1,
        bytes
    );
    let nodes = axes(g, cells, extent)?;
    let axes: Vec<_> = nodes
        .iter()
        .map(|n| UniformAxis::try_from_coordinates(n))
        .collect::<std::result::Result<_, _>>()?;
    let mut values = Vec::new();
    values.try_reserve_exact(count)?;
    for x in 0..=cells {
        for y in 0..=cells {
            for z in 0..=cells {
                for xi in 0..=cells {
                    let i = [x, y, z, xi];
                    // Values belong to reconstructed canonical nodes, not supplied decimals.
                    values
                        .push(g.scalar(std::array::from_fn(|a| axes[a].coordinate(i[a]).unwrap())));
                }
            }
        }
    }
    Ok(EnvelopeGrid::try_new(
        [&nodes[0], &nodes[1], &nodes[2], &nodes[3]],
        values,
        EnvelopeMetadata::try_new(g.params.wavelength, g.pol)?,
    )?)
}

#[cfg(test)]
mod tests {
    use super::super::gaussian::Parameters;
    use super::*;
    use crate::field::Polarization;

    #[test]
    fn small_grid_gaussian_comparison() {
        let g = Gaussian::new(
            Parameters {
                offsets: [1e-6, -2e-6, 3e-6, -1e-6],
                ..Parameters::default()
            },
            Polarization::Linear,
        )
        .unwrap();
        // Eight cells over +/-0.5 characteristic lengths: h/scale=0.125.
        let grid = build(&g, 8, 0.5, 1024 * 1024).unwrap();
        for factors in [[0.17, -0.21, 0.31, 0.09], [-0.11, 0.19, -0.23, -0.17]] {
            let q = std::array::from_fn(|i| g.params.offsets[i] + factors[i] * g.lengths()[i]);
            let exact = g.stored(q);
            let sample = grid.sample_stored(q).unwrap();
            // Hessian-based O(h^2) value bound (<0.05); O(h) gradients (<0.5)
            // on this central box. Coarse checks, not production accuracy gates.
            assert!((sample.a_sqd - exact[0]).abs() / g.peak < 0.05);
            for i in 0..4 {
                assert!((sample.derivatives[i] - exact[i + 1]).abs() / g.scales()[i + 1] < 0.5);
            }
            let r = [q[2] + q[3], q[0], q[1], q[2]];
            let lab = grid.sample_lab(r.into()).unwrap();
            let exact = g.lab(r);
            for i in 0..4 {
                assert!((lab.grad_a_sqd[i as i32] - exact[i + 1]).abs() / g.scales()[5 + i] < 0.5);
            }
        }
    }

    #[test]
    fn sampler_symmetry_planes_at_cell_centres() {
        let g = Gaussian::new(Parameters::default(), Polarization::Circular).unwrap();
        let grid = build(&g, 3, 0.5, 1024 * 1024).unwrap(); // odd cells: origin is not a knot
        for i in 0..4 {
            let a = &grid.axes()[i];
            assert!(a.coordinate(1).unwrap() < 0.0 && a.coordinate(2).unwrap() > 0.0);
        }
        let s = grid.sample_stored([0.0; 4]).unwrap();
        for i in 0..4 {
            assert!(s.derivatives[i].abs() / g.scales()[i + 1] < 128.0 * f64::EPSILON);
        }
    }

    #[test]
    fn allocation_limits_and_canonical_values() {
        assert!(allocation(usize::MAX, usize::MAX).is_err());
        assert!(allocation(64, 1024).is_err());
        assert!(allocation(0, 1024).is_err());
        assert_eq!(allocation(2, 1024).unwrap(), (81, 648));
        let g = Gaussian::new(Parameters::default(), Polarization::Linear).unwrap();
        let grid = build(&g, 3, 2.0, 4096).unwrap();
        let q = std::array::from_fn(|a| grid.axes()[a].coordinate(2).unwrap());
        assert_eq!(grid.node([2; 4]).unwrap(), g.scalar(q));
    }
}
