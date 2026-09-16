use super::{
    gaussian::{Gaussian, Parameters},
    grid,
    report::Report,
    statistics::Errors,
    Result,
};
use crate::field::{envelope::UniformAxis, Polarization};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::{env, path::PathBuf, str::FromStr};

pub struct Config {
    pub params: Parameters,
    pub cells: Vec<usize>,
    pub samples: usize,
    pub seed: u64,
    pub limit: usize,
    pub extent: f64,
    pub output: PathBuf,
}

fn setting<T: FromStr>(key: &str, default: T) -> Result<T> {
    match env::var(key) {
        Ok(s) => s
            .parse()
            .map_err(|_| format!("invalid {}: {}", key, s).into()),
        Err(env::VarError::NotPresent) => Ok(default),
        Err(e) => Err(e.into()),
    }
}

impl Config {
    pub fn from_env(mode: &str) -> Result<Self> {
        let p = Parameters::default();
        let cells = env::var("ENVELOPE_CELLS")
            .unwrap_or_else(|_| "16,32,64".into())
            .split(',')
            .map(|s| s.trim().parse())
            .collect::<std::result::Result<Vec<usize>, _>>()?;
        if cells.is_empty() || cells[0] == 0 || cells.windows(2).any(|w| w[1] <= w[0]) {
            return Err("ENVELOPE_CELLS must be strictly increasing positive cell counts".into());
        }
        let params = Parameters {
            wavelength: setting("ENVELOPE_WAVELENGTH_M", p.wavelength)?,
            waist: setting("ENVELOPE_WAIST_M", p.waist)?,
            duration: setting("ENVELOPE_DURATION_S", p.duration)?,
            a0: setting("ENVELOPE_A0", p.a0)?,
            offsets: [
                setting("ENVELOPE_X_C_M", 0.0)?,
                setting("ENVELOPE_Y_C_M", 0.0)?,
                setting("ENVELOPE_Z_F_M", 0.0)?,
                setting("ENVELOPE_XI_C_M", 0.0)?,
            ],
        };
        let samples = setting("ENVELOPE_SAMPLES", 1024usize)?;
        let extent = setting("ENVELOPE_EXTENT", 2.0f64)?;
        if samples == 0 || !extent.is_finite() || extent <= 0.0 {
            return Err("samples and finite extent must be positive".into());
        }
        let limit = setting("ENVELOPE_MAX_MIB", 256usize)?
            .checked_mul(1024 * 1024)
            .ok_or("memory limit overflow")?;
        for pol in [Polarization::Linear, Polarization::Circular] {
            Gaussian::new(params, pol)?;
        }
        Ok(Self {
            params,
            cells,
            samples,
            extent,
            limit,
            seed: setting("ENVELOPE_SEED", 20260916)?,
            output: env::var_os("ENVELOPE_OUTPUT")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from(format!("output/envelope-{}", mode))),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Query {
    pub stored: [f64; 4],
    pub lab: [f64; 4],
}

pub fn queries(
    g: &Gaussian,
    cells: &[usize],
    extent: f64,
    count: usize,
    seed: u64,
) -> Result<Vec<Query>> {
    // Axes only: never allocate a scalar grid while choosing reusable query points.
    let mut all_axes = Vec::new();
    for &n in cells {
        let nodes = grid::axes(g, n, extent)?;
        let axes = nodes
            .iter()
            .map(|a| UniformAxis::try_from_coordinates(a))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        all_axes.push(axes);
    }
    let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
    let mut points = Vec::new();
    points.try_reserve_exact(count)?;
    let attempts = count
        .checked_mul(1000)
        .ok_or("sample attempt limit overflow")?;
    for _ in 0..attempts {
        let q: [f64; 4] = std::array::from_fn(|i| {
            g.params.offsets[i] + 0.95 * extent * g.lengths()[i] * (2.0 * rng.gen::<f64>() - 1.0)
        });
        let lab = [q[2] + q[3], q[0], q[1], q[2]];
        // Use the effective round-tripped xi for BOTH scalar-grid and reference
        // queries; do not count coordinate reconstruction as interpolation error.
        let stored = [lab[1], lab[2], lab[3], lab[0] - lab[3]];
        if lab.iter().chain(stored.iter()).any(|v| !v.is_finite()) {
            continue;
        }
        let valid = (0..4).all(|i| {
            let away_from_symmetry =
                ((stored[i] - g.params.offsets[i]) / g.lengths()[i]).abs() > 1e-6;
            away_from_symmetry
                && all_axes.iter().all(|axes| {
                    let a = &axes[i];
                    let min = a.coordinate(0).unwrap();
                    let max = a.coordinate(a.len() - 1).unwrap();
                    // Exact canonical comparisons plus a benchmark-only exclusion
                    // band around knots. The production domain is never broadened.
                    stored[i] > min
                        && stored[i] < max
                        && (0..a.len()).all(|j| {
                            (stored[i] - a.coordinate(j).unwrap()).abs() > 1e-7 * a.spacing()
                        })
                })
        });
        if valid {
            points.push(Query { stored, lab });
        }
        if points.len() == count {
            return Ok(points);
        }
    }
    Err("unable to generate requested interior points; check offsets/scales".into())
}

pub fn evaluate(
    g: &Gaussian,
    grid: &crate::field::envelope::EnvelopeGrid,
    points: &[Query],
) -> Result<[Errors; 9]> {
    let mut stats = [Errors::default(); 9];
    for q in points {
        let stored = grid.sample_stored(q.stored)?;
        let lab = grid.sample_lab(q.lab.into())?;
        let reference = g.stored(q.stored);
        let reference_lab = g.lab(q.lab);
        stats[0].add(stored.a_sqd, reference[0])?;
        for i in 0..4 {
            stats[i + 1].add(stored.derivatives[i], reference[i + 1])?;
            stats[i + 5].add(lab.grad_a_sqd[i as i32], reference_lab[i + 1])?;
        }
    }
    Ok(stats)
}

pub fn run() -> Result<()> {
    let c = Config::from_env("gaussian")?;
    for &n in &c.cells {
        grid::allocation(n, c.limit)?;
    }
    let geometry = Gaussian::new(c.params, Polarization::Linear)?;
    let points = queries(&geometry, &c.cells, c.extent, c.samples, c.seed)?;
    let mut report = Report::new(&c, "gaussian_convergence")?;
    report.queries(&points)?;
    let mut csv = report.file("errors.csv")?;
    report.error_header(&mut csv)?;
    let mut axes_csv = report.file("axes.csv")?;
    use std::io::Write;
    writeln!(
        axes_csv,
        "polarization,cells,axis,origin_m,spacing_m,nodes,canonical_upper_m,estimated_scalar_bytes"
    )?;
    for pol in [Polarization::Linear, Polarization::Circular] {
        let g = Gaussian::new(c.params, pol)?;
        report.gaussian(&g)?;
        let mut previous = None;
        for &cells in &c.cells {
            let (_, bytes) = grid::allocation(cells, c.limit)?;
            let stats = {
                let grid = grid::build(&g, cells, c.extent, c.limit)?;
                report.axes(&mut axes_csv, &g, &grid, cells, bytes)?;
                evaluate(&g, &grid, &points)?
            }; // grid dropped BEFORE constructing the next resolution
            report.errors(&mut csv, &g, cells, bytes, &stats, previous.as_ref())?;
            previous = Some((cells, stats));
        }
    }
    csv.flush()?;
    axes_csv.flush()?;
    report.finish()
}

#[test]
fn deterministic_queries_avoid_knots_and_roundtrip() {
    let g = Gaussian::new(Parameters::default(), Polarization::Linear).unwrap();
    let a = queries(&g, &[2, 3, 4], 2.0, 16, 7).unwrap();
    assert_eq!(a, queries(&g, &[2, 3, 4], 2.0, 16, 7).unwrap());
    assert_ne!(a, queries(&g, &[2, 3, 4], 2.0, 16, 8).unwrap());
    for q in a {
        assert_eq!(q.lab[0] - q.lab[3], q.stored[3]);
    }
}
