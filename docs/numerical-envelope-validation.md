# Numerical envelope validation

## Step 1 checks

Seven unit tests in `src/field/envelope/grid.rs` cover:

- Shifted origins, unequal lengths/spacings, all 120 nodes of a `[2,3,4,5]`
  sentinel grid, xi-fastest indexing, out-of-range indices and retained allocation.
- Empty/singleton, repeated/decreasing, nonfinite/nonuniform and poorly resolved
  axes; overflow/underflow spacing; realistic metre-scale rounded coordinates.
- Shape mismatch, axis identification and negative/nonfinite scalar values.
- Checked element-count, byte-count and address-space limits without huge allocations.
- Invalid wavelength, both polarization metadata values, zero S, and no rescaling.
- Informative error messages.
- Canonical reconstruction of a decimal, micrometre-scale axis: coordinates and
  spacing stay close to independently specified physical values, remain increasing,
  and the reconstructed upper endpoint differs from the supplied endpoint.

Run results are recorded in [progress](numerical-envelope-progress.md). The node
tests alone do not verify interpolation; separate sampling tests follow.

## Stored-coordinate interpolation checks

Nine tests in `src/field/envelope/multilinear.rs` cover:

- Constant positive S at off-grid points and endpoints; exactly zero derivatives
  from paired identical nodal values.
- Independent affine references with distinct coefficients, shifted origins,
  unequal spacings and `[2,3,4,5]` axis lengths.
- A positive product of four affine factors, with derivatives formed analytically
  from products of the other three factors, exercising mixed-axis weights.
- Micrometre-scale affine data with coefficients of order `1e5` to `1e6 /m`.
- Separate xi-only and stored-z-only data, checking all derivative components.
- Piecewise nodal slopes differing at each knot on each axis: exact canonical
  knots select the positive side; adjacent floating-point queries retain the
  appropriate side; the upper endpoint selects the final cell at fraction one.
- All 16 endpoint corners, faces, and immediately inside/outside representable
  coordinates along every axis of a shifted, rounded grid.
- Supplied upper endpoints inside and outside reconstructed canonical bounds,
  with both movement directions explicitly checked by the fixtures.
- NaN and both infinities on every axis; error variants, axis identification,
  and informative out-of-domain messages.

Tolerance helper: `128 * f64::EPSILON * scale`. The fixtures have well-resolved
cells and moderate values; this allowance covers accumulated corner arithmetic
and cancellation. Scale is a bound on S for values and that bound divided by
spacing for derivatives, reflecting subtraction conditioning. Constant
derivatives and cell indices/fractions at exact knots are tested exactly.
No relative division by a zero derivative is used. These tests verify exact
functions and boundary semantics; they do not measure physics accuracy,
convergence order or performance. See progress for actual commands/results.

## Laboratory-coordinate checks

Ten focused tests in `laboratory.rs` cover:

- Constant S with zero gradient; shifted unequal axes with stored affine
  derivatives `(2,3,5,7)` and raised lab gradient `(7,-2,-3,2)`.
- Xi-only values, including equal simultaneous ct/z shifts; stored-z-only values
  and gradients, including time independence at fixed laboratory z.
- Known micrometre-scale laboratory and independently specified stored points.
- Independent central differences of the sampled scalar in each laboratory
  coordinate, at two interior points of a positive multi-affine fixture. Perturbing
  laboratory z holds ct fixed. Three steps `h=1e-2,1e-3,1e-4 m` stay in the same
  cells of metre-scale mathematical grids. These are small mathematical fixtures,
  not physical simulations or convergence studies.
- All 16 transformed canonical endpoint corners, with binary-exact boundary
  arithmetic checked independently; outside queries on every stored axis and a
  representable point immediately beyond the transformed xi upper bound.
- NaN and both infinities in every laboratory input; both signs of ct-z overflow;
  error variants and diagnostic context.
- End-to-end derivative overflow from finite nodal inputs, and longitudinal
  conversion overflow despite finite stored derivatives.
- Direct result-guard tests using synthetic nonfinite scalars and derivatives,
  clearly separate from the end-to-end overflow tests; a large finite scalar is
  accepted without imposing an arbitrary physical limit.

Direct analytical checks use `128*EPSILON*scale`, with S or S/cell-width scales
to allow accumulated interpolation/subtraction roundoff. The finite-difference
fixture is linear along laboratory ct/x/y and quadratic along laboratory z;
central differences have zero truncation error in exact arithmetic. Its error
bound is `128*EPSILON*40/h`, where 40 bounds the scalar magnitude in the fixture,
to allow roundoff amplification as h shrinks. These checks verify the chain rule
and sign/unit conventions without calling the derivative transformation in the
reference calculation. They establish no Gaussian or trajectory accuracy.

## Small Gaussian checks actually run

Nine additional ordinary tests cover reference peak normalization (both
polarizations), intensity FWHM, transverse exp(-2) radius, parameter/derived-scale
validation, independent scalar finite differences in all stored and laboratory
directions, symmetry and translations, coarse sampler agreement, cell-centred
sampler symmetry, checked allocations, canonical-node values, error normalization,
rate edge cases and deterministic reusable queries. Related assertions are grouped
in tests; there are nine test functions, not one function per assertion.

Reference finite differences use steps 1e-3, 1e-4, 1e-5 times characteristic
lengths, at generic translated off-axis/off-focus points. The normalized tolerance
`20*epsilon^2 + 128*machine_epsilon/epsilon` allows second-order truncation and
roundoff amplification. At fixed ct, perturbing laboratory z changes both stored
z and xi. No production derivative helper is used in the reference.

The coarse sampler fixture uses eight cells per axis over +/-0.5 characteristic
lengths (9^4 nodes, 52,488 scalar bytes). Cell widths are 0.125 in characteristic
units. Values use a 0.05 peak-normalized bound: the tensor linear interpolation
remainder scales as sum(h_i^2 * max|S_ii|)/8. Gradients use a 0.5 characteristic-
normalized bound allowing the leading half-cell times curvature error plus
transverse interpolation error. These central-box coarse bounds are deliberately
distinct from the provisional 1e-3/1e-2 engineering targets. Three cells per axis
place symmetry planes between nodes for cancellation tests; no zero-derivative
claim is made at a one-sided knot.

The ordinary tests and a tiny invocation/serialization smoke run passed; counts
and actual parameters are recorded in progress. This is not evidence of asymptotic
Gaussian convergence or trajectory accuracy.

## Gaussian tooling commands

Run from the repository root with the existing Homebrew toolchain and lockfile:

```bash
export PATH="$(brew --prefix rustup)/bin:$PATH"
```

Small focused tests (both diagnostic runners stay ignored):

```bash
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
```

Tiny smoke example, solely for invocation/output serialization:

```bash
ENVELOPE_CELLS=2,3,4 ENVELOPE_SAMPLES=16 ENVELOPE_MAX_MIB=1 \
ENVELOPE_SEED=20260916 ENVELOPE_OUTPUT=output/envelope-smoke \
cargo test --locked --release -p ptarmigan --no-default-features \
field::envelope::validation::gaussian_convergence -- --ignored --exact --nocapture --test-threads=1
```

Full three-resolution Gaussian study, **prepared but not run**:

```bash
ENVELOPE_CELLS=16,32,64 ENVELOPE_SAMPLES=1024 ENVELOPE_MAX_MIB=256 \
ENVELOPE_SEED=20260916 ENVELOPE_OUTPUT=output/envelope-gaussian \
cargo test --locked --release -p ptarmigan --no-default-features \
field::envelope::validation::gaussian_convergence -- --ignored --exact --nocapture --test-threads=1
```

Separate native diagnostic, **prepared but not run** (leave offsets zero):

```bash
ENVELOPE_X_C_M=0 ENVELOPE_Y_C_M=0 ENVELOPE_Z_F_M=0 ENVELOPE_XI_C_M=0 \
ENVELOPE_OUTPUT=output/envelope-native \
cargo test --locked --release -p ptarmigan --no-default-features \
field::envelope::validation::native_gaussian_diagnostic -- --ignored --exact --nocapture --test-threads=1
```

Use fresh output directories; existing report files are not overwritten. These
commands use no HDF5 feature and require no HDF5_DIR. Do not run all ignored
tests together: invoke the exact diagnostic selected above.

### Configuration and regeneration example

All physical values are SI. Both polarizations are always evaluated. Defaults:

| Environment variable | Default / meaning |
| --- | --- |
| `ENVELOPE_WAVELENGTH_M` | `0.8e-6` |
| `ENVELOPE_WAIST_M` | `5e-6` |
| `ENVELOPE_DURATION_S` | `30e-15`, intensity FWHM |
| `ENVELOPE_A0` | `1` |
| `ENVELOPE_X_C_M`, `ENVELOPE_Y_C_M` | `0`, transverse centre |
| `ENVELOPE_Z_F_M` | `0`, focal z |
| `ENVELOPE_XI_C_M` | `0`, retarded-coordinate offset |
| `ENVELOPE_EXTENT` | `2`, domain +/-extent times w0,w0,zR,cT about offsets |
| `ENVELOPE_CELLS` | `16,32,64`, positive strictly increasing counts; nodes=cells+1 |
| `ENVELOPE_SAMPLES` | `1024`, fixed interior query count |
| `ENVELOPE_SEED` | `20260916`, u64 RNG seed |
| `ENVELOPE_MAX_MIB` | `256`, scalar-allocation ceiling only |
| `ENVELOPE_OUTPUT` | `output/envelope-gaussian` or `output/envelope-native` |

Example physical overrides, applied before the Gaussian command:

```bash
export ENVELOPE_WAVELENGTH_M=0.8e-6 ENVELOPE_WAIST_M=5e-6 ENVELOPE_DURATION_S=30e-15
export ENVELOPE_A0=1 ENVELOPE_X_C_M=1e-6 ENVELOPE_Y_C_M=-2e-6
export ENVELOPE_Z_F_M=3e-6 ENVELOPE_XI_C_M=-3e-6
```

That xi offset places the focus/temporal peak at ct=0. Unset overrides to restore
defaults; metadata records the actual settings. The native diagnostic rejects
nonzero offsets and does not use configured grid counts or random samples.

Gaussian outputs: `metadata.txt`, `axes.csv`, `queries.csv`, `errors.csv`.
Native outputs: `metadata.txt`, `native.csv` (all four components and all three
finite-difference steps; scalar comparisons included). Estimated scalar memory is
8*(cells+1)^4 bytes, printed before each allocation and recorded per resolution.
It is not measured peak process memory. Each resolution's grid is released before
the next; queries and small axis/report structures persist.

## Pending benchmarks (not run at production resolutions)

Simulation runs will be performed separately by the user. None of these items
is marked complete by the exact-function or finite-difference tests:

- [ ] Broad analytical Gaussian value/gradient error study (small regression fixtures passed).
- [ ] Three-resolution Gaussian convergence.
- [ ] Independent numerical check of the native Gaussian time-gradient discrepancy.
- [ ] Sampling cost and memory measurements.
- [ ] Radiation-free trajectories.
- [ ] LMA radiation.
- [ ] LASY Gaussian and matched-history flying-focus comparisons.
- [ ] Later LCFA and finite-beam validation.

Matched-history comparisons isolate local radiation calculations; they do not
reproduce spatial gradients or finite-beam dynamics.

## Remaining validation design (not yet verified)

Gaussian tooling now provides reproducible off-grid points, both polarizations,
translations and independent scalar checks. Larger measurements remain pending.
At symmetry knots the multilinear one-sided derivative need not vanish. The
native time-gradient finding is a **source-level discrepancy
awaiting numerical verification**. The source suggests
`g_native[0] - g_reference[0] = -S_z`; numerical confirmation and any trajectory
consequences remain unestablished.

Prepared native diagnostic: choose nonzero amplitude, `x=y=0`, `z=z_R/2 != 0`, and the
temporal envelope peak `xi=0` (thus `ct=z`). Here `S_xi=0`, while the independent
stored-z derivative is generally nonzero. Evaluate the native scalar directly:

```text
[FocusedLaser::a_sqd((ct+h,x,y,z)) - FocusedLaser::a_sqd((ct-h,x,y,z))] / (2h)
```

Hold laboratory x, y and z fixed; h is in metres. Compare this independent
estimate of `dS/d(ct)` with `FocusedLaser::grad_a_sqd(r)[0]` and the analytic
`B F'`. Vary h to check roundoff/truncation, and add an off-peak control point
where the time derivative is nonzero. Report scale-normalized absolute errors
near the peak's zero derivative. This diagnostic is implemented as an ignored
runner, but has not been run; native field physics remains unchanged.

Files: missing/wrong metadata, unknown version, wrong rank/type/shape/units,
invalid axes and scalar values, resource limits and multi-rank rejection.

Convergence: begin with 16/32/64 cells per axis on a fixed domain
`x,y in [-2w0,2w0]`, `z in [-2z_R,2z_R]`, `xi in [-2cT,2cT]`. Reuse deterministic
off-grid points across resolutions. Report max and RMS absolute error, normalized
by S_peak, S_peak/w0, S_peak/z_R and S_peak/(cT) as appropriate; use
`S_peak/z_R + S_peak/(cT)` for the combined lab longitudinal component. Report
`log2(error_coarse/error_fine)`. Refine if provisional value `1e-3` and gradient
`1e-2` targets are missed; these are engineering targets, not physics tolerances.

Performance: a future sampling-cost runner should use warm-up, precomputed
coherent/random queries, black_box, repeated batches and median ns/sample.
Measure generation/loading separately. No sampling-cost or peak-memory runner
has been implemented or run; the Gaussian runner reports estimated scalar storage
only. Full-resolution Gaussian convergence is prepared but unrun.
