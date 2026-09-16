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

Tiny smoke example, solely for invocation/output serialization and loading the
previous 16-query smoke file (use an existing saved file at this path):

```bash
ENVELOPE_CELLS=2,3,4 ENVELOPE_SAMPLES=16 ENVELOPE_MAX_MIB=1 \
ENVELOPE_QUERIES=/private/tmp/ptarmigan-gaussian-smoke-20260916/queries.csv \
ENVELOPE_OUTPUT=output/envelope-cubic-smoke \
cargo test --locked --release -p ptarmigan --no-default-features \
field::envelope::validation::gaussian_convergence -- --ignored --exact --nocapture --test-threads=1
```

Full three-resolution method comparison, **prepared but not run**. This loads
the exact queries from the completed multilinear study into a new output directory.
Use its default physical parameters, zero offsets and extent=2; clear conflicting
ENVELOPE_* overrides first:

```bash
ENVELOPE_CELLS=16,32,64 ENVELOPE_MAX_MIB=256 \
ENVELOPE_QUERIES=output/envelope-gaussian/queries.csv \
ENVELOPE_OUTPUT=output/envelope-cubic-comparison \
cargo test --locked --release -p ptarmigan --no-default-features \
field::envelope::validation::gaussian_convergence -- --ignored --exact --nocapture --test-threads=1
```

Separate native diagnostic command (already run by the user; **not rerun for
cubic work**). A repeat needs a fresh directory and zero offsets:

```bash
ENVELOPE_X_C_M=0 ENVELOPE_Y_C_M=0 ENVELOPE_Z_F_M=0 ENVELOPE_XI_C_M=0 \
ENVELOPE_OUTPUT=output/envelope-native-repeat \
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
| `ENVELOPE_CELLS` | `16,32,64`, strictly increasing counts >=2 for comparison; nodes=cells+1 |
| `ENVELOPE_SAMPLES` | `1024`, generated interior query count; ignored when loading |
| `ENVELOPE_QUERIES` | unset; optional saved numeric CSV, all rows used unchanged |
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

Gaussian outputs: `metadata.txt`, `axes.csv`, `queries.csv`, `errors.csv`,
`samples.csv`. Errors and signed-sample summaries separate both methods; axes are
shared. `ENVELOPE_QUERIES` requires finite required columns and exact parsed
`xi == ct-z`, with signed zero equality. All rows must fit every canonical grid;
no tolerance, relocation, dropping or regeneration is allowed. Saved 17-decimal
scientific numbers round-trip exactly. Loaded points need not avoid knots.
Reusing a seed alone is not a substitute for loading the original physical points.
Native outputs: `metadata.txt`, `native.csv` (all four components and all three
finite-difference steps; scalar comparisons included). Estimated scalar memory is
8*(cells+1)^4 bytes, printed before each allocation and recorded per resolution.
It is not measured peak process memory. Each resolution's grid is released before
the next; queries and small axis/report structures persist.

## Existing user-run results and pending benchmarks

Read-only inspection of `output/envelope-gaussian` confirms the completed
16/32/64-cell multilinear study at commit `6caf275`, 1024 saved queries, both
polarizations and default physical parameters. The supplied interpretation is
approximately second-order scalar and first-order gradient RMS convergence.
At 64 cells, sampled maximum normalized errors are about 2.329e-3 for S and
1.024e-1 for S_xi: provisional 1e-3/1e-2 sampled targets were missed. This
is convergence evidence, not satisfactory production accuracy or a domain bound.
The existing output directory was preserved. No full comparison ran in this task.

The separately completed `output/envelope-native` diagnostic confirms the native
time-gradient discrepancy at the tested points. On-axis at z=z_R/2 and xi=0,
the linear native time component is about 3259.493 /m while all three fixed-position
scalar finite differences and the analytical time derivative are zero. The
circular result is twice that value. Native physics and trajectory effects are
outside this change; no native diagnostic was rerun.


Simulation runs will be performed separately by the user. Completed items below
refer to the existing user-run outputs, not the small regression tests:

- [x] Initial multilinear three-resolution Gaussian study (user run, targets missed).
- [x] Native Gaussian time-gradient diagnostic at prepared points (user run).
- [ ] Full multilinear/cubic comparison on the saved physical queries.
- [ ] Positivity policy and broader error validation before particle use.
- [ ] Sampling cost and memory measurements.
- [ ] Radiation-free trajectories.
- [ ] LMA radiation.
- [ ] LASY Gaussian and matched-history flying-focus comparisons.
- [ ] Later LCFA and finite-beam validation.

Matched-history comparisons isolate local radiation calculations; they do not
reproduce spatial gradients or finite-beam dynamics.

## Remaining validation design (not yet verified)

Gaussian tooling now provides reproducible off-grid points, both polarizations,
translations and independent scalar checks. The larger cubic comparison remains
pending. At symmetry knots the multilinear one-sided derivative need not vanish.
The native time-gradient source prediction `g_native[0] - g_reference[0] = -S_z`
is supported by the existing diagnostic at its tested points; trajectory
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
runner and was run separately by the user; native field physics remains unchanged.

Files: missing/wrong metadata, unknown version, wrong rank/type/shape/units,
invalid axes and scalar values, resource limits and multi-rank rejection.

Convergence: begin with 16/32/64 cells per axis on a fixed domain
`x,y in [-2w0,2w0]`, `z in [-2z_R,2z_R]`, `xi in [-2cT,2cT]`. Reuse deterministic
off-grid points across resolutions. Report max and RMS absolute error, normalized
by S_peak, S_peak/w0, S_peak/z_R and S_peak/(cT) as appropriate; use
`S_peak/z_R + S_peak/(cT)` for the combined lab longitudinal component. Report
`ln(error_coarse/error_fine)/ln(cells_fine/cells_coarse)`, separately for each
method and polarization. Provisional sampled value `1e-3` and gradient `1e-2`
targets are engineering targets, not physics tolerances. Report misses; decide
refinement separately, with no automatic refinement in this runner.

Performance: a future sampling-cost runner should use warm-up, precomputed
coherent/random queries, black_box, repeated batches and median ns/sample.
Measure generation/loading separately. No sampling-cost or peak-memory runner
has been implemented or run; the Gaussian runner reports estimated scalar storage
only. Full-resolution cubic comparison is prepared but unrun. The nominal 16
versus up to 256 distinct nodal values is not a measured runtime ratio.


## Small cubic validation

Six new mathematical tests cover constants (including a large common offset),
shifted affine fields on unequal axes, positive tensor-product quadratics,
interior/boundary cells, all nodal combinations, three-node operation, independent
left/right polynomials on every shared face, and micrometre finite differences
of a nonquadratic scalar. Laboratory z differences hold ct fixed. Other checks
cover immediate outside floats, NaN/infinity, arithmetic overflow, two-node
rejection only for cubic, signed overshoot and unchanged multilinear defaults.
Continuity tests use rounded canonical separations instead of assuming exact h.
Algebraic comparisons allow 512 machine epsilons times scalar or S/h scales for
stencil summation and coordinate roundoff. Finite differences use characteristic
micrometre scales with O(step^2) truncation and O(machine_epsilon/step) roundoff.

One new loader test checks saved-query round trips, invalid/missing/nonfinite/
inconsistent data, and a rounded-domain point accepted at one resolution but
rejected at another. No points are repaired. Existing multilinear tests remain.
The small fixture verifies S=-0.125 and an off-centre negative sample with its
unchanged derivative; these demonstrate why signed diagnostics cannot become a
particle-field policy without separate work. Counts and the tiny smoke's actual
signed samples are recorded in progress. Cubic orders near three/two remain
expected behaviour, not verified asymptotic rates.
