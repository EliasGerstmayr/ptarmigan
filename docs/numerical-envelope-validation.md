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

## Pending benchmarks (not run)

Simulation runs will be performed separately by the user. None of these items
is marked complete by the exact-function or finite-difference tests:

- [ ] Analytical Gaussian values and all gradient components.
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

Gaussian: reproducible off-grid points, both polarizations, symmetry planes,
translated profiles and nonzero origins. Independently finite-difference the
analytic scalar. At symmetry knots the multilinear one-sided derivative need
not vanish; also test symmetry planes at cell centres. Compare native values and
spatial gradients. The native time-gradient finding is a **source-level discrepancy
awaiting numerical verification**. The source suggests
`g_native[0] - g_reference[0] = -S_z`; numerical confirmation and any trajectory
consequences remain unestablished.

Future diagnostic: choose nonzero amplitude, `x=y=0`, `z=z_R != 0`, and the
temporal envelope peak `xi=0` (thus `ct=z`). Here `S_xi=0`, while the independent
stored-z derivative is generally nonzero. Evaluate the native scalar directly:

```text
[FocusedLaser::a_sqd((ct+h,x,y,z)) - FocusedLaser::a_sqd((ct-h,x,y,z))] / (2h)
```

Hold laboratory x, y and z fixed; h is in metres. Compare this independent
estimate of `dS/d(ct)` with `FocusedLaser::grad_a_sqd(r)[0]` and the analytic
`B F'`. Vary h to check roundoff/truncation, and add an off-peak control point
where the time derivative is nonzero. Report scale-normalized absolute errors
near the peak's zero derivative. This diagnostic is planned, not implemented or
run by the storage tests; native field physics remains unchanged.

Files: missing/wrong metadata, unknown version, wrong rank/type/shape/units,
invalid axes and scalar values, resource limits and multi-rank rejection.

Convergence: begin with 16/32/64 cells per axis on a fixed domain
`x,y in [-2w0,2w0]`, `z in [-2z_R,2z_R]`, `xi in [-2cT,2cT]`. Reuse deterministic
off-grid points across resolutions. Report max and RMS absolute error, normalized
by S_peak, S_peak/w0, S_peak/z_R and S_peak/(cT) as appropriate; use
`S_peak/z_R + S_peak/(cT)` for the combined lab longitudinal component. Report
`log2(error_coarse/error_fine)`. Refine if provisional value `1e-3` and gradient
`1e-2` targets are missed; these are engineering targets, not physics tolerances.

Performance: ignored release tests for `gaussian_convergence` and `sampling_cost`
will use warm-up, precomputed coherent/random queries, black_box, repeated batches,
and median ns/sample. Measure generation/loading separately. Report allocation
and peak memory, releasing each resolution before generating the next. No field
sampling cost or Gaussian convergence measurements exist yet.
