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

Run results are recorded in [progress](numerical-envelope-progress.md). Tests of
node indices do not verify continuous out-of-domain sampling, which is unimplemented.

## Planned gates (not yet verified)

Interpolation: constant S; affine derivatives `(2,3,5,7)` yielding raised lab
gradient `(7,-2,-3,2)`; positive multi-affine mixed terms; exact cell endpoints
and internal knots. Test xi-only and stored-z-only profiles separately.
Use laboratory finite differences with ct held fixed while changing laboratory z.

The future sampling domain uses canonical endpoints. Test both supplied and
canonical endpoints, including fixtures where reconstruction moves the upper
bound inward and outward, and immediately adjacent representable coordinates.
Canonical bounds are inclusive; supplied endpoints are classified by their actual
position relative to those bounds (boundary, interior, or out of domain). Do not
reuse axis-validation tolerances to expand the domain or silently clamp queries.

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
sampling or convergence measurements exist in step 1.
