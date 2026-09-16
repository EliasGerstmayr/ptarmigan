# Numerical envelope design

## Implemented: in-memory grid contract

`src/field/envelope/` holds an envelope-specific numerical core, independent of
HDF5 and of the `Field` trait. Existing field implementations are unchanged.
Validated construction, read-only metadata/storage access, node indexing and
stored-coordinate multilinear sampling and laboratory-coordinate sampling are
implemented. There is no simulation input or Field integration yet.

`UniformAxis::try_from_coordinates` accepts at least two finite, strictly
increasing, uniformly spaced nodes. Coordinates are metres. It stores the origin,
endpoint-derived spacing, and length, rather than retaining the coordinate array.
Canonical coordinates use `origin + index * spacing` with fused multiply-add.
Residuals must be within `1e-10 * spacing + 8 * EPSILON * max(abs(endpoints))`.
The roundoff term must not exceed `1e-3 * spacing`; this deliberately rejects
poorly resolved grids at very large origins. Spacing and its reciprocal must be
finite and positive. These are numerical validation rules, not physics tolerances.

The canonical grid is authoritative: supplied nodes are discarded after validation.
Even the reconstructed final coordinate can differ slightly from the supplied
endpoint through floating-point rounding. Interpolation bounds are
`coordinate(0)` and `coordinate(len - 1)`, inclusive, on each axis. Supplied scalar
samples are associated with these canonical nodes; preserving the original node
coordinates exactly is not part of this representation's contract.

`EnvelopeGrid::try_new` accepts four coordinate slices, an owned `Vec<f64>`, and
validated `EnvelopeMetadata`. Axis order is always `(x,y,z,xi)`, with xi fastest:

```text
index = (((ix * ny) + iy) * nz + iz) * nxi + ixi
```

Shape is derived from axes, preventing disagreement with separate shape metadata.
The constructor checks element and byte counts (including the Rust `isize::MAX`
allocation limit), exact value count, and finite nonnegative values. It consumes
the vector without copying. A future file reader must check resource limits
before allocating; this in-memory constructor receives already allocated data.
All fields are private and exposed through immutable accessors. Invalid node
indices return `None`. Errors identify the offending axis or flat scalar index.

## Implemented: stored-coordinate sampling

`EnvelopeGrid::sample_stored(point: [f64; 4])` returns
`Result<StoredEnvelopeSample, SampleError>`:

```rust
let sample = grid.sample_stored([x, y, z, xi])?;
let s = sample.a_sqd;
let [s_x, s_y, s_z, s_xi] = sample.derivatives;
```

The query is already in `(x,y,z,xi)` coordinates, in metres. Derivatives hold
the other stored coordinates fixed and have units of inverse metres. They are
not a laboratory four-vector. S is already polarization-averaged and is not
rescaled. Sampling borrows existing validated storage without copying arrays or
allocating on the heap; the enclosing cell's 16 values use stack storage.

`SampleError` is separate from malformed-grid `GridError`. Stored query variants are
`NonFiniteCoordinate { axis, coordinate }` and
`OutOfDomain { axis, coordinate, min, max }`, with axis indices in `(x,y,z,xi)`
order and axis names in error messages. If several axes are invalid, the first
in this order is reported.

For each cell, `u=(q-q_lower)/spacing` defines lower/upper weights `1-u,u`.
The value is the sum of the 16 corner values times four weights. Implementation
uses nested linear combinations, algebraically the same product-weight sum.
For each derivative, replacing its axis weight with `-1/spacing,+1/spacing`
gives eight weighted edge differences. Pairing the terms before division removes
a common offset and makes constant-data derivatives exactly zero. Derivatives
use the same corners and fractions as the value, not finite-differenced samples
or separately interpolated derivative arrays.

The canonical domain is closed. Bounds are checked with strict comparisons
before any division; nonfinite coordinates and out-of-domain queries are rejected.
At an exact upper bound use the last cell with fraction one. Otherwise estimate
the cell by division, then check and adjust against canonical node coordinates.
An exact interior knot uses the cell on its positive side at fraction zero;
neighbouring representable queries retain their own cells. There is no tolerance
that expands coverage or snaps query coordinates to knots.

After strict bounds and cell checks, the computed fraction is capped at one
only if canonical reconstruction/subtraction/division rounded it above one.
This arithmetic safeguard changes neither query coordinates nor cell selection.
It prevents negative lower weights from roundoff; it is not an out-of-domain
clamp. The exact upper endpoint is handled explicitly rather than relying on
division to produce one. Floating-point reconstruction and weight arithmetic
still limit endpoint/continuity accuracy; no bit-exact continuity is promised.

A supplied endpoint equal to a canonical bound is a boundary sample. One lying
inside is an interior sample; one outside returns `OutOfDomain`, even if the
difference arose during reconstruction. Tests cover inward/outward upper-bound
movement and adjacent representable queries on every axis.

Multilinear values are continuous in exact arithmetic and gradients generally
jump at cell boundaries. Nonnegative nodes give convex interpolated values up
to floating-point arithmetic. Extremely large values or very small spacings
can overflow derivative arithmetic despite finite input validation; these tests
cover well-conditioned physical scales, not arbitrary f64 extremes. No claim of
Gaussian convergence, trajectory accuracy or production suitability follows
from exact-function tests. Coverage does not determine particle termination.

## Implemented: laboratory sampling

`EnvelopeGrid::sample_lab(r: FourVector) -> Result<LabEnvelopeSample, SampleError>`
accepts `(ct,x,y,z)` in metres and delegates to `sample_stored([x,y,z,ct-z])`.
It does not duplicate or modify interpolation. The result type is distinct from
`StoredEnvelopeSample`:

```rust
let r = FourVector::new(ct, x, y, z);
let sample = grid.sample_lab(r)?;
let s = sample.a_sqd;
let grad = sample.grad_a_sqd; // (d/d(ct), -d/dx, -d/dy, -d/dz) S
```

The chain rule at fixed laboratory position gives `dS/d(ct)=S_xi`. Varying
laboratory z at fixed ct changes both stored z and xi, giving `dS/dz=S_z-S_xi`.
Raising the derivative index with the `(+---)` metric gives
`(S_xi,-S_x,-S_y,-S_z+S_xi)`. All components are inverse metres; the first is a
ct derivative, not a time derivative, so no factor of c is introduced. S remains
dimensionless and already cycle-averaged. This preserves independent stored-z
evolution rather than imposing a plane-wave relation.

The API validates every laboratory coordinate before transforming or sampling.
Errors distinguish:

- `NonFiniteLabCoordinate { component, coordinate }`: invalid input, with
  component order `(ct,x,y,z)` and laboratory names in messages.
- `NonFiniteXi { ct, z }`: overflow of `ct-z` despite finite operands.
- `TransformedOutOfDomain { position, axis, coordinate, min, max }`: the original
  laboratory position maps outside the canonical stored domain. Axis order here
  is `(x,y,z,xi)`; messages explicitly identify xi as derived from ct-z.
- `NonFiniteResult { quantity, value }`: a nonfinite scalar, stored derivative
  or raised-gradient component. In particular, finite stored derivatives can
  overflow in `-S_z+S_xi`. Checks run before returning a laboratory sample.

No arbitrary amplitude/gradient thresholds are imposed. Stored sampling's existing
arithmetic is unchanged; the laboratory wrapper rejects nonfinite results rather
than repairing them. Rounding in ct-z follows ordinary f64 arithmetic, and its
result must satisfy strict canonical bounds. There is no domain expansion,
position clamp, zero padding, or particle-termination rule. The wrapper inherits
the stored interpolant's knot policy and discontinuous cellwise derivatives.

## Physics conventions

Stored values are dimensionless `S = a_rms^2`. Cycle averaging is already applied:
for the intended peak-amplitude convention, linear polarization has
`S_peak = a0^2 / 2`, circular polarization `S_peak = a0^2`. Metadata carries
polarization and a positive finite reference wavelength in metres. Construction
does not rescale S. Coordinate and normalization conventions are fixed by the
type contract; unsupported alternatives cannot be selected in memory.

Ptarmigan positions are `(ct,x,y,z)` in metres, with metric `(+---)`. The implemented
sampling coordinate is `(x,y,z,ct-z)`. The chain rule gives the raised gradient

```text
grad_a_sqd = (S_xi, -S_x, -S_y, -S_z + S_xi)  [1/m]
```

Here `S_z` holds stored xi fixed. It must not be eliminated using a plane-wave
relation. This matches the LMA pusher's `c * grad_a_sqd / 2` force convention.
The future constant axial wavevector is physical `2*pi/lambda * (1,0,0,1)` in
inverse metres; event routines use `kappa = c * COMPTON_TIME * k` instead.
Laboratory-gradient conversion is implemented; wavevector evaluation remains pending.

The native Gaussian time-gradient discrepancy is separate from interpolation.
Source inspection predicts `g_native[0] - S_xi = -S_z`: native `grad_a_sqd`
includes independent longitudinal beam evolution in its time component.
The user's completed diagnostic in `output/envelope-native` numerically confirms
a discrepancy at the tested points. At on-axis z=z_R/2, xi=0, fixed-position
finite differences and the analytical time derivative are zero, while the linear
native time component is about 3259.493 /m. Trajectory effects remain unestablished.
These existing outputs were inspected, not regenerated in the cubic task; native
physics is unchanged. See the validation guide for scope and conventions.

## Remaining milestone 1 work

1. Versioned HDF5 input behind proposed `hdf5-input = ["dep:hdf5-writer"]`.

Proposed v1 HDF5 group `/envelope` contains scalar datasets `format`, `version`,
`coordinates`, `axis_order`, `storage_order`, `normalization`, `polarization`,
and `reference_wavelength`; rank-1 `axes/{x,y,z,xi}`; and rank-4
`a_rms_squared` with shape `[nx,ny,nz,nxi]`. Conventions are respectively
`ptarmigan-envelope`, unsigned 32-bit `1`, `x,y,z,xi;xi=ct-z`, `x,y,z,xi`,
`C;xi-fastest`, and `cycle-averaged-a-rms-squared`. Polarization is `linear` or
`circular`. Numeric data are little-endian float64 except the version. Unit
attributes are `m` for axes/wavelength and `1` for S. Validate all required
metadata, rank, shape, datatype and values. Ignore unrelated datasets; reject
unknown versions. Future complex data do not imply carrier-resolved fields.

Reuse the local HDF5 reader. Its array read scatters the first dimension over
the communicator, so the initial loader must require a single-rank communicator.
No new HDF5 binding, MPI distribution scheme, or CLI input is needed yet.

## Implemented Gaussian validation tooling (test-only)

`src/field/envelope/validation/` is compiled only with `cfg(test)`. It separates
the reference (`gaussian.rs`), canonical grid builder (`grid.rs`), error statistics
(`statistics.rs`), CSV/metadata generation (`report.rs`), configurable study
(`runner.rs`) and native diagnostic (`native.rs`). Ignored tests expose the two
runners without changing the binary interface or adding dependencies.

The independent stationary Gaussian reference uses
`X=x-x_c`, `Y=y-y_c`, `Z=z-z_f`, `Xi=xi-xi_c`, `W=1+Z^2/z_R^2`,
`z_R=pi*w0^2/lambda`, `L=c*T`, with T the intensity FWHM:

```text
S = C_pol*a0^2/W * exp(-2*(X^2+Y^2)/(w0^2*W) - 4*ln(2)*Xi^2/L^2)
S_x  = -4*X*S/(w0^2*W)
S_y  = -4*Y*S/(w0^2*W)
S_xi = -8*ln(2)*Xi*S/L^2
S_z  = 2*Z/(z_R^2*W) * (2*(X^2+Y^2)/(w0^2*W)-1) * S
```

The reference independently constructs `(S_xi,-S_x,-S_y,-S_z+S_xi)` without
calling the production laboratory conversion. Its scalar is independently
finite-differenced in small tests. Source inspection confirms native factors
`C_pol=1/2` for linear and `1` for circular, with `n_cycles=c*T/lambda` matching
the native Gaussian intensity duration. The reference requires positive finite
wavelength, waist, duration and a0, finite offsets and representable derived
scales. Zero a0 is excluded specifically because normalized validation requires
a positive S_peak; this is not a new production-field restriction.

`xi_c` is a retarded-coordinate offset, not a laboratory time. The pulse peak at
the focus occurs at `ct=z_f+xi_c`. Moving z_f alone does not generally keep that
peak at ct=0; choose `xi_c=-z_f` to do so.

Grid values are evaluated at canonical reconstructed coordinates. Count/byte
arithmetic and a configurable scalar-memory budget are checked before allocation;
estimated scalar storage is printed before construction. Each grid is dropped
before constructing the next resolution. The budget covers the scalar Vec, not
process peak memory, coordinate vectors or the reusable query list.

The convergence runner generates deterministic interior points once using the
existing Xoshiro256StarStar RNG. Candidates lie in 95% of the configured domain,
exclude symmetry planes and knots of every requested resolution, and are reused
for both polarizations. It round-trips through ct=z+xi then xi=ct-z and uses that
effective stored point for the analytical and numerical evaluations. Production
domain rules are unchanged. `queries.csv` records the actual physical queries.

`errors.csv` reports max/RMS absolute and characteristic-scale-normalized errors
for S, four stored derivatives and four laboratory components. Rates use the
actual successive cell-count ratio; zero/nonfinite errors or missing previous
data give `NA`. Nonfinite measurements abort a study instead of producing a
misleading success. Expected multilinear asymptotic orders (value 2, gradients 1) are not
assertions. Provisional maximum normalized targets are 1e-3 for S and 1e-2 for
gradients. Misses are reported without automatic refinement.

The separate native diagnostic uses an unshifted reference, both polarizations,
on-axis z=z_R/2 at xi=0, and a generic nonzero-xi/off-axis point. It writes scalar,
full-gradient and native-scalar finite-difference comparisons at three steps.
Disagreement is diagnostic output, not a failing regression assertion. This
diagnostic was run separately by the user; its existing output confirms the
time-gradient discrepancy at the tested points. It was not rerun for cubic work.

Output includes human-readable summaries, CSVs, and `metadata.txt` containing
runtime Git commit/dirty state, compile-time Git/features, physical parameters,
domain/offsets, seed/algorithm, normalization, allocation estimates and diagnostic
identity. Status is appended: only a final `run_status=completed` indicates success.
Files are created exclusively; reuse requires a fresh output directory. Default
outputs live under the already ignored `/output/`. No output data is committed.

See [validation commands and configuration](numerical-envelope-validation.md#gaussian-tooling-commands)
and [progress](numerical-envelope-progress.md) for actual runs. The existing
multilinear 16/32/64 study shows expected convergence but misses sampled targets.
The full cubic comparison and measured sampling cost/peak memory remain pending.

## Optional shared-slope cubic interpolation

`InterpolationMethod::{Multilinear,Cubic}` selects
`sample_stored_with_method(q, method)` or `sample_lab_with_method(r, method)`.
Existing `sample_stored(q)` and `sample_lab(r)` remain multilinear. Both methods
share the laboratory transform and finite-result checks; derivatives retain the
stored/laboratory distinction and the same normalization. The cubic stored sampler
also rejects nonfinite results. `cubic.rs` contains the separate implementation.

On a cell, u=(q-q_i)/h, use

```text
P = H00 f_i + H10 h m_i + H01 f_(i+1) + H11 h m_(i+1)
H00=2u^3-3u^2+1; H10=u^3-2u^2+u
H01=-2u^3+3u^2; H11=u^3-u^2
m_i=(f_(i+1)-f_(i-1))/(2h)                 interior
m_0=(-3f_0+4f_1-f_2)/(2h)                 lower endpoint
m_N=(3f_N-4f_(N-1)+f_(N-2))/(2h)          upper endpoint
```

These slopes belong to shared nodes, not individual cells. Value weights and
analytically differentiated weights define a four-dimensional tensor product;
mixed terms follow from the linear operators without derivative arrays. The
implementation subtracts a common stencil value before applying weights, reducing
common-offset cancellation and preserving constants exactly. It uses stack arrays,
no per-sample heap allocation and no coefficient grids. Each axis uses at most
four nodes (three in boundary cells); cubic requires at least three nodes on every
axis or returns `UnsupportedCubicAxis`. Two-node grids remain valid for multilinear.
There is no fallback. Nominal distinct nodal reads grow from 16 to at most 256;
this is not a measured runtime ratio.

The canonical locator, strict domain and positive-side knot convention are shared
unchanged. The upper endpoint selects the final cell. No extrapolation or domain
expansion is introduced. As in multilinear sampling, the locator caps an in-domain
fraction that rounds above one; it never admits or moves an outside coordinate.
Canonical FMA reconstruction can slightly change node separations from nominal h.
Shared slopes give C1 continuity in exact arithmetic; numerical continuity and
quadratic reproduction are expected only to floating-point accuracy on
well-conditioned grids, not bit-exact equality. Centred approximate slopes generally
limit smooth-function convergence to about order three for values and two for
gradients, not four/three. Those orders remain expectations for the pending study.

**Signed overshoot is exposed unchanged.** Nonnegative nodes can produce negative
S; no clipping, limiter, log interpolation or hidden method switch is applied.
The small fixture `[0,0,1,1]` gives S=-0.125 at x=0.5 and S=-0.09375,
S_x=-0.25 at x=0.25 (unit spacing). This option is for diagnostic evaluation and
is not ready for particle/radiation use until positivity is addressed separately.

The Gaussian runner evaluates both methods on each grid before releasing it.
`ENVELOPE_QUERIES` optionally loads an existing `queries.csv`; otherwise the
original deterministic generator is retained. Loading requires exactly one each
of `x_m,y_m,z_m,xi_m,ct_m` (column order may vary), finite numeric values and exact
numerical `xi == ct-z` after f64 parsing, with signed zeros equal. The writer's
17-decimal scientific representation round-trips f64. Rounded/truncated inconsistent
external CSVs fail; points are never adjusted, dropped or regenerated. Every
point must lie inside the closed canonical domain at every requested resolution.
Loaded knots/symmetry planes are retained. Metadata records the source and actual
count; configured seed/count do not override a loaded file.

`errors.csv` has a method column and independent per-method/per-polarization rates.
Maxima and target decisions explicitly refer to sampled points. `samples.csv`
separately records sampled minimum S/S_peak, negative count and maximum negative
excursion `max(0,-min(S/S_peak))`. These diagnostics establish no full-domain
positivity guarantee. Nonfinite values, derivatives or normalized errors abort.

Tracking, radiation, LASY propagation/import, flying-focus optics, general local
wavevectors and LCFA remain out of scope. Scalar S cannot reconstruct consistent
carrier-resolved E and B. Later matched-history benchmarks isolate local radiation
calculations, not spatial gradients or finite-beam dynamics.

See [validation](numerical-envelope-validation.md) and
[progress](numerical-envelope-progress.md).
