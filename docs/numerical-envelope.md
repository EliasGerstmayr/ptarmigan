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

Reinspection of `FocusedLaser::a_sqd`, `envelope_and_grad` and `grad_a_sqd`
identifies a **source-level discrepancy awaiting numerical verification**.
For `S = B(x,y,z) F(xi)`, varying ct at fixed laboratory x, y and z gives
`dS/d(ct) = B F' = S_xi`. The native gradient's time component also includes
`-grad_beam[2] * envelope`, apparently contributing `-S_z` from independent
longitudinal beam evolution. The source therefore suggests a time-component
difference of `-S_z`; this is not a numerically validated bug, and its effect on
trajectories has not been established. A future diagnostic must finite-difference
the native scalar at fixed laboratory position and compare it to the native
gradient (see validation). Do not alter native physics or force agreement here.

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

## Milestone 2 and later

Generate an independent analytical Gaussian using `W=1+z^2/z_R^2`,
`z_R=pi*w0^2/lambda`, and intensity FWHM T:

```text
S = C_pol*a0^2/W * exp(-2*(x^2+y^2)/(w0^2*W) - 4*ln(2)*xi^2/(c*T)^2)
S_x  = -4*x*S/(w0^2*W)
S_y  = -4*y*S/(w0^2*W)
S_xi = -8*ln(2)*xi*S/(c*T)^2
S_z  = 2*z/(z_R^2*W) * (2*(x^2+y^2)/(w0^2*W)-1) * S
```

Use native Gaussian evaluation only as an additional comparison, with
`n_cycles=c*T/lambda`. Measure three-resolution convergence and sampling cost.
Expected multilinear orders are two for values and one for gradients; they are
not measured results yet. Scalar storage is `8*nx*ny*nz*nxi` bytes.

Defer smooth tensor-product cubic interpolation until after the reference is
validated. A shared-slope cubic candidate uses four nodes per axis (256 samples),
with an explicit one-node halo or subsequently validated boundary slopes.
On-demand coefficients avoid persistent derivative arrays. Overshoot can cause
negative S; never clamp a value while retaining the original gradient.

Tracking, radiation, LASY propagation/import, flying-focus optics, general local
wavevectors and LCFA remain out of scope. Scalar S cannot reconstruct consistent
carrier-resolved E and B. Later matched-history benchmarks isolate local radiation
calculations, not spatial gradients or finite-beam dynamics.

See [validation](numerical-envelope-validation.md) and
[progress](numerical-envelope-progress.md).
