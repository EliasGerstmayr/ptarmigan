# Numerical envelope design

## Implemented: in-memory grid contract

`src/field/envelope/` holds an envelope-specific numerical core, independent of
HDF5 and of the `Field` trait. Existing field implementations are unchanged.
Only validated construction, read-only metadata/storage access, and node indexing
are implemented. There is no interpolated sampling or simulation input yet.

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
endpoint through floating-point rounding. Future interpolation bounds will be
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

## Physics conventions

Stored values are dimensionless `S = a_rms^2`. Cycle averaging is already applied:
for the intended peak-amplitude convention, linear polarization has
`S_peak = a0^2 / 2`, circular polarization `S_peak = a0^2`. Metadata carries
polarization and a positive finite reference wavelength in metres. Construction
does not rescale S. Coordinate and normalization conventions are fixed by the
type contract; unsupported alternatives cannot be selected in memory.

Ptarmigan positions are `(ct,x,y,z)` in metres, with metric `(+---)`. The planned
sampling coordinate is `(x,y,z,ct-z)`. The chain rule gives the raised gradient

```text
grad_a_sqd = (S_xi, -S_x, -S_y, -S_z + S_xi)  [1/m]
```

Here `S_z` holds stored xi fixed. It must not be eliminated using a plane-wave
relation. This matches the LMA pusher's `c * grad_a_sqd / 2` force convention.
The future constant axial wavevector is physical `2*pi/lambda * (1,0,0,1)` in
inverse metres; event routines use `kappa = c * COMPTON_TIME * k` instead.
Neither gradient conversion nor wavevector evaluation is implemented in step 1.

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

1. Multilinear interpolation of the same 16 corners for value and derivatives.
2. Laboratory-coordinate conversion, tested with independent stored-z evolution.
3. Versioned HDF5 input behind proposed `hdf5-input = ["dep:hdf5-writer"]`.

Planned coverage is the closed rectangular domain of canonical coordinates. At the upper endpoint use the
last cell; at interior knots choose the cell on the positive side. Derivatives
are one-sided there. Outside and nonfinite queries return explicit errors;
coverage does not determine particle termination. These policies are not yet
implemented by the node-only API.

Test supplied and canonical endpoints separately in the future interpolator.
A supplied endpoint equal to a canonical bound is a boundary sample. One lying
strictly inside is an ordinary interior sample; one outside is out of domain,
even if the difference arose during reconstruction. Test immediately adjacent
representable coordinates on both sides of the canonical bounds as well. The
axis-construction tolerance is not a sampling tolerance: do not widen coverage
or silently clamp genuinely out-of-domain coordinates.

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
