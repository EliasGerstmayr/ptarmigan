# Numerical envelope progress

## 2026-09-15: step 1

Implemented validated in-memory axes, metadata, scalar storage, node indexing and
errors in `src/field/envelope/{mod,grid,error}.rs`. Registered the module in
`src/field/mod.rs`. Added design, validation and progress documentation. No
dependencies or existing field implementations changed.

Branch at start: `numerical-field-lma`, commit `9e7caf7`, clean working tree.
The sibling `ptarmigan-baseline` worktree remains outside the implementation scope.
No commit or push requested or performed.

Verification completed. Cargo 1.98.1 is available through the direct stable
toolchain binary. During that run the cargo proxy reported rustup instead;
the historical workaround is recorded below. The first sandboxed Cargo test attempt encountered DNS
failures downloading locked dependencies (exit 101); the escalated retry
downloaded the dependencies and passed. No dependency versions or lockfile changed.

Actual checks and results:

| Check | Result |
| --- | --- |
| Focused release tests, no default features | 6 passed, 0 failed; 132 unrelated tests filtered out. Repeated after the final warning-annotation change: same result. |
| Default release build, no default features | Passed, exit 0. |
| rustfmt check on all three new Rust files | Passed, exit 0. |
| git diff --check | Passed, exit 0. |
| Read-only baseline status | Clean detached HEAD; no baseline edits. |

Build/test output includes warnings in existing code and a nonfatal toolchain
warning: rust-objcopy cannot load `libLLVM.dylib` while stripping debug information.
Both the executable build and focused test process completed successfully.
No full-suite, HDF5/MPI, physics convergence or sampling benchmark runs were made.
The six tests verify the grid contract only; derivative signs, Gaussian accuracy,
convergence orders and sampling performance remain planned, unverified behaviour.

Historical commands below reproduce the checks using the same direct toolchain as the run
(the actual Cargo invocations used its absolute path and inline RUSTC):

```bash
export PATH="/Users/eliasgerstmayr/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"
export RUSTC="/Users/eliasgerstmayr/.rustup/toolchains/stable-aarch64-apple-darwin/bin/rustc"
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
cargo build --locked --release -p ptarmigan --no-default-features
rustfmt --edition 2018 --check src/field/envelope/mod.rs src/field/envelope/grid.rs src/field/envelope/error.rs
git diff --check
```

Next bounded task: multilinear value and grid-coordinate derivatives from one
interpolant. Milestone 1 is not complete: interpolation, coordinate conversion
and file input remain. Milestone 2 Gaussian and performance gates are pending.

Suggested commit message: `feat(field): add validated in-memory envelope grid`.

## Step 1 review follow-up

Inspected the working tree before editing and preserved all pending step 1 changes.
This follow-up changes documentation and one small grid regression test only.

- Rechecked native Gaussian source; the gradient finding is a source-level
  discrepancy awaiting numerical verification. Documented a future native-scalar
  ct finite-difference diagnostic at nonzero z and the temporal peak, with fixed
  laboratory x/y/z. Neither the diagnostic nor native physics was implemented.
- Made canonical coordinates authoritative for future interpolation coverage;
  supplied endpoints are classified relative to canonical bounds without clamping.
  Added a rounded-input regression test that actually exercises endpoint movement.
- Added portable macOS/Homebrew setup to [build instructions](build.md#macos-with-homebrew).
  Historical commands/results above are retained. HDF5_DIR applies to HDF5-enabled
  builds, not focused storage tests. No dependencies or Cargo.lock changed.

Review verification results:

| Check | Result |
| --- | --- |
| Focused release storage tests | 7 matched, 7 passed, 0 failed; 132 unrelated tests filtered out. |
| rustfmt check on all three envelope Rust files | Passed. |
| git diff --check and explicit trailing-whitespace check of follow-up files (including untracked files) | Passed. |
| Diff check of Cargo.toml, Cargo.lock and native focused-field implementations | No changes. |
| Read-only baseline status | Clean detached HEAD. |

Actual test and formatting commands for this follow-up:

```bash
export PATH="$(brew --prefix rustup)/bin:$PATH"
rustfmt --edition 2018 src/field/envelope/grid.rs
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
rustfmt --edition 2018 --check src/field/envelope/mod.rs src/field/envelope/grid.rs src/field/envelope/error.rs
git diff --check
git diff --exit-code -- Cargo.lock Cargo.toml src/field/focused_laser.rs src/field/fast_focused_laser.rs
git --no-optional-locks -C ../ptarmigan-baseline status --short --branch
```

The run reported existing code warnings and the previously recorded nonfatal
rust-objcopy/libLLVM build-script warning. No native Gaussian numerical diagnostic,
trajectory test, interpolation test, full-suite test or HDF5 build was run in
this follow-up. The default release build success above remains historical evidence.
The regression verifies reconstruction, not future continuous-domain sampling.

Step 1 is ready to commit within its bounded storage-only scope; no commit or push
has been performed. Suggested message: `feat(field): add validated in-memory envelope grid`.

Next task: four-dimensional multilinear interpolation returning S and all four
stored-coordinate derivatives from the same interpolant. Laboratory-gradient
conversion remains a separate subsequent step.

## Stored-coordinate multilinear interpolation

Started from clean `numerical-field-lma` at `72d0d70` (step 1), tracking the fork's
branch. Inspected the envelope code, design/validation/progress documents and
repository guidance; no applicable AGENTS.md was found. Baseline remains at
`9e7caf7` and is not modified.

Changes:

- Added `multilinear.rs` with `sample_stored([x,y,z,xi])`, returning dimensionless
  `a_sqd` and `[S_x,S_y,S_z,S_xi]` in inverse metres, with other stored coordinates
  fixed. Sampling uses 16 stack-resident corner values and existing grid storage.
- Added distinct sampling errors and exported the sample/error types. Updated
  the axis contract comment from future to implemented sampling bounds.
- Enforced strict canonical bounds and canonical-node verification of estimated
  cells. Exact knots select the positive side; the upper endpoint uses the last
  cell. Any post-check fraction cap corrects arithmetic roundoff only.
- Added nine interpolation tests and updated design/validation documentation.
  Native-gradient numerical verification remains an independently pending issue.

The first focused run matched 16 tests: 15 passed, one failed because the proposed
inward-endpoint test fixture did not actually move inward. Replaced that fixture
with independently checked rounded decimal coordinates and corrected a Rust 2018
test-message formatting warning. This failure did not require a sampler change.

Final verification:

| Check | Result |
| --- | --- |
| Focused release envelope tests | 16 matched, 16 passed, 0 failed; 132 unrelated tests filtered out. Includes 7 grid and 9 sampling tests. |
| Default package release build | Passed, exit 0. |
| rustfmt check on all four envelope Rust files | Passed. |
| git diff --check | Passed. |
| Explicit whitespace check including new untracked sampler | Passed. |
| Cargo.toml, Cargo.lock, native focused fields and existing numerical fields | No changes. |
| Read-only baseline status | Clean detached HEAD at 9e7caf7. |

Actual commands (Homebrew rustup; no HDF5 feature or HDF5_DIR required):

```bash
export PATH="$(brew --prefix rustup)/bin:$PATH"
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
cargo build --locked --release -p ptarmigan --no-default-features
rustfmt --edition 2018 --check src/field/envelope/mod.rs src/field/envelope/grid.rs src/field/envelope/error.rs src/field/envelope/multilinear.rs
git diff --check
git diff --exit-code -- Cargo.lock Cargo.toml src/field/focused_laser.rs src/field/fast_focused_laser.rs src/field/numerical/
git --no-optional-locks -C ../ptarmigan-baseline status --short --branch
```

Output retains existing code warnings and the historical nonfatal build-script
rust-objcopy/libLLVM warning. No check was blocked. No full-suite, HDF5/MPI,
Gaussian, trajectory or performance checks were run. Exact-function tests verify
interpolation identities and boundary semantics, not demonstrated physics accuracy
or convergence. Sampling allocates no heap storage by source inspection; allocation
instrumentation and runtime cost measurements were not performed.

Suggested commit message: `feat(field): add multilinear stored-envelope sampling`.

No laboratory conversion, native physics change, HDF5, particle tracking,
radiation, Gaussian convergence, benchmark, cubic interpolation or wavevector
evaluation is included. No dependencies or lockfile changed; no commit or push.
Next bounded step is laboratory-coordinate sampling and raised-gradient conversion,
after review of this interpolation diff.

## Laboratory-coordinate sampling and raised gradient

Started from clean `numerical-field-lma` at `8f14e1f`: the stored sampler had
been committed. Inspected the envelope implementation, documentation, worktrees
and repository instructions; no applicable AGENTS.md was found. Verified the
`(+---)` metric in `FourVector::Mul` and metre-valued `(ct,x,y,z)` positions from
the existing `c*u*dtau` position updates. Baseline remains at `9e7caf7`.

Added `laboratory.rs` and exported `LabEnvelopeSample`. `sample_lab` reuses the
stored interpolant at `[x,y,z,ct-z]` and returns `(S_xi,-S_x,-S_y,-S_z+S_xi)`.
Extended `SampleError` to distinguish invalid laboratory inputs, transformed xi
overflow, transformed-domain failures, and nonfinite output arithmetic. Updated
design and validation documentation, including the explicitly pending benchmark
checklist. No interpolation arithmetic or existing field physics was changed.

Verification completed:

| Check | Actual result |
| --- | --- |
| Focused release envelope tests | 26 matched, 26 passed, 0 failed; 132 unrelated tests filtered out. Includes 7 grid, 9 stored sampling and 10 laboratory tests. |
| Default package release build | Passed, exit 0. |
| rustfmt check on all five envelope Rust files | Passed. |
| git diff --check and explicit whitespace check including the new file | Passed. |
| Dependency/lockfile, native field and stored interpolation diffs | Unchanged. |
| Read-only baseline status | Clean detached HEAD. |

The laboratory tests include independent finite differences at three steps,
coordinate mappings and finite-input overflow cases. Synthetic nonfinite samples
also test the result guard directly. Tests use mathematical fixtures, not particle
simulations. All requested checks ran; none was blocked. Output includes existing
code warnings and the historical nonfatal rust-objcopy/libLLVM build-script warning.

Actual verification commands:

```bash
export PATH="$(brew --prefix rustup)/bin:$PATH"
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
cargo build --locked --release -p ptarmigan --no-default-features
rustfmt --edition 2018 --check src/field/envelope/mod.rs src/field/envelope/grid.rs src/field/envelope/error.rs src/field/envelope/multilinear.rs src/field/envelope/laboratory.rs
git diff --check
git diff --exit-code -- Cargo.lock Cargo.toml src/field/focused_laser.rs src/field/fast_focused_laser.rs src/field/numerical/ src/field/envelope/grid.rs src/field/envelope/multilinear.rs
git --no-optional-locks -C ../ptarmigan-baseline status --short --branch
```

No physical simulation, large convergence study or benchmark was launched. Native
Gaussian discrepancy verification remains pending, as do the Gaussian, memory/cost,
trajectory, radiation, LASY, flying-focus, LCFA and finite-beam benchmarks listed
in the validation document. No dependency or lockfile changes, commit or push.

Suggested commit message: `feat(field): add laboratory envelope sampling and raised gradient`.

## 2026-09-16: analytical Gaussian validation tooling

Started from clean `numerical-field-lma` at `e670c3f`. The laboratory sampler was
already committed. Inspected envelope source, design/validation/progress documents,
Cargo/features and native Gaussian source; no applicable AGENTS.md was found.
Baseline remains at `9e7caf7`. Production samplers, native physics, dependencies
and Cargo.lock are unchanged.

Added test-only `validation/{mod,gaussian,grid,statistics,runner,report,native}.rs`
and its cfg(test) registration. The reference supports independent transverse,
focus and retarded-time offsets; normalized validation requires positive a0.
Native source confirms linear/circular averaging 1/2 and 1 and the duration
mapping n_cycles=c*T/lambda. No native diagnostic was run to verify its pending
time-gradient discrepancy numerically.

Implemented two ignored runners with ENVELOPE_* configuration, exclusive output
files, CSVs plus metadata, deterministic query reuse, canonical grid evaluation,
checked memory budget and sequential grid release. Full study and diagnostic
commands are documented in the validation guide. Expected orders and engineering
targets are reports, not forced assertions or automatic refinement triggers.

Actual checks:

| Check | Result |
| --- | --- |
| Focused release envelope tests | 37 matched: 35 passed, 0 failed, 2 ignored; 132 unrelated tests filtered out. Includes 9 new small tests. |
| Tiny Gaussian runner smoke | 1 matched, 1 passed; 168 filtered out. Native diagnostic not invoked. |
| Default package release build | Passed, exit 0. |
| CSV/metadata serialization inspection | Passed: 54 error rows, 24 axis rows, 16 queries, finite metrics, NA first rates and actual ratios 1.5 and 4/3. Metadata records completion and native-not-run. |
| Formatting and whitespace | rustfmt --check and git diff --check passed; explicit whitespace check includes every new validation file. |
| Production/dependency/lockfile diff and baseline status | Production samplers, native fields and dependencies unchanged; baseline clean. |

Smoke invocation used cells=2,3,4, samples=16, seed=20260916, memory ceiling=1 MiB,
default physical parameters (lambda=0.8 um, w0=5 um, T=30 fs, a0=1, zero offsets),
extent=2 and both polarizations. Largest scalar allocation was 5,000 bytes.
Output directory: `/private/tmp/ptarmigan-gaussian-smoke-20260916`. This was solely
an invocation/serialization smoke; missed engineering targets and its measured
coarse rates are not production accuracy/convergence evidence. No automatic
refinement followed. The report's exact runtime source state is captured in its
metadata; documentation was finalized after the smoke.

Actual runner invocation:

```bash
export PATH="$(brew --prefix rustup)/bin:$PATH"
ENVELOPE_CELLS=2,3,4 ENVELOPE_SAMPLES=16 ENVELOPE_MAX_MIB=1 ENVELOPE_SEED=20260916 ENVELOPE_OUTPUT=/private/tmp/ptarmigan-gaussian-smoke-20260916 cargo test --locked --release -p ptarmigan --no-default-features field::envelope::validation::gaussian_convergence -- --ignored --exact --nocapture --test-threads=1
cargo test --locked --release -p ptarmigan --no-default-features field::envelope:: -- --test-threads=1
cargo build --locked --release -p ptarmigan --no-default-features
rustfmt --edition 2018 --check src/field/envelope/mod.rs src/field/envelope/validation/*.rs
git diff --check
```

Existing code warnings and the nonfatal rust-objcopy/libLLVM build-script warning
remain. No requested check was blocked. Full 16/32/64-cell study, native diagnostic, performance/peak-memory
measurements and all particle/physics simulations were not run. No commit or push.

Suggested commit message: `test(field): add analytical Gaussian validation tools`.
