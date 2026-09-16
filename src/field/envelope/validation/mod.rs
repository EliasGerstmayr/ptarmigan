//! Test-only analytical references and opt-in diagnostics; no simulation integration.
mod gaussian;
mod grid;
mod native;
mod report;
mod runner;
mod statistics;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[test]
#[ignore = "opt-in Gaussian study; configured with ENVELOPE_* environment variables"]
fn gaussian_convergence() {
    runner::run().expect("Gaussian validation runner failed");
}

#[test]
#[ignore = "opt-in native field diagnostic, independent of regression pass/fail"]
fn native_gaussian_diagnostic() {
    native::run().expect("native Gaussian diagnostic failed");
}
