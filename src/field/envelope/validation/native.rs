use super::{
    gaussian::{polarization_name, Gaussian},
    report::Report,
    runner::Config,
    Result,
};
use crate::{
    constants::SPEED_OF_LIGHT,
    field::{Envelope, FocusedLaser, Polarization},
};
use std::io::Write;

/// Diagnostic only. Numerical disagreement is output data, never a failure gate.
pub fn run() -> Result<()> {
    let c = Config::from_env("native")?;
    if c.params.offsets != [0.0; 4] {
        return Err("native diagnostic requires all offsets zero to match FocusedLaser".into());
    }
    let mut report = Report::new(&c, "native_gaussian_diagnostic")?;
    let mut f = report.file("native.csv")?;
    writeln!(f,"polarization,point,ct_m,x_m,y_m,z_m,component,h_m,native_scalar,analytical_scalar,scalar_error_normalized,native_gradient,analytical_gradient,native_scalar_fd,native_minus_analytic_normalized,native_minus_fd_normalized")?;
    for pol in [Polarization::Linear, Polarization::Circular] {
        let g = Gaussian::new(c.params, pol)?;
        report.gaussian(&g)?;
        // FocusedLaser duration = n_cycles * wavelength / c, Gaussian intensity FWHM.
        let native = FocusedLaser::new(
            c.params.a0,
            c.params.wavelength,
            c.params.waist,
            SPEED_OF_LIGHT * c.params.duration / c.params.wavelength,
            pol,
            0.0,
        )
        .with_envelope(Envelope::Gaussian);
        let z = 0.5 * g.zr;
        for (label, r) in [
            ("on_axis_temporal_peak", [z, 0.0, 0.0, z]),
            (
                "generic",
                [
                    z + 0.23 * g.length,
                    0.31 * c.params.waist,
                    -0.27 * c.params.waist,
                    z,
                ],
            ),
        ] {
            let s = native.a_sqd(r.into());
            let exact = g.lab(r);
            let gradient = native.grad_a_sqd(r.into());
            let steps = [g.length, c.params.waist, c.params.waist, g.length.min(g.zr)];
            println!(
                "{} {}: scalar normalized difference={:.3e}; time native={:.6e} analytic={:.6e}",
                polarization_name(pol),
                label,
                (s - exact[0]) / g.peak,
                gradient[0],
                exact[1]
            );
            for epsilon in [1e-3, 1e-4, 1e-5] {
                for i in 0..4 {
                    let h = epsilon * steps[i];
                    let mut a = r;
                    let mut b = r;
                    a[i] += h;
                    b[i] -= h;
                    let fd = (native.a_sqd(a.into()) - native.a_sqd(b.into())) / (2.0 * h)
                        * if i == 0 { 1.0 } else { -1.0 };
                    let v = gradient[i as i32];
                    let scale = g.scales()[5 + i];
                    writeln!(f,"{},{},{:.17e},{:.17e},{:.17e},{:.17e},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",polarization_name(pol),label,r[0],r[1],r[2],r[3],["ct","x","y","z"][i],h,s,exact[0],(s-exact[0])/g.peak,v,exact[i+1],fd,(v-exact[i+1])/scale,(v-fd)/scale)?;
                }
            }
        }
    }
    f.flush()?;
    report.finish()
}
