use super::Result;
use crate::{constants::SPEED_OF_LIGHT, field::Polarization};
use std::f64::consts::{LN_2, PI};

#[derive(Clone, Copy, Debug)]
pub struct Parameters {
    pub wavelength: f64,
    pub waist: f64,
    pub duration: f64,
    pub a0: f64,
    /// x_c, y_c, z_f, xi_c; all metres. xi_c is a retarded-coordinate offset.
    pub offsets: [f64; 4],
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            wavelength: 0.8e-6,
            waist: 5e-6,
            duration: 30e-15,
            a0: 1.0,
            offsets: [0.0; 4],
        }
    }
}

pub struct Gaussian {
    pub params: Parameters,
    pub pol: Polarization,
    pub zr: f64,
    pub length: f64,
    pub peak: f64,
}

impl Gaussian {
    pub fn new(params: Parameters, pol: Polarization) -> Result<Self> {
        for (name, value) in [
            ("wavelength", params.wavelength),
            ("waist", params.waist),
            ("duration", params.duration),
            ("a0", params.a0),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(format!(
                    "{} must be finite and positive for normalized validation",
                    name
                )
                .into());
            }
        }
        if params.offsets.iter().any(|x| !x.is_finite()) {
            return Err("offsets must be finite metres".into());
        }
        let zr = PI * params.waist.powi(2) / params.wavelength;
        let length = SPEED_OF_LIGHT * params.duration;
        let peak = match pol {
            Polarization::Linear => 0.5,
            Polarization::Circular => 1.0,
        } * params.a0.powi(2);
        let result = Self {
            params,
            pol,
            zr,
            length,
            peak,
        };
        // Reject unrepresentable derived scales before evaluating exponent/gradients.
        for value in [params.waist.powi(2), zr.powi(2), length.powi(2), peak]
            .iter()
            .copied()
            .chain(result.scales().iter().copied())
        {
            if !value.is_finite() || value <= 0.0 || !value.recip().is_finite() {
                return Err("Gaussian derived scales are not representable".into());
            }
        }
        Ok(result)
    }

    pub fn lengths(&self) -> [f64; 4] {
        [self.params.waist, self.params.waist, self.zr, self.length]
    }
    pub fn scales(&self) -> [f64; 9] {
        let [w, _, z, l] = self.lengths();
        [
            self.peak,
            self.peak / w,
            self.peak / w,
            self.peak / z,
            self.peak / l,
            self.peak / l,
            self.peak / w,
            self.peak / w,
            self.peak / z + self.peak / l,
        ]
    }

    pub fn scalar(&self, q: [f64; 4]) -> f64 {
        let [x, y, z, xi] = std::array::from_fn(|i| q[i] - self.params.offsets[i]);
        let w = 1.0 + (z / self.zr).powi(2);
        self.peak / w
            * (-2.0 * (x * x + y * y) / (self.params.waist.powi(2) * w)
                - 4.0 * LN_2 * (xi / self.length).powi(2))
            .exp()
    }

    /// Independent reference: S followed by four stored derivatives. No sampler calls.
    pub fn stored(&self, q: [f64; 4]) -> [f64; 5] {
        let [x, y, z, xi] = std::array::from_fn(|i| q[i] - self.params.offsets[i]);
        let w = 1.0 + (z / self.zr).powi(2);
        let transverse = self.params.waist.powi(2) * w;
        let s = self.scalar(q);
        [
            s,
            -4.0 * x * s / transverse,
            -4.0 * y * s / transverse,
            2.0 * z / (self.zr.powi(2) * w) * (2.0 * (x * x + y * y) / transverse - 1.0) * s,
            -8.0 * LN_2 * xi * s / self.length.powi(2),
        ]
    }

    /// Independent chain rule, not the production laboratory conversion helper.
    pub fn lab(&self, r: [f64; 4]) -> [f64; 5] {
        let [s, sx, sy, sz, sxi] = self.stored([r[1], r[2], r[3], r[0] - r[3]]);
        [s, sxi, -sx, -sy, sxi - sz]
    }
    pub fn scalar_lab(&self, r: [f64; 4]) -> f64 {
        self.scalar([r[1], r[2], r[3], r[0] - r[3]])
    }
}

pub fn polarization_name(pol: Polarization) -> &'static str {
    match pol {
        Polarization::Linear => "linear",
        Polarization::Circular => "circular",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn reference() -> Gaussian {
        Gaussian::new(Parameters::default(), Polarization::Linear).unwrap()
    }

    #[test]
    fn peak_fwhm_and_transverse_radius() {
        for pol in [Polarization::Linear, Polarization::Circular] {
            let g = Gaussian::new(
                Parameters {
                    a0: 2.0,
                    ..Parameters::default()
                },
                pol,
            )
            .unwrap();
            let expected = if pol == Polarization::Linear {
                2.0
            } else {
                4.0
            };
            assert_eq!(g.scalar([0.0; 4]), expected);
            for sign in [-1.0, 1.0] {
                assert!(
                    (g.scalar([0.0, 0.0, 0.0, sign * g.length / 2.0]) / expected - 0.5).abs()
                        < 8.0 * f64::EPSILON
                );
            }
            assert!(
                (g.scalar([g.params.waist, 0.0, 0.0, 0.0]) / expected - (-2.0f64).exp()).abs()
                    < 8.0 * f64::EPSILON
            );
        }
    }

    #[test]
    fn parameter_validation() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            for field in 0..4 {
                let mut p = Parameters::default();
                match field {
                    0 => p.wavelength = bad,
                    1 => p.waist = bad,
                    2 => p.duration = bad,
                    _ => p.a0 = bad,
                }
                assert!(Gaussian::new(p, Polarization::Linear).is_err());
            }
        }
        let mut p = Parameters::default();
        p.offsets[2] = f64::NAN;
        assert!(Gaussian::new(p, Polarization::Linear).is_err());
        p = Parameters::default();
        p.waist = f64::MAX;
        assert!(Gaussian::new(p, Polarization::Linear).is_err());
    }

    #[test]
    fn independent_scalar_finite_differences() {
        let mut p = Parameters::default();
        p.offsets = [1e-6, -2e-6, 3e-6, -1e-6];
        let g = Gaussian::new(p, Polarization::Circular).unwrap();
        for factors in [[0.31, -0.27, 0.43, 0.22], [-0.42, 0.19, -0.37, -0.31]] {
            let q: [f64; 4] = std::array::from_fn(|i| p.offsets[i] + factors[i] * g.lengths()[i]);
            let r = [q[2] + q[3], q[0], q[1], q[2]];
            let expected = g.lab(r);
            let hs = [g.length, p.waist, p.waist, g.length.min(g.zr)];
            for epsilon in [1e-3, 1e-4, 1e-5] {
                for i in 0..4 {
                    let h = epsilon * hs[i];
                    let mut a = r;
                    let mut b = r;
                    a[i] += h;
                    b[i] -= h;
                    let fd = (g.scalar_lab(a) - g.scalar_lab(b)) / (2.0 * h)
                        * if i == 0 { 1.0 } else { -1.0 };
                    // O(h^2) truncation plus O(eps_machine/h) roundoff, scaled.
                    let tol = 20.0 * epsilon * epsilon + 128.0 * f64::EPSILON / epsilon;
                    assert!((fd - expected[i + 1]).abs() / g.scales()[5 + i] < tol);
                    // Also independently differentiate each stored coordinate.
                    let h = epsilon * g.lengths()[i];
                    let mut a = q;
                    let mut b = q;
                    a[i] += h;
                    b[i] -= h;
                    let fd = (g.scalar(a) - g.scalar(b)) / (2.0 * h);
                    assert!((fd - g.stored(q)[i + 1]).abs() / g.scales()[i + 1] < tol);
                }
            }
        }
    }

    #[test]
    fn symmetries_translations_and_delay_convention() {
        let g = reference();
        let q = [
            0.3 * g.params.waist,
            -0.2 * g.params.waist,
            0.4 * g.zr,
            0.1 * g.length,
        ];
        for i in 0..4 {
            let mut reflected = q;
            reflected[i] = -reflected[i];
            assert_eq!(g.scalar(q), g.scalar(reflected));
            assert_eq!(g.stored(q)[i + 1], -g.stored(reflected)[i + 1]);
            let mut plane = q;
            plane[i] = 0.0;
            assert_eq!(g.stored(plane)[i + 1], 0.0);
        }
        let offsets = [1e-6, -2e-6, 3e-6, -4e-6];
        let translated = Gaussian::new(
            Parameters {
                offsets,
                ..g.params
            },
            g.pol,
        )
        .unwrap();
        let moved = std::array::from_fn(|i| q[i] + offsets[i]);
        for i in 0..5 {
            assert!(
                (translated.stored(moved)[i] - g.stored(q)[i]).abs() / g.scales()[i]
                    < 32.0 * f64::EPSILON
            );
        }
        let peak_r = [offsets[2] + offsets[3], offsets[0], offsets[1], offsets[2]];
        assert!((translated.scalar_lab(peak_r) / g.peak - 1.0).abs() < 8.0 * f64::EPSILON);
        // z_f alone shifts the time at which the pulse peak reaches focus.
        let shifted = Gaussian::new(
            Parameters {
                offsets: [0.0, 0.0, g.length, 0.0],
                ..g.params
            },
            g.pol,
        )
        .unwrap();
        assert!(shifted.scalar_lab([0.0, 0.0, 0.0, g.length]) < g.peak);
    }
}
