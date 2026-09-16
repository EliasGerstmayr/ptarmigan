use super::Result;

pub const COMPONENTS: [&str; 9] = [
    "S",
    "stored_x",
    "stored_y",
    "stored_z",
    "stored_xi",
    "lab_ct",
    "lab_x",
    "lab_y",
    "lab_z",
];

#[derive(Clone, Copy, Default)]
pub struct Errors {
    pub max: f64,
    norm: f64,
    count: usize,
}
impl Errors {
    pub fn add(&mut self, actual: f64, reference: f64) -> Result<()> {
        let error = (actual - reference).abs();
        if !error.is_finite() {
            return Err("nonfinite error measurement".into());
        }
        self.max = self.max.max(error);
        self.norm = self.norm.hypot(error);
        self.count += 1;
        if !self.norm.is_finite() {
            return Err("error norm overflow".into());
        }
        Ok(())
    }
    pub fn rms(&self) -> f64 {
        if self.count == 0 {
            f64::NAN
        } else {
            self.norm / (self.count as f64).sqrt()
        }
    }
    pub fn normalized(&self, scale: f64) -> Result<(f64, f64)> {
        if !scale.is_finite() || scale <= 0.0 {
            return Err("invalid error normalization scale".into());
        }
        let values = (self.max / scale, self.rms() / scale);
        if !values.0.is_finite() || !values.1.is_finite() {
            return Err("nonfinite normalized errors or empty samples".into());
        }
        Ok(values)
    }
}

pub fn rate(coarse: f64, fine: f64, ratio: f64) -> Option<f64> {
    if !coarse.is_finite()
        || !fine.is_finite()
        || coarse <= 0.0
        || fine <= 0.0
        || !ratio.is_finite()
        || ratio <= 1.0
    {
        return None;
    }
    let p = (coarse.ln() - fine.ln()) / ratio.ln();
    if p.is_finite() {
        Some(p)
    } else {
        None
    }
}

#[test]
fn normalization_and_rates() {
    let mut e = Errors::default();
    e.add(2.0, 0.0).unwrap();
    e.add(4.0, 0.0).unwrap();
    let (max, rms) = e.normalized(2.0).unwrap();
    assert_eq!(max, 2.0);
    assert!((rms - 2.5f64.sqrt()).abs() < 4.0 * f64::EPSILON);
    assert!((rate(0.09, 0.01, 3.0).unwrap() - 2.0).abs() < 8.0 * f64::EPSILON);
    assert_eq!(rate(0.0, 0.0, 2.0), None);
    assert_eq!(rate(1.0, 0.0, 2.0), None);
    assert_eq!(rate(f64::INFINITY, 1.0, 2.0), None);
    assert_eq!(rate(1.0, f64::NAN, 2.0), None);
    assert_eq!(rate(1.0, 1.0, 1.0), None);
    assert!(e.add(f64::NAN, 0.0).is_err());
    assert!(e.normalized(0.0).is_err());
    assert!(Errors::default().normalized(1.0).is_err());
}
