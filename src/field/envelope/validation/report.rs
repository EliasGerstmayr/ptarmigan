use super::{
    gaussian::{polarization_name, Gaussian},
    runner::{method_name, Config, Query},
    statistics::{rate, Errors, COMPONENTS},
    Result,
};
use crate::field::envelope::InterpolationMethod;
use std::{
    fs::{File, OpenOptions},
    io::{BufWriter, Write},
    path::PathBuf,
    process::Command,
};

pub struct Report {
    directory: PathBuf,
    metadata: BufWriter<File>,
}
fn git(args: &[&str]) -> String {
    Command::new("git")
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned())
        .unwrap_or_else(|| "unavailable".into())
}
fn optional(x: Option<f64>) -> String {
    x.map(|v| format!("{:.17e}", v))
        .unwrap_or_else(|| "NA".into())
}

impl Report {
    pub fn new(c: &Config, mode: &str) -> Result<Self> {
        // Runtime provenance describes the working tree immediately before report
        // creation, including any pre-existing untracked files in output directories.
        let commit = git(&["rev-parse", "HEAD"]);
        let status = git(&["status", "--porcelain"]);
        std::fs::create_dir_all(&c.output)?;
        let metadata = BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(c.output.join("metadata.txt"))?,
        );
        let mut r = Self {
            directory: c.output.clone(),
            metadata,
        };
        writeln!(r.metadata,"format=ptarmigan-envelope-validation-v2\nrun_status=started\ndiagnostic={}\nnative_diagnostic_run={}\ngaussian_study_run={}",mode,mode=="native_gaussian_diagnostic",mode=="gaussian_convergence")?;
        writeln!(
            r.metadata,
            "git_commit={}\ngit_dirty={}\ngit_status={:?}",
            commit,
            if status == "unavailable" {
                "unavailable"
            } else if status.is_empty() {
                "false"
            } else {
                "true"
            },
            status
        )?;
        writeln!(r.metadata,"compiled_git={}\npackage_version={}\nfeatures={}\ndebug_assertions={}\narchitecture={}\noperating_system={}",option_env!("VERGEN_GIT_SHA").unwrap_or("unavailable"),env!("CARGO_PKG_VERSION"),option_env!("PTARMIGAN_ACTIVE_FEATURES").unwrap_or("unavailable"),cfg!(debug_assertions),std::env::consts::ARCH,std::env::consts::OS)?;
        writeln!(r.metadata,"parameters={:?}\nunits=metres;duration_seconds;dimensionless_a0\npolarizations=linear,circular\naxis_order=x,y,z,xi\nxi=ct-z\nxi_c_convention=retarded_offset;focus_peak_ct=z_f+xi_c",c.params)?;
        writeln!(r.metadata,"configured_cells={:?}\nconfigured_nodes={:?}\ndomain_half_extent_in_characteristic_lengths={:.17e}\nscalar_memory_limit_bytes={}\nmemory_kind=estimated_scalar_allocation_not_measured_peak_process_memory",c.cells,c.cells.iter().map(|n| n.checked_add(1)).collect::<Vec<_>>(),c.extent,c.limit)?;
        writeln!(r.metadata,"sampling_seed={}\nconfigured_sample_count={}\nsampling_algorithm=Xoshiro256StarStar_seed_from_u64_rand0.7;uniform_95percent_domain;reject_symmetry_1e-6_and_knots_1e-7cell;lab_roundtrip_xi\nquery_reuse=same_physical_points_all_resolutions_methods_and_polarizations",c.seed,c.samples)?;
        writeln!(r.metadata,"normalization=S_peak;S_peak/w0;S_peak/w0;S_peak/zR;S_peak/L;S_peak/L;S_peak/w0;S_peak/w0;S_peak/zR+S_peak/L\nrate=ln(error_coarse/error_fine)/ln(cells_fine/cells_coarse);NA_for_zero_nonfinite_or_no_previous\nexpected_multilinear_orders=value_2_gradient_1_not_asserted\nprovisional_sampled_targets=max_normalized_value<1e-3;max_normalized_gradient<1e-2;not_physics_tolerances\nauto_refinement=false\nparticles_run=false\nperformance_measured=false")?;
        if mode == "native_gaussian_diagnostic" {
            writeln!(r.metadata,"sampling_override=two_fixed_points;configured_random_sampling_and_grids_not_used\nfd_relative_steps=1e-3,1e-4,1e-5")?;
        }
        if mode == "gaussian_convergence" {
            if c.query_file.is_some() {
                writeln!(r.metadata,"sampling_override=loaded_file;configured_count_seed_and_generation_algorithm_not_used")?;
            }
            writeln!(r.metadata, "interpolation_methods=multilinear,cubic\nquery_source={:?}\nloaded_query_policy=exact_f64_xi_equals_ct_minus_z;strict_canonical_bounds_all_resolutions;no_modification\nsample_scope=sampled_points_only;no_full_domain_positivity_guarantee\nexpected_cubic_orders=value_3_gradient_2_not_asserted", c.query_file)?;
        }
        r.metadata.flush()?;
        println!("{} -> {}", mode, c.output.display());
        Ok(r)
    }
    pub fn file(&self, name: &str) -> Result<BufWriter<File>> {
        Ok(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(self.directory.join(name))?,
        ))
    }
    pub fn gaussian(&mut self, g: &Gaussian) -> Result<()> {
        let bounds: [[f64; 2]; 4] = std::array::from_fn(|i| [g.params.offsets[i], g.lengths()[i]]);
        writeln!(self.metadata,"{}: zR_m={:.17e} L_m={:.17e} S_peak={:.17e} normalization_scales={:?} domain_centres_and_characteristic_lengths={:?}",polarization_name(g.pol),g.zr,g.length,g.peak,g.scales(),bounds)?;
        Ok(())
    }
    pub fn queries(&mut self, points: &[Query]) -> Result<()> {
        writeln!(self.metadata, "actual_sample_count={}", points.len())?;
        let mut f = self.file("queries.csv")?;
        writeln!(f, "x_m,y_m,z_m,xi_m,ct_m")?;
        for q in points {
            writeln!(
                f,
                "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                q.stored[0], q.stored[1], q.stored[2], q.stored[3], q.lab[0]
            )?;
        }
        f.flush()?;
        Ok(())
    }
    pub fn axes(
        &self,
        f: &mut impl Write,
        g: &Gaussian,
        grid: &crate::field::envelope::EnvelopeGrid,
        cells: usize,
        bytes: usize,
    ) -> Result<()> {
        for (i, a) in grid.axes().iter().enumerate() {
            writeln!(
                f,
                "{},{},{},{:.17e},{:.17e},{},{:.17e},{}",
                polarization_name(g.pol),
                cells,
                ["x", "y", "z", "xi"][i],
                a.origin(),
                a.spacing(),
                a.len(),
                a.coordinate(a.len() - 1).unwrap(),
                bytes
            )?;
        }
        Ok(())
    }
    pub fn error_header(&self, f: &mut impl Write) -> Result<()> {
        writeln!(f,"method,polarization,cells,nodes,estimated_scalar_bytes,component,scale,sampled_max_abs,rms_abs,sampled_max_normalized,rms_normalized,refinement_ratio,sampled_max_rate,rms_rate,provisional_sampled_target,sampled_target_met")?;
        Ok(())
    }
    pub fn errors(
        &self,
        f: &mut impl Write,
        g: &Gaussian,
        method: InterpolationMethod,
        cells: usize,
        bytes: usize,
        stats: &[Errors; 9],
        previous: Option<&(usize, [Errors; 9])>,
    ) -> Result<()> {
        for i in 0..9 {
            let scale = g.scales()[i];
            let (max, rms) = stats[i].normalized(scale)?;
            let ratio = previous.map(|p| cells as f64 / p.0 as f64);
            let max_rate = previous.and_then(|p| rate(p.1[i].max, stats[i].max, ratio.unwrap()));
            let rms_rate =
                previous.and_then(|p| rate(p.1[i].rms(), stats[i].rms(), ratio.unwrap()));
            let target = if i == 0 { 1e-3 } else { 1e-2 };
            writeln!(
                f,
                "{},{},{},{},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{},{},{:.17e},{}",
                method_name(method),
                polarization_name(g.pol),
                cells,
                cells + 1,
                bytes,
                COMPONENTS[i],
                scale,
                stats[i].max,
                stats[i].rms(),
                max,
                rms,
                optional(ratio),
                optional(max_rate),
                optional(rms_rate),
                target,
                max < target
            )?;
            println!(
                "{} {} n={} {:>9}: sampled max={:.3e} rms={:.3e} p(max)={} sampled target={}",
                method_name(method),
                polarization_name(g.pol),
                cells,
                COMPONENTS[i],
                max,
                rms,
                optional(max_rate),
                if max < target { "met" } else { "MISSED" }
            );
        }
        Ok(())
    }
    pub fn finish(mut self) -> Result<()> {
        writeln!(self.metadata, "run_status=completed")?;
        self.metadata.flush()?;
        Ok(())
    }
}
