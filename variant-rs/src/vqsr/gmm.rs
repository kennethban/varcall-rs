//! GATK's VQSR Gaussian mixture model: `MultivariateGaussian`, `GaussianMixtureModel` and the variational Bayes EM
//! of `VariantRecalibratorEngine`, ported operation for operation. Every random draw comes from GATK's one
//! `Utils.getRandomGenerator()`, in GATK's order. The per-datum steps (k-means assignment, the E step, evaluation)
//! run in parallel: each datum's values depend only on that datum. The sums over data (the M step, the k-means means)
//! run in data order, as GATK's do, one Gaussian per thread.

use super::jama::Matrix;
use super::Datum;
use hc_assembly::jmath;
use hc_genotype::javamath::JavaRandom;
use rayon::prelude::*;

const EPSILON: f64 = 1e-200;
const COVARIANCE_REGULARIZATION_EPSILON: f64 = 1e-6;

/// `Math.log(10.0)`
fn ln10() -> f64 {
    jmath::log(10.0)
}

/// commons-math `Gamma.digamma`
fn digamma(x: f64) -> f64 {
    if x.is_nan() || x.is_infinite() {
        return x;
    }
    hc_assembly::m2math::digamma(x)
}

/// `MathUtils.maxElementIndex`: the first index of the maximum.
pub fn max_element_index(v: &[f64]) -> usize {
    let mut m = 0;
    for i in 1..v.len() {
        if v[i] > v[m] {
            m = i;
        }
    }
    m
}

/// `MathUtils.log10sumLog10` (a sum that is NaN or +Infinity is an IllegalArgumentException in GATK)
pub fn log10_sum_log10(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return if v.is_empty() { f64::NEG_INFINITY } else { v[0] };
    }
    let mi = max_element_index(v);
    let max = v[mi];
    if max == f64::NEG_INFINITY {
        return max;
    }
    let mut s = 0.0;
    for (i, &x) in v.iter().enumerate() {
        s += if i == mi { 0.0 } else { jmath::pow(10.0, x - max) };
    }
    let sum = 1.0 + s;
    assert!(!sum.is_nan() && sum != f64::INFINITY, "log10p values must be non-infinite and non-NAN");
    max + jmath::log10(sum)
}

/// `GaussianMixtureModel.nanTolerantLog10SumLog10`
fn nan_tolerant_log10_sum_log10(v: &[f64]) -> f64 {
    if v.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }
    log10_sum_log10(v)
}

/// `MathUtils.normalizeLog10DeleteMePlease(array, false)`: linear probabilities summing to one.
fn normalize_log10_linear(v: &[f64]) -> Vec<f64> {
    let max = v[max_element_index(v)];
    let mut n: Vec<f64> = v.iter().map(|&x| jmath::pow(10.0, x - max)).collect();
    let mut sum = 0.0;
    for &x in &n {
        sum += x;
    }
    for x in &mut n {
        *x /= sum;
    }
    n
}

/// `MathUtils.normalizeLog10DeleteMePlease(array, true)`, in place.
fn normalize_log10_in_place(v: &mut [f64]) {
    let max = v[max_element_index(v)];
    let mut sum = 0.0;
    for &x in v.iter() {
        sum += jmath::pow(10.0, x - max);
    }
    let log10_sum = jmath::log10(sum);
    for x in v.iter_mut() {
        *x = *x - max - log10_sum;
    }
}

/// `MathUtils.normalDistributionLog10` (GATK passes the variance as `sd`)
fn normal_distribution_log10(mean: f64, sd: f64, x: f64) -> f64 {
    assert!(sd >= 0.0, "sd: Standard deviation of normal must be > 0");
    assert!(mean.is_finite() && sd.is_finite() && x.is_finite(), "mean, sd, or, x : Normal parameters must be well formatted (non-INF, non-NAN)");
    let root_two_pi = (2.0 * std::f64::consts::PI).sqrt();
    let a = -1.0 * jmath::log10(sd * root_two_pi);
    let d = x - mean;
    let b = -1.0 * ((d * d) / (2.0 * (sd * sd))) / ln10();
    a + b
}

/// `MultivariateGaussian`
#[derive(Clone)]
pub struct Gaussian {
    pub p_mixture_log10: f64,
    pub sum_prob: f64,
    pub mu: Vec<f64>,
    pub sigma: Matrix,
    hyper_a: f64,
    hyper_b: f64,
    hyper_lambda: f64,
    cached_denom_log10: f64,
    cached_sigma_inverse: Matrix,
    /// probability of each training datum in this Gaussian (`pVarInGaussian`)
    p_var: Vec<f64>,
}

impl Gaussian {
    fn new(num_variants: usize, d: usize) -> Gaussian {
        Gaussian {
            p_mixture_log10: 0.0,
            sum_prob: 0.0,
            mu: vec![0.0; d],
            sigma: Matrix::zeros(d),
            hyper_a: 0.0,
            hyper_b: 0.0,
            hyper_lambda: 0.0,
            cached_denom_log10: 0.0,
            cached_sigma_inverse: Matrix::zeros(d),
            p_var: vec![0.0; num_variants],
        }
    }

    fn initialize_random_mu(&mut self, rng: &mut JavaRandom) {
        for m in &mut self.mu {
            *m = -4.0 + 8.0 * rng.next_double();
        }
    }

    fn initialize_random_sigma(&mut self, rng: &mut JavaRandom) {
        let d = self.mu.len();
        let mut r = Matrix::zeros(d);
        for i in 0..d {
            for j in i..d {
                let mut v = 0.55 + 1.25 * rng.next_double();
                if rng.next_boolean() {
                    v *= -1.0;
                }
                r.set(j, i, v);
                if i != j {
                    r.set(i, j, 0.0);
                }
            }
        }
        self.sigma = r.times(&r.transpose());
    }

    /// `MathUtils.distanceSquared(datum.annotations, mu)`
    fn distance_squared(&self, ann: &[f64]) -> f64 {
        let mut s = 0.0;
        for (x, m) in ann.iter().zip(&self.mu) {
            let d = x - m;
            s += d * d;
        }
        s
    }

    fn precompute_inverse(&mut self) -> Result<(), String> {
        self.cached_sigma_inverse = self.sigma.inverse().map_err(|e| {
            format!(
                "Error during clustering. Most likely there are too few variants used during Gaussian mixture modeling. Please consider raising the number of variants used to train the negative model (via --percentBadVariants 0.05, for example) or lowering the maximum number of Gaussians to use in the model (via --maxGaussians 4, for example). ({e})"
            )
        })?;
        Ok(())
    }

    /// `precomputeDenominatorForEvaluation`; an Err is GATK's exception, which marks the model as failed.
    fn precompute_denominator_for_evaluation(&mut self) -> Result<(), String> {
        if self.p_mixture_log10 == f64::NEG_INFINITY {
            return Ok(());
        }
        self.precompute_inverse()?;
        let d = self.mu.len() as f64;
        let det = self.sigma.det();
        self.cached_denom_log10 = jmath::log10(jmath::pow(2.0 * std::f64::consts::PI, (-1.0 * d) / 2.0)) + jmath::log10(jmath::pow(det, -0.5));
        if self.cached_denom_log10.is_nan() || det < EPSILON {
            return Err(format!("Denominator for gaussian evaluation cannot be computed. Covariance determinant is {det}."));
        }
        Ok(())
    }

    fn precompute_denominator_for_variational_bayes(&mut self, sum_hyper_lambda: f64) -> Result<(), String> {
        self.precompute_inverse()?;
        self.cached_sigma_inverse.times_equals(self.hyper_a);
        let d = self.mu.len();
        let mut sum = 0.0;
        for j in 1..=d {
            sum += digamma((self.hyper_a + 1.0 - j as f64) / 2.0);
        }
        sum -= jmath::log(self.sigma.det());
        sum += jmath::log(2.0) * d as f64;
        let lambda = 0.5 * sum;
        let pi = digamma(self.hyper_lambda) - digamma(sum_hyper_lambda);
        let beta = (-1.0 * d as f64) / (2.0 * self.hyper_b);
        self.cached_denom_log10 = (pi / ln10()) + (lambda / ln10()) + (beta / ln10());
        Ok(())
    }

    /// `evaluateDatumLog10`
    fn evaluate_log10(&self, ann: &[f64]) -> f64 {
        if self.p_mixture_log10 == f64::NEG_INFINITY {
            return f64::NEG_INFINITY;
        }
        let d = self.mu.len();
        let inv = &self.cached_sigma_inverse;
        let mut cross = [0.0f64; 16];
        for i in 0..d {
            for j in 0..d {
                cross[i] += (ann[j] - self.mu[j]) * inv.get(j, i);
            }
        }
        let mut sum_kernel = 0.0;
        for i in 0..d {
            sum_kernel += cross[i] * (ann[i] - self.mu[i]);
        }
        ((-0.5 * sum_kernel) / ln10()) + self.cached_denom_log10
    }

    /// sigma += sum over the data of prob * (x - mu)(x - mu)' + regularization, in data order
    fn add_weighted_scatter(&mut self, data: &[Datum], idx: &[usize]) {
        let d = self.mu.len();
        for (k, &di) in idx.iter().enumerate() {
            let prob = self.p_var[k];
            let a = &data[di].ann;
            for i in 0..d {
                for j in 0..d {
                    let reg = if i == j { COVARIANCE_REGULARIZATION_EPSILON } else { 0.0 };
                    let v = prob * (a[i] - self.mu[i]) * (a[j] - self.mu[j]) + reg;
                    self.sigma.a[i * d + j] += v;
                }
            }
        }
    }

    /// sum_prob and mu = the probability-weighted mean of the data, in data order
    fn weighted_mean(&mut self, data: &[Datum], idx: &[usize], start: f64) {
        self.sum_prob = start;
        self.mu.iter_mut().for_each(|m| *m = 0.0);
        for (k, &di) in idx.iter().enumerate() {
            let prob = self.p_var[k];
            self.sum_prob += prob;
            for (m, a) in self.mu.iter_mut().zip(&data[di].ann) {
                *m += prob * a;
            }
        }
        let s = self.sum_prob;
        self.mu.iter_mut().for_each(|m| *m /= s);
    }

    /// `maximizeGaussian`
    fn maximize(&mut self, data: &[Datum], idx: &[usize], emp_mu: &[f64], emp_sigma: &Matrix, shrinkage: f64, dirichlet: f64, dof: f64) {
        let d = self.mu.len();
        self.sigma = Matrix::zeros(d);
        self.weighted_mean(data, idx, 1e-10);
        let shrinkage_factor = (shrinkage * self.sum_prob) / (shrinkage + self.sum_prob);
        let mut wishart = Matrix::zeros(d);
        for i in 0..d {
            let delta_mu = shrinkage_factor * (self.mu[i] - emp_mu[i]);
            for j in 0..d {
                wishart.set(i, j, delta_mu * (self.mu[j] - emp_mu[j]));
            }
        }
        self.add_weighted_scatter(data, idx);
        self.sigma.plus_equals(emp_sigma);
        self.sigma.plus_equals(&wishart);
        for i in 0..d {
            self.mu[i] = (self.sum_prob * self.mu[i] + shrinkage * emp_mu[i]) / (self.sum_prob + shrinkage);
        }
        self.hyper_a = self.sum_prob + dof;
        self.hyper_b = self.sum_prob + shrinkage;
        self.hyper_lambda = self.sum_prob + dirichlet;
        self.p_var.iter_mut().for_each(|p| *p = 0.0);
    }

    /// `evaluateFinalModelParameters`
    fn evaluate_final(&mut self, data: &[Datum], idx: &[usize]) {
        let d = self.mu.len();
        self.sigma = Matrix::zeros(d);
        self.weighted_mean(data, idx, 0.0);
        self.add_weighted_scatter(data, idx);
        self.sigma.times_equals(1.0 / self.sum_prob);
        self.p_var.iter_mut().for_each(|p| *p = 0.0);
    }
}

/// `GaussianMixtureModel`
pub struct Gmm {
    pub gaussians: Vec<Gaussian>,
    shrinkage: f64,
    dirichlet: f64,
    prior_counts: f64,
    emp_mu: Vec<f64>,
    emp_sigma: Matrix,
    ready: bool,
    pub failed_to_converge: bool,
}

/// The engine's settings (`VariantRecalibratorArgumentCollection`).
#[derive(Clone)]
pub struct Settings {
    pub max_iterations: usize,
    pub kmeans_iterations: usize,
    pub shrinkage: f64,
    pub dirichlet: f64,
    pub prior_counts: f64,
}

/// `VariantRecalibratorEngine.MIN_ACCEPTABLE_LOD_SCORE`
const MIN_ACCEPTABLE_LOD_SCORE: f64 = -20000.0;
const MIN_PROB_CONVERGENCE: f64 = 2e-3;

impl Gmm {
    fn new(k: usize, n: usize, d: usize, s: &Settings) -> Gmm {
        Gmm {
            gaussians: (0..k).map(|_| Gaussian::new(n, d)).collect(),
            shrinkage: s.shrinkage,
            dirichlet: s.dirichlet,
            prior_counts: s.prior_counts,
            emp_mu: vec![0.0; d],
            emp_sigma: Matrix::identity(d).times_scalar(200.0).inverse().expect("200 * I is invertible"),
            ready: false,
            failed_to_converge: false,
        }
    }

    fn initialize_random_model(&mut self, data: &[Datum], idx: &[usize], kmeans_iterations: usize, rng: &mut JavaRandom) {
        for g in &mut self.gaussians {
            g.initialize_random_mu(rng);
        }
        self.kmeans(data, idx, kmeans_iterations, rng);
        let k = self.gaussians.len() as f64;
        for g in &mut self.gaussians {
            g.p_mixture_log10 = jmath::log10(1.0 / k);
            g.sum_prob = 1.0 / k;
            g.initialize_random_sigma(rng);
            g.hyper_a = self.prior_counts;
            g.hyper_b = self.shrinkage;
            g.hyper_lambda = self.dirichlet;
        }
    }

    /// `initializeMeansUsingKMeans`
    fn kmeans(&mut self, data: &[Datum], idx: &[usize], iterations: usize, rng: &mut JavaRandom) {
        for _ in 0..iterations {
            let gs = &self.gaussians;
            let assignment: Vec<usize> = idx
                .par_iter()
                .map(|&di| {
                    let mut min_distance = f64::MAX;
                    let mut min_g = usize::MAX;
                    for (gi, g) in gs.iter().enumerate() {
                        let dist = g.distance_squared(&data[di].ann);
                        if dist < min_distance {
                            min_distance = dist;
                            min_g = gi;
                        }
                    }
                    assert!(min_g != usize::MAX, "k-means: a datum is at no finite distance from any Gaussian (GATK fails here too)");
                    min_g
                })
                .collect();
            // the assigned data's sum per Gaussian, in data order
            let sums: Vec<(usize, Vec<f64>)> = (0..gs.len())
                .into_par_iter()
                .map(|gi| {
                    let mut mu = vec![0.0; gs[gi].mu.len()];
                    let mut n = 0usize;
                    for (k, &di) in idx.iter().enumerate() {
                        if assignment[k] == gi {
                            n += 1;
                            for (m, a) in mu.iter_mut().zip(&data[di].ann) {
                                *m += 1.0 * a;
                            }
                        }
                    }
                    (n, mu)
                })
                .collect();
            for (g, (n, mu)) in self.gaussians.iter_mut().zip(sums) {
                if n != 0 {
                    g.mu = mu;
                    g.mu.iter_mut().for_each(|m| *m /= n as f64);
                } else {
                    g.initialize_random_mu(rng);
                }
            }
        }
    }

    fn sum_hyper_lambda(&self) -> f64 {
        hc_stream_sum(self.gaussians.iter().map(|g| g.hyper_lambda))
    }

    /// `expectationStep`
    fn expectation_step(&mut self, data: &[Datum], idx: &[usize]) -> Result<(), String> {
        let sl = self.sum_hyper_lambda();
        for g in &mut self.gaussians {
            g.precompute_denominator_for_variational_bayes(sl)?;
        }
        let gs = &self.gaussians;
        let probs: Vec<Vec<f64>> = idx
            .par_iter()
            .map(|&di| {
                let p: Vec<f64> = gs.iter().map(|g| g.evaluate_log10(&data[di].ann)).collect();
                normalize_log10_linear(&p)
            })
            .collect();
        for (k, p) in probs.iter().enumerate() {
            for (g, &v) in self.gaussians.iter_mut().zip(p) {
                g.p_var[k] = v;
            }
        }
        Ok(())
    }

    /// `maximizationStep`: each Gaussian on its own thread
    fn maximization_step(&mut self, data: &[Datum], idx: &[usize]) {
        let (emp_mu, emp_sigma) = (&self.emp_mu, &self.emp_sigma);
        let (s, dir, pc) = (self.shrinkage, self.dirichlet, self.prior_counts);
        self.gaussians.par_iter_mut().for_each(|g| g.maximize(data, idx, emp_mu, emp_sigma, s, dir, pc));
    }

    fn evaluate_final_model_parameters(&mut self, data: &[Datum], idx: &[usize]) {
        self.gaussians.par_iter_mut().for_each(|g| g.evaluate_final(data, idx));
        self.normalize_p_mixture_log10();
    }

    /// `normalizePMixtureLog10`: the summed change of the mixture coefficients.
    fn normalize_p_mixture_log10(&mut self) -> f64 {
        let sum_pk = hc_stream_sum(self.gaussians.iter().map(|g| g.sum_prob));
        let log10_sum_pk = jmath::log10(sum_pk);
        let mut p: Vec<f64> = self.gaussians.iter().map(|g| jmath::log10(g.sum_prob) - log10_sum_pk).collect();
        normalize_log10_in_place(&mut p);
        let mut sum_diff = 0.0;
        for (g, &v) in self.gaussians.iter_mut().zip(&p) {
            sum_diff += (v - g.p_mixture_log10).abs();
            g.p_mixture_log10 = v;
        }
        sum_diff
    }

    /// `precomputeDenominatorForEvaluation`
    fn precompute_denominator_for_evaluation(&mut self) -> Result<(), String> {
        for g in &mut self.gaussians {
            g.precompute_denominator_for_evaluation()?;
        }
        self.ready = true;
        Ok(())
    }

    /// `evaluateDatum` for a datum without missing annotations
    fn evaluate_full(&self, ann: &[f64]) -> f64 {
        let p: Vec<f64> = self.gaussians.iter().map(|g| g.p_mixture_log10 + g.evaluate_log10(ann)).collect();
        nan_tolerant_log10_sum_log10(&p)
    }

    /// `evaluateDatumMarginalized`, given the standard-normal draws GATK makes for it (20 per missing annotation, in
    /// annotation order). The datum's missing annotations end at the last draw, as in GATK.
    fn evaluate_marginalized(&self, ann: &mut [f64], is_null: &[bool], draws: &[f64]) -> f64 {
        let mut sum = 0.0;
        let mut n = 0usize;
        for i in 0..ann.len() {
            if is_null[i] {
                for _ in 0..20 {
                    ann[i] = draws[n];
                    sum += jmath::pow(10.0, self.evaluate_full(ann));
                    n += 1;
                }
            }
        }
        jmath::log10(sum / n as f64)
    }

    /// `evaluateDatumInOneDimension`
    pub fn evaluate_in_one_dimension(&self, d: &Datum, i: usize) -> Option<f64> {
        if d.is_null[i] {
            return None;
        }
        let p: Vec<f64> = self
            .gaussians
            .iter()
            .map(|g| {
                let mut v = g.p_mixture_log10;
                if g.p_mixture_log10 != f64::NEG_INFINITY {
                    v += normal_distribution_log10(g.mu[i], g.sigma.get(i, i), d.ann[i]);
                }
                v
            })
            .collect();
        Some(nan_tolerant_log10_sum_log10(&p))
    }
}

/// `DoubleStream.sum()`
fn hc_stream_sum(v: impl IntoIterator<Item = f64>) -> f64 {
    crate::mutect::somatic::stream_sum(v)
}

/// `VariantRecalibratorEngine.generateModel`: k-means initialization, then variational Bayes EM.
pub fn generate_model(data: &[Datum], idx: &[usize], max_gaussians: usize, s: &Settings, rng: &mut JavaRandom) -> Result<Gmm, String> {
    if idx.is_empty() {
        return Err("No data found.".into());
    }
    let d = data[idx[0]].ann.len();
    let mut m = Gmm::new(max_gaussians, idx.len(), d, s);
    m.initialize_random_model(data, idx, s.kmeans_iterations, rng);
    m.normalize_p_mixture_log10();
    m.expectation_step(data, idx)?;
    let mut iteration = 0;
    while iteration < s.max_iterations {
        iteration += 1;
        m.maximization_step(data, idx);
        let change = m.normalize_p_mixture_log10();
        m.expectation_step(data, idx)?;
        if iteration > 2 && change < MIN_PROB_CONVERGENCE {
            break;
        }
    }
    eprintln!("[vqsr] model of {} Gaussians from {} variants: {iteration} iterations", max_gaussians, idx.len());
    m.evaluate_final_model_parameters(data, idx);
    Ok(m)
}

/// `VariantRecalibratorEngine.evaluateData` over all data. Positive model (`contrastive` false): lod = the model's
/// log10 likelihood. Negative model: lod = prior + positive lod - negative lod, and a positive lod of -Infinity becomes
/// `MIN_ACCEPTABLE_LOD_SCORE * (1 + U)`. Random draws (marginalization, that replacement) are made first, sequentially
/// in data order as GATK makes them; the evaluation itself is parallel.
pub fn evaluate_data(data: &mut [Datum], model: &mut Gmm, contrastive: bool, rng: &mut JavaRandom) {
    if !model.ready {
        if let Err(e) = model.precompute_denominator_for_evaluation() {
            eprintln!("[vqsr] WARN Model could not pre-compute denominators. {e}");
            model.failed_to_converge = true;
            return;
        }
    }
    // draws, per datum in order: 20 normals per missing annotation, then (contrastive, positive lod infinite) one
    // uniform for the replacement
    let mut draws: Vec<(Vec<f64>, f64)> = Vec::with_capacity(data.len());
    for d in data.iter() {
        let nulls = d.is_null.iter().filter(|&&x| x).count();
        let g: Vec<f64> = (0..20 * nulls).map(|_| rng.next_gaussian()).collect();
        let u = if contrastive && d.lod.is_infinite() { rng.next_double() } else { 0.0 };
        draws.push((g, u));
    }
    let model_ref = &*model;
    let this_lods: Vec<f64> = data
        .par_iter_mut()
        .zip(draws.par_iter())
        .map(|(d, (g, _))| {
            if g.is_empty() {
                model_ref.evaluate_full(&d.ann)
            } else {
                let is_null = d.is_null.clone();
                model_ref.evaluate_marginalized(&mut d.ann, &is_null, g)
            }
        })
        .collect();
    for (k, d) in data.iter_mut().enumerate() {
        let this_lod = this_lods[k];
        if this_lod.is_nan() {
            // GATK stops here: this and later data keep their lods; the run then fails
            model.failed_to_converge = true;
            return;
        }
        d.lod = if contrastive {
            if d.lod.is_infinite() {
                MIN_ACCEPTABLE_LOD_SCORE + draws[k].1 * MIN_ACCEPTABLE_LOD_SCORE
            } else {
                d.prior + d.lod - this_lod
            }
        } else {
            this_lod
        };
    }
}
