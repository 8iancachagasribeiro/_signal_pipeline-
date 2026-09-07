# BRM analysis specification

## Purpose

This document converts the strategic revision request into executable analyses without fabricating results. It should be used with the manuscript's original simulation code so that the exact estradiol trajectory, estimators, preprocessing, and calibration are preserved.

## A. Monte Carlo architecture

Use at least 1,000 replications per cell. Prefer 2,000 for cells used to estimate Type I error or 95% interval coverage when computationally feasible. For every cell save:

- seed;
- generative parameters;
- sampling schedule;
- missingness mechanism;
- estimator status/failure flag;
- true individual parameters;
- estimated individual parameters;
- interval limits if available;
- pooled effect;
- heterogeneity statistic and p value;
- SSF and comparator outputs.

Report Monte Carlo standard errors:
- For a proportion p over R replications: MCSE = sqrt[p(1-p)/R].
- For a mean metric: MCSE = SD(metric across replications)/sqrt(R).
- For quantiles, use simulation-level bootstrap intervals or a sufficiently large R.

## B. Core condition grid

### Sampling density
n = 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 48, 56 observations/person.

### Number of cycles
1, 2, 3, 4 cycles.

Do not hold only total n constant. Cross total n with cycle count because identifiability can differ when the same number of observations is spread across more cycles.

### Directional heterogeneity
Evaluate at least:
1. balanced unimodal;
2. 60/40 positive-negative;
3. 70/30;
4. 80/20;
5. shifted normal sensitivity distribution;
6. symmetric bimodal;
7. asymmetric bimodal.

The manuscript should distinguish changes in the *mean* of the individual-effect distribution from changes in its *shape*.

### Measurement quality
Use granular SNR/error levels, including:
- no added measurement error;
- low, moderate, high white noise;
- AR(1) measurement error with at least rho = .3, .6;
- if predictor and outcome error are both simulated, vary them independently.

## C. Model violations

### Non-stationarity
- linear drift in baseline outcome;
- slow random-walk or low-frequency drift;
- optional person-specific drift.

### Residual serial dependence
AR(1) residuals with rho = 0, .2, .4, .6, .8.

### Heteroscedasticity
Variance depends on cycle phase or E2 level. Use at least:
- constant variance;
- 1.5x variance near periovulatory peak;
- 2x variance near peak.

### Functional form
- inverted-U reference;
- sigmoid;
- asymmetric linear;
- piecewise linear;
- optional saturating Hill/logistic response.

### Hormonal gain
Replace fixed K with K_i:
- small CV;
- moderate CV;
- large CV;
- optionally correlate K_i with baseline b_i to test confounding of mechanisms.

### Outliers
- 1%, 3%, 5% contaminated observations;
- additive extreme values and leverage-like predictor contamination.

### Missingness
- MCAR at 5%, 10%, 20%;
- weekend/systematic missingness;
- phase-dependent missingness;
- outcome-dependent logistic missingness approximating MNAR.

## D. Recoverability metrics

For each cell report:
1. F = cor(beta_hat, beta_true);
2. RMSE;
3. mean bias;
4. median absolute error;
5. 95% interval coverage, if intervals exist;
6. sign accuracy;
7. rank-order reliability/ICC across repeated realizations.

Do not define a design as "adequate" from F alone. A recommended design region should satisfy:
- F >= .70;
- low bias relative to the true slope SD;
- acceptable RMSE;
- sign accuracy reported explicitly;
- near-nominal interval coverage when intervals are used.

## E. SSF validation study

### Truth definition
Define the target quantity before comparing methods. A defensible target is:
true smooth fraction = Var(s_t) / Var(s_t + epsilon_t)
for a known latent smooth signal s_t and specified stochastic high-frequency component.

### Signal families
- sinusoid;
- two-peak cycle;
- narrow pulse;
- asymmetric cycle;
- smooth trend + periodic component;
- smooth signal plus genuine high-frequency biological component.

### Series lengths
Include lengths corresponding to 1-4 cycles and the observed empirical range.

### SSF cutoff sensitivity
Evaluate fixed cutoffs around the current .25 threshold, for example:
.15, .20, .25, .30, .35 cycles/unit,
plus an adaptive plateau rule.

### Comparator methods
1. state-space/Kalman smoothing;
2. wavelet denoising;
3. penalized spline/GAM smoothing with CV/GCV;
4. autocorrelation-based structured-variance estimator.

### Outcomes
- bias;
- MSE;
- median absolute error;
- Monte Carlo SD;
- estimator failure rate;
- computation time as a secondary practical metric.

The conclusion should be conditional. "SSF is useful when..." is preferable to "SSF is superior."

## F. Sampling-design study

Compare:
- uniform;
- phase-targeted;
- derivative-adaptive;
- peak-targeted with incorrect assumed ovulation day.

Peak uncertainty:
shift the true peak relative to the planned schedule by 0, ±1, ±2, ±3, ±4 days.

Missing schedule:
- random missing days;
- weekend non-response;
- clustered missed blocks.

Output:
for each design, identify the smallest n/cycles region reaching F >= .70 and report the accompanying RMSE and sign accuracy.

## G. Surrogate-test calibration and heterogeneity power

### Type I error
Generate no true between-person slope variance while preserving:
- realistic autocorrelation;
- observed-like sampling gaps;
- series length distribution.

Estimate empirical alpha at nominal .05.

### Power
Cross:
- N participants;
- observations/person;
- true slope SD;
- residual rho;
- measurement error.

Report 50%, 80%, and 90% power contours. Keep this separate from person-specific recovery.

## H. Recommended manuscript figures

1. **Figure 1:** SSF benchmark: bias/MSE by estimator and series length.
2. **Figure 2:** Cancellation and directional-distribution boundary conditions.
3. **Figure 3:** Recovery response surfaces: n × cycles × measurement quality.
4. **Figure 4:** Robustness heat map across model violations.
5. **Figure 5:** Sampling design under peak uncertainty/missingness.
6. **Figure 6:** Surrogate Type I error and power contours.
7. **Figure 7:** Empirical mcPHASES illustration.

Move secondary grids to supplement rather than compressing all conditions into the main text.

## I. BRM-specific editorial requirements

Behavior Research Methods currently states that:
- the abstract should be no more than 250 words;
- manuscripts should provide variability/effect-size information and statistical power when relevant;
- an Open Practices Statement must appear immediately before the References;
- Level 2 transparency means relevant validation materials, analyses, and code should be placed in a trusted repository.

The revised manuscript should use these requirements as structural constraints, not merely mention reproducibility in prose.

## J. Claims that must not be made until the new analyses exist

Do not state that:
- SSF is unbiased;
- SSF outperforms Kalman/wavelet/spline methods;
- the .25 cutoff is optimal;
- the surrogate test is conservative or correctly calibrated;
- 20 observations is robust across violations;
- fidelity >= .70 is reached under MNAR, AR(1), non-stationarity, or alternative functional forms;
- interval coverage is nominal;
- a particular number of participants gives 80% power for heterogeneity.

Those are empirical simulation results and must be calculated.
