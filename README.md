# nidesigns

`nidesigns` implements the general framework for designing and evaluating
active-controlled trials with non-inferiority objectives described by
Olivas-Martinez, Gao, and Janes (2026). It supports preservation-of-effect and
inferred-efficacy criteria, fixed-margin and synthesis methods, conditional and
unconditional power, and sensitivity analyses for non-constancy.

Within the parameter region considered in Section 5 of the paper, conditional
operating characteristics are more optimistic: conditional power is higher
than unconditional power, while conditional type-I error is lower than
unconditional type-I error.

The package preserves the `ni_design()` interface and the numerical results in
Appendix B of the paper. It also extends sample-size planning to relative risks,
risk differences, and mean differences.

> Olivas-Martinez A, Gao F, Janes H. A General Framework for Designing and
> Evaluating Active-Controlled Trials with Non-Inferiority Objectives.
> *Statistics in Medicine*. 2026;45:e70618.
> [doi:10.1002/sim.70618](https://doi.org/10.1002/sim.70618)

## Background

Active-controlled non-inferiority trials are used when an effective
intervention already exists and a placebo-controlled design is not feasible or
ethical. The objective is to show that an experimental intervention is not
unacceptably worse than the active control while it may offer advantages in
safety, tolerability, cost, or feasibility.

Because these trials do not include a placebo arm, assessing non-inferiority
relies on historical estimates of the active-control effect relative to placebo
and assumptions about how well that evidence applies to the target trial. The
traditional preservation-of-effect criterion requires the experimental
intervention to preserve a specified fraction of that effect. When the active
control is highly effective, an inferred-efficacy criterion can instead require
the experimental intervention to exceed a prespecified efficacy threshold
relative to a hypothetical placebo.

The framework makes the historical assumptions explicit, represents common
fixed-margin and synthesis methods through `(u, lambda1)`, and evaluates both
conditional and unconditional operating characteristics.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("aolivasm/ni-designs")
```

Load the package with:

```r
library(nidesigns)
```

## Quick start

The example below compares the traditional synthesis and 95-95 methods for a
time-to-event outcome. Effects are analyzed on the log-hazard-ratio scale, with
hazard ratios below one indicating benefit.

```r
hist_hr <- 0.45
hist_se <- 0.25

design <- ni_design_hr(
  u.list = c(1, 0),
  l1.list = c(0, 1.96 * hist_se / log(hist_hr)),
  f.preserv = 0.5,
  null.hr = NULL,
  design.alternative.hr = 0.30,
  hist.ac.hr = hist_hr,
  hist.ac.effect.se = hist_se,
  target.on.unconditional.power = TRUE,
  power = 0.90
)

summary(design)
```

The principal method settings are:

| Method | `u` | `l1` |
|---|---:|---:|
| Traditional synthesis | `1` | `0` |
| Bias-adjusted synthesis | `1` | prespecified relative effect deviation |
| Odem-Davis | `1 / (1 + l1)` | prespecified relative effect deviation |
| 95-95 fixed margin | `0` | `1.96 * hist.ac.effect.se / log(hist.ac.hr)` |
| 0-95 fixed margin | `0` | `0` |

The primary user-facing functions are:

- `ni_design()` for prevention efficacy and log hazard-ratio designs.
- `ni_design_hr()` when hazard ratios are the natural inputs.
- `ni_design_rr()`, `ni_design_rd()`, and `ni_design_md()` for relative risks,
  risk differences, and mean differences.
- `explore_max_uncond_power()` for examining the maximum attainable
  unconditional power across design alternatives.

Use `?ni_design` and the corresponding function help pages for complete input
and output definitions.

For every requested method and success criterion, the design output reports the
non-inferiority margin, required number of events or target variance, sample
size, controlled non-constancy value, and unconditional power. Set
`target.on.unconditional.power = FALSE` to reproduce conventional planning that
targets conditional power instead.

## Constructing the historical input

The package expects a historical active-control-versus-placebo estimate and its
standard error on the analysis scale. It does not select historical studies or
fit a meta-analysis. Those steps must be completed before calling a design
function.

For independent study estimates `gamma[j]` with standard errors `se[j]`, a
conventional fixed-effect inverse-variance summary is

```r
weight <- 1 / se^2
hist_effect <- sum(weight * gamma) / sum(weight)
hist_effect_se <- 1 / sqrt(sum(weight))
```

For ratio measures, `gamma` must be on the log scale and oriented as active
control versus placebo. For example, negate a reported log placebo-versus-active
control ratio before using it. Then pass `hist_effect` and `hist_effect_se`, or
their natural-scale equivalent where requested, to the appropriate design
function.

A fixed-effect summary is appropriate only when its common-effect and
transportability assumptions are scientifically defensible. When historical
effects are heterogeneous or not directly transportable to the target trial,
derive an appropriate historical estimate and uncertainty externally and
document that choice before using `nidesigns`.

## Reproduce the paper example

This is the conditional-power design in Appendix B. The output matches the
published non-inferiority margins, required event counts, sample sizes,
controlled non-constancy values, and unconditional powers.

```r
paper_design <- ni_design(
  u.list = c(1, 1, 1 / (1 - 0.23), 0, 0),
  l1.list = c(
    0,
    -0.23,
    -0.23,
    1.96 * 0.61 / log(1 - 0.928),
    0
  ),
  f.preserv = 0.5,
  null.pe = 0.3,
  design.alternative.pe = 0.95,
  hist.ac.pe = 0.928,
  hist.ac.effect.se = 0.61,
  lambda0.for.design = 0,
  target.on.unconditional.power = FALSE,
  allocation.ratio = 1,
  power = 0.9,
  sign.level = 0.025,
  lambda0.sens.analysis = 0.12,
  placebo.incidence.rate = 0.03,
  loss.to.followup = 0.075,
  trial.duration = 2,
  correction = FALSE
)

summary(paper_design)
```

Use `ni_design_hr()` when it is more natural to provide hazard ratios directly.
The paper functions retain their event-driven, fixed-follow-up approximation
for translating required events into randomized sample size.

The `correction = TRUE` option retains the manuscript's one-event-per-arm
adjustment. It is not a complete group-sequential monitoring procedure.

## Additional effect measures

The framework requires an additive analysis scale on which smaller values
indicate benefit. The new functions use the following scales:

| Function | Analysis scale | Sample-size variance |
|---|---|---|
| `ni_design_rr()` | log relative risk | `(1-pX)/(nX*pX) + (1-pC)/(nC*pC)` |
| `ni_design_rd()` | risk difference | `pX*(1-pX)/nX + pC*(1-pC)/nC` |
| `ni_design_md()` | mean difference | `sdX^2/nX + sdC^2/nC` |

Here, `X` denotes the experimental arm and `C` the active-control arm. For the
binary measures, risks refer to a prespecified analysis time. The functions
first solve the paper's power equation for the required variance and then find
integer arm sizes whose achieved variance is no larger than that target.

### Relative risk

The scientific and operational hypotheses are evaluated on the log-RR scale.
The returned non-inferiority margin is transformed back to a relative risk.

```r
rr_design <- ni_design_rr(
  u.list = c(1, 0),
  l1.list = c(0, 0),
  f.preserv = 0.5,
  null.rr = 0.80,
  design.alternative.rr = 0.50,
  hist.ac.rr = 0.65,
  hist.ac.effect.se = 0.12,  # standard error of log(RR)
  placebo.risk = 0.20,
  allocation.ratio = 1,
  power = 0.90,
  dropout = 0.05
)

summary(rr_design)
```

### Risk difference

Risk differences are experimental minus placebo or active control. Negative
values therefore represent lower event risk and benefit.

```r
rd_design <- ni_design_rd(
  u.list = c(1, 0),
  l1.list = c(0, 0),
  f.preserv = 0.5,
  null.rd = -0.02,
  design.alternative.rd = -0.10,
  hist.ac.rd = -0.07,
  hist.ac.effect.se = 0.02,
  placebo.risk = 0.20,
  power = 0.90
)
```

### Mean difference

Mean differences are experimental minus placebo or active control. If larger
outcomes are favorable, multiply all means and mean differences by `-1` before
using this function so that smaller values indicate benefit.

```r
md_design <- ni_design_md(
  u.list = c(1, 0),
  l1.list = c(0, 0),
  f.preserv = 0.5,
  null.md = -0.5,
  design.alternative.md = -2.0,
  hist.ac.md = -1.5,
  hist.ac.effect.se = 0.4,
  sd.experimental = 5,
  sd.control = 5,
  power = 0.90
)
```

For a target variance obtained elsewhere, use `ni_sample_size()` directly:

```r
ni_sample_size(
  required.variance = 0.04,
  effect.measure = "mean_difference",
  allocation.ratio = 2,
  dropout = c(0.05, 0.10),
  sd.experimental = 2,
  sd.control = 2.5
)
```

## Assumptions for the new sample-size calculations

- Historical and active-trial effect estimators are independent and
  approximately normal on the stated additive scale.
- Binary outcomes use independent binomial sampling and risks fixed at a
  prespecified analysis time.
- Mean-difference calculations treat the planning standard deviations as
  known inputs and use the large-sample normal approximation.
- `allocation.ratio` is `n_experimental / n_control`.
- Dropout is non-informative. It inflates randomized counts but does not change
  the analyzable-sample variance.
- The returned sample sizes are planning approximations. Very rare outcomes,
  small samples, covariate adjustment, clustering, repeated measures, or
  informative missingness require a model-specific calculation or simulation.

## Compatibility and intentional improvements

The package has automated regression tests for the complete Appendix B example.
Default calls used in the paper reproduce the published results. The following
improvements intentionally affect only cases outside that published example:

- `sign.level` and `power` are now propagated consistently to variance, event,
  margin, and controlled-non-constancy calculations. The original script used
  hard-coded defaults in some downstream calls.
- Solving for required variance no longer changes the user's global warning
  option.
- Inputs now receive explicit range, length, and feasibility checks instead of
  failing later with obscure numerical errors.
- Nonstandard fixed-margin method labels now report the supplied `lambda1`
  value correctly.
- The maximum-unconditional-power calculation now uses the manuscript's
  definition of the historical variance for fixed-margin methods (`u = 0`).
- The new effect-measure functions report both target and achieved variance and
  mark unattainable designs explicitly.

## Legacy script interface

The original interface remains available for users who previously cloned the
repository:

- `ni_designs_fns.R` contains the original design and auxiliary functions.
- `example_ni_designs.R` contains the original examples.

Existing scripts may continue to use:

```r
source("ni_designs_fns.R")
```

The package interface is recommended for new work because it adds documented
functions, input validation, regression tests, and standard installation.

The original root-level scripts remain in the repository for provenance and
backward compatibility. Existing users may continue to source those scripts;
the package implementation lives under `R/`, documentation under `man/`, and
regression tests under `tests/testthat/`.

Pulling this repository does not modify an already installed copy of the
package or a running R session. Reinstall `nidesigns` to use a newly pulled
package version. Git will also stop and request conflict resolution rather than
silently overwrite a user's local edits to tracked files.

## Development

Run the test suite and package checks with:

```r
devtools::test()
devtools::check()
```
