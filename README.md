# Non-Inferiority Trial Design R Functions

This repository provides R code for generating non-inferiority trial design specifications and computing operating characteristics under the preservation-of-effect and inferred-efficacy criteria when the outcome is time to event and the effect is measured in the log hazard ratio scale. The main function is `ni_design()`, which computes non-inferiority margins, required numbers of events, sample sizes, and unconditional power for a set of analytical methods defined by `(u, lambda1)` pairs.

The examples used in this repository reproduce the results in the manuscript: [Olivas-Martinez et al., 2025, arXiv:2510.22071](https://arxiv.org/abs/2510.22071).

## Files

- `ni_designs_fns.R` – Contains the two main functions for non-inferiority trial design and all auxiliary functions needed to compute the operating characteristics described in the manuscript.
- `example_ni_designs.R` – Provides example scripts illustrating how to use the functions.

## Functions

### `ni_design()`

Generates non-inferiority trial designs and computes operating characteristics.

**Required inputs:**

- `u.list` – Vector of `u` parameters specifying the desired analytical methods.
- `l1.list` – Vector of `lambda1` parameters specifying the desired analytical methods. Must have the same length as `u.list`. Examples include:
  - `(u, lambda1) = (1, 0)` for the traditional synthesis method
  - `(u, lambda1) = (0, 1.96*sqrt(V_CP,H)/gamma_CP,H)` for the 95–95 method
- `design.alternative.pe` – True prevention efficacy (PE) of the experimental treatment relative to placebo under the design alternative (PE scale).
- `hist.ac.pe` – Historical estimate of active control PE.
- `hist.ac.effect.se` – Standard error of the historical active control effect estimate (log hazard ratio scale).

**Optional inputs and defaults:**

- `f.preserv` (default: 0.5) – Fraction of active control effect to preserve.
- `null.pe` (default: 0.3) – Null prevention efficacy for inferred-efficacy criterion.
- `lambda0.for.design` (default: 0) – Relative effect deviation assumed for the primary design.
- `target.on.unconditional.power` (default: TRUE) – If TRUE, designs target unconditional power.
- `allocation.ratio` (default: 1) – Experimental:control allocation ratio.
- `power` (default: 0.9) – Desired power.
- `sign.level` (default: 0.025) – One-sided significance level.
- `lambda0.sens.analysis` – Optional lambda0 for sensitivity analysis of unconditional power.
- `placebo.incidence.rate` (default: 0.03) – Annual incidence rate in the placebo population.
- `loss.to.followup` (default: 0.075) – Annual loss-to-follow-up proportion.
- `trial.duration` (default: 2) – Planned trial duration in years.
- `correction` (default: FALSE) – Apply correction for interim monitoring if TRUE.

**Outputs:**

- Returns an object of class `ni.design`, which is a list containing:
  - `Specifications` – Character vector summarizing the design approach and sensitivity analysis assumptions.
  - Data frames for each success criterion (e.g., preservation-of-effect, inferred efficacy), including:
    - `Method`, `NI margin`, `RNE`, `Exp`, `Ctr`, `Sample size`, `Exp.arm`, `Ctr.arm`, `CNC`, `U.power`, `U.power (SA)`.

- A custom `summary()` method displays the trial specifications followed by the design result tables.

### `explore_max_uncond_power()`

Explores the maximum unconditional power across a range of design alternatives and produces plots.

## Usage

All example code is in `example_ni_designs.R`. To run the examples:

```R
# Load the functions
source("ni_designs_fns.R")
```

Then run the example designs and plots in the `example_ni_designs.R` script.

## Requirements
- R (version >= 4.0)
- Any additional packages will be loaded within the scripts.
