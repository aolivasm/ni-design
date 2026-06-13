# Non-Inferiority Trial Design R Functions

## Background

Active-controlled non-inferiority trials are commonly used when an effective intervention already exists and a placebo-controlled design is not feasible or ethical. In these trials, a new experimental intervention is compared against an active control rather than placebo. The goal is to show that the experimental intervention is not unacceptably worse than the active control, while potentially offering other advantages such as improved safety, tolerability, cost, or feasibility.

A fundamental challenge in active-controlled trials is the absence of a placebo arm. As a result, assessing non-inferiority necessarily relies on historical estimates of the active control effect relative to placebo, along with assumptions about how well those estimates apply to the current trial population.

The traditional non-inferiority criterion is based on _preservation of effect_, which requires the experimental intervention to preserve a specified fraction (often 50%) of the active control’s effect relative to placebo.

When the active control is highly effective, preservation-of-effect criteria may impose an unrealistically high bar, even for interventions with meaningful public health value. This has motivated alternative success criteria, such as _inferred efficacy_, which assess whether the experimental intervention would have met a predefined efficacy threshold relative to a hypothetical placebo.

This repository provides R code for computing non-inferiority trial design operating characteristics under the preservation-of-effect and inferred-efficacy criteria when the outcome is time to event and intervention effects are measured in the log hazard ratio scale. The main function, `ni_design()`, computes non-inferiority margins, required numbers of events, sample sizes, and unconditional power for a set of analytical methods defined by `(u, lambda1)` pairs.

The examples included in this repository reproduce the results presented in the manuscript [Olivas-Martinez, Gao, and Janes (Stat Med, 2026)](https://onlinelibrary.wiley.com/doi/10.1002/sim.70618), which develops the underlying methodological framework and how the non-inferiority criteria and different analytical approaches are parameterized and interpreted.

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

## Citation

If you use the `ni_design()` function or this repository in your research, please cite the associated paper:

A. Olivas-Martinez, F. Gao, and H. Janes, “A General Framework for Designing and Evaluating Active-Controlled Trials with Non-Inferiority Objectives,” Statistics in Medicine 45, no. 13-14 (2026): e70618, https://doi.org/10.1002/sim.70618.
