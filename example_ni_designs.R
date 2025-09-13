rm(list = ls())
source("ni_designs_fns.R")

### Main function for trial design: ni_designs()
### Inputs:
###  - v.list: list of v parameters for the desired analytical methods
###  - w.list: list of w parameters for the desired analytical methods
###  - l1.list: list of l1 parameters for the desired analytical methods
###    + (v, w, l1) = (1, 1, 0) for the 95-95 method
###    + (v, w, l1) = (0, 1, 0) for the traditional synthesis method
###  - f.preserv: fraction for the preservation of effect criterion
###  - null.pe: null efficacy for the inference criterion
###  - design.alternative.pe: true efficacy of experimental relative to placebo (design alternative).
###      Note: This must be given in the prevention efficacy (PE) scale.
###  - hist.ac.pe: historical AC effect estimate
###      Note: This must be given in the prevention efficacy (PE) scale.
###  - hist.ac.effect.se: standard error of the historical AC effect estimate
###      Note: This is for the estimate in the log HR scale.
###  - lambda0.for.design: true fraction of the historical AC effect that has been lost.
###  - target.on.unconditional.power: TRUE if the design is based on unconditional power.
###  - allocation.ratio: allocation ratio of participants of experimental versus control groups
###  - power: desired power. Default is 0.9.
###  - sign.level: significance level. Default is 0.025.
###  - lambda0.sens.analysis: lambda0 for sensitivity analysis of unconditional power.
###  - placebo.incidence.rate: incidence rate of disease in untreated placebo population
###  - loss.to.followup: % of loss to follow-up per year.
###  - trial.duration: duration of the trial in years.
###  - correction: TRUE if the sample size computation is corrected to account for interim monitoring.

### Output: A data frame with the trial design specification for each method.
###         The columns are:
###         - Method: name of the method
###         - NI margin: non-inferiority margin
###         - RNE: required number of events
###         - Exp: required number of events in experimental arm
###         - Ctr: required number of events in control arm
###         - Sample size: required sample size
###         - Exp.arm: required sample size in experimental arm
###         - Ctr.arm: required sample size in control arm
###         - CNC: control non-constancy or minimum active control PE for which TIE is controlled
###         - U.power: unconditional power
###         - U.power (SA): unconditional power for sensitivity analysis of lambda0

### Examples used in the paper
##### Controlling conditional power and assuming CAB-LA PE = 92.8% (lambda0 = 0)
##### Sensitivity analysis of lambda0 = 12% (lambda0 = 0.12)

obj.designs.app1 = ni_design(u.list = c(1, 1, 1/(1-0.23), 0, 0),
                              l1.list = c(0, -0.23, -0.23, 1.96*0.61/log(1-0.928), 0),
                              f.preserv = 0.5,
                              null.pe = 0.3,
                              design.alternative.pe = 0.95,
                              hist.ac.pe = 0.928,
                              hist.ac.effect.se = 0.61,
                              lambda0.for.design = 0,
                              target.on.unconditional.power = FALSE,
                              allocation.ratio = 1,
                              power = 0.9, sign.level = 0.025,
                              lambda0.sens.analysis = 0.12,
                              placebo.incidence.rate = 0.03,
                              loss.to.followup = 0.075,
                              trial.duration = 2,
                              correction = FALSE)

##### Controlling unconditional power and assuming CAB-LA PE = 92.8% (lambda0 = 0)
##### Sensitivity analysis of lambda0 = 12% (lambda0 = 0.12)

obj.designs.app2 = ni_design(u.list = c(1, 1, 1/(1-0.23), 0, 0),
                              l1.list = c(0, -0.23, -0.23, 1.96*0.61/log(1-0.928), 0),
                              f.preserv = 0.5,
                              null.pe = 0.3,
                              design.alternative.pe = 0.95,
                              hist.ac.pe = 0.928,
                              hist.ac.effect.se = 0.61,
                              lambda0.for.design = 0,
                              target.on.unconditional.power = TRUE,
                              allocation.ratio = 1,
                              power = 0.9, sign.level = 0.025,
                              lambda0.sens.analysis = 0.12,
                              placebo.incidence.rate = 0.03,
                              loss.to.followup = 0.075,
                              trial.duration = 2,
                              correction = FALSE)

##### Controlling conditional power and assuming CAB-LA PE = 94.7% (lambda0 = 0.12)
##### Sensitivity analysis of lambda0 = 0% (lambda0 = 0)

obj.designs.app3 = ni_design(u.list = c(1, 1, 1/(1-0.23), 0, 0),
                              l1.list = c(0, -0.23, -0.23, 1.96*0.61/log(1-0.928), 0),
                              f.preserv = 0.5,
                              null.pe = 0.3,
                              design.alternative.pe = 0.95,
                              hist.ac.pe = 0.928,
                              hist.ac.effect.se = 0.61,
                              lambda0.for.design = 0.12,
                              target.on.unconditional.power = FALSE,
                              allocation.ratio = 1,
                              power = 0.9, sign.level = 0.025,
                              lambda0.sens.analysis = 0,
                              placebo.incidence.rate = 0.03,
                              loss.to.followup = 0.075,
                              trial.duration = 2,
                              correction = FALSE)

### Outputs
obj.designs.app1
obj.designs.app2
obj.designs.app3

obj.designs.app3$`NI criterion: Preserving 50% of active control effect`$`Sample size`/obj.designs.app1$`NI criterion: Preserving 50% of active control effect`$`Sample size`
obj.designs.app3$`NI criterion: Inferred efficacy of 30% relative to hypothetical placebo`$`Sample size`/obj.designs.app1$`NI criterion: Inferred efficacy of 30% relative to hypothetical placebo`$`Sample size`
#####################################################################################
#### Function to explore maximum unconditional power for different design alternatives

### Inputs:
###  - v.list: list of v parameters for the desired analytical methods
###  - w.list: list of w parameters for the desired analytical methods
###  - l1.list: list of l1 parameters for the desired analytical methods
###    + (v, w, l1) = (1, 1, 0) for the 95-95 method
###    + (v, w, l1) = (0, 1, 0) for the traditional synthesis method
###  - f.preserv: fraction for the preservation of effect criterion
###  - null.pe: null efficacy for the inference criterion
###  - design.alternatives.pe: range of design alternatives to explore.
###      Note: This must be given in the prevention efficacy (PE) scale.
###  - hist.ac.pe: historical AC effect estimate
###      Note: This must be given in the prevention efficacy (PE) scale.
###  - hist.ac.effect.se: standard error of the historical AC effect estimate
###      Note: This is for the estimate in the log HR scale.
###  - lambda0.for.design: true fraction of the historical AC effect that has been lost.
###  - sign.level: significance level. Default is 0.025.
###  - opt.power: desired optimal unconditional power. Default is 0.9.
###      Note: This is only to add an horizontal line to the plots.
###  - legend.position: position of the legend in the plots. Default is "bottom".

obj.max.up = explore_max_uncond_power(u.list = c(1, 1, 1/(1-0.23), 0, 0),
                                      l1.list = c(0, -0.23, -0.23, -0.45, 0),
                                      f.preserv = 0.5,
                                      null.pe = 0.3,
                                      design.alt.pe = seq(0.65, 0.98, by = 0.01),
                                      hist.ac.pe = 0.928,
                                      hist.ac.effect.se = 0.61,
                                      lambda0.for.design = 0,
                                      sign.level = 0.025, opt.power = 0.9,
                                      legend.position = "bottom")

obj.max.up$plot.preserv.eff
obj.max.up$plot.inf.eff
