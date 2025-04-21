###############################################################################
######################### Functions for ni_designs ############################

# Install the following package if not already installed
#install.packages("nleqslv") # This package is required for the function solve_v_xs_uncond and solve_v_xs_cond
library(nleqslv)
library(ggplot2)
library(dplyr)

###############################################################################
########## Type I error and tolerable level of non-constancy ##################

### Inputs:
###  - g_sp_h: historical AC effect estimate
###  - se_sp: standard error of the historical AC effect estimate
###  - l0: fraction of the historical AC effect that has been lost.
###        l0=0 means no loss of AC effect.
###  - se_xs: standard error of the effect estimate of experimental relative to
###           active control.
###  - f: fraction of the AC effect that want to be preserved. If not provided,
###       f = 0.5 by default.
###       For inferred efficacy margins, set f=0.
###  - v: tool parameter to navigate between synthesis (v=0) and fixed-margin (v=1) methods.
###  - w: tool parameter to choose a particular fixed-margin method.
###       w=1 for the 95-95 method.
###  - l1: Assumed fraction of the AC effect that has been lost.
###       l1=0 for the 95-95 method.
###  - alpha: significance level. Default is 0.025.

### Function to calculate the unconditional type I error rate
uncond_t1e = function(g_sp_h, se_sp, l0 = 0, se_xs, f = 0.5, v =1, w = 1, l1 = 0, alpha = 0.025){
  if(is.numeric(se_sp)&is.numeric(se_xs)&is.numeric(f)){
    bias = (1-f)*(l0-l1)*g_sp_h
    v_tilde = se_xs^2+2*v*w*(1-f)*(1+l1)*se_xs*se_sp+w^2*(1-f)^2*(1+l1)^2*se_sp^2
    numer = bias-sqrt(v_tilde)*qnorm(1-alpha)
    denom = sqrt(se_xs^2+(1-f)^2*(1+l1)^2*se_sp^2)
    pnorm(numer/denom)
  }
  else{
    return("A numerical value must be provided for se_sp, se_xs, delta, g_xs, and f")
  }
}

### Function to calculate the conditional type I error rate
cond_t1e = function(g_sp_h, se_sp, l0 = 0, se_xs, f = 0.5, v =1, w = 1, l1 = 0, alpha = 0.025){
  if(is.numeric(se_sp)&is.numeric(se_xs)&is.numeric(f)){
    bias = (1-f)*(l0-l1)*g_sp_h
    v_tilde = se_xs^2+2*v*w*(1-f)*(1+l1)*se_xs*se_sp+w^2*(1-f)^2*(1+l1)^2*se_sp^2
    numer = bias-sqrt(v_tilde)*qnorm(1-alpha)
    pnorm(numer/se_xs)
  }
  else{
    return("A numerical value must be provided for se_sp, se_xs, delta, g_xs, and f")
  }
}

### Function to calculate the tolerable level of non-constancy
min_l0 = function(se_sp, g_sp_h, se_xs, f = 0.5, v =1, w = 1, l1 = 0, alpha = 0.025){
  if(f==1){
    return(0)
  }else{
    v_tilde = se_xs^2+2*v*w*(1-f)*(1+l1)*se_xs*se_sp+w^2*(1-f)^2*(1+l1)^2*se_sp^2
    numer = sqrt(v_tilde)-sqrt(se_xs^2+(1-f)^2*(1+l1)^2*se_sp^2)
    l1 + numer*qnorm(1-alpha)/(1-f)/g_sp_h
  }
}

#####################################################################################
################# Power and maximum unconditional power ########################

### Inputs:
###  - g_sp_h: historical AC effect estimate
###  - se_sp: standard error of the historical AC effect estimate
###  - l0: fraction of the historical AC effect that has been lost.
###        l0=0 means no loss of AC effect.
###  - g_xp: True effect of experimental relative to placebo (design alternative).
###  - se_xs: standard error of the effect estimate of experimental relative to
###           active control.
###  - f: fraction of the AC effect that want to be preserved. If not provided,
###       f = 0.5 by default.
###       For inferred efficacy margins, set f=0.
###  - delta0: maximum null efficacy. delta0 = 0 for preservation of effect margins.
###  - v: tool parameter to navigate between synthesis (v=0) and fixed-margin (v=1) methods.
###  - w: tool parameter to choose a particular fixed-margin method.
###       c=1 for the 95-95 method.
###  - l1: Assumed fraction of the AC effect that has been lost.
###       l1=0 for the 95-95 method.
###  - alpha: significance level. Default is 0.025.

uncond_power = function(g_sp_h, se_sp, l0 = 0, g_xp, se_xs, f = 0.5, delta0 =0, v =1, w = 1, l1 = 0, alpha = 0.025){
  v_tilde = se_xs^2+2*v*w*(1-f)*(1+l1)*se_xs*se_sp+w^2*(1-f)^2*(1+l1)^2*se_sp^2
  numer = delta0+((1+l0)-(1-f)*(1+l1))*g_sp_h-g_xp-sqrt(v_tilde)*qnorm(1-alpha)
  denom = sqrt(se_xs^2+(1-f)^2*(1+l1)^2*se_sp^2)
  pnorm(numer/denom)
}

cond_power = function(g_sp_h, se_sp, l0 = 0, g_xp, se_xs, f = 0.5, delta0 =0, v =1, w = 1, l1 = 0, alpha = 0.025){
  v_tilde = se_xs^2+2*v*w*(1-f)*(1+l1)*se_xs*se_sp+w^2*(1-f)^2*(1+l1)^2*se_sp^2
  numer = delta0+((1+l0)-(1-f)*(1+l1))*g_sp_h-g_xp-sqrt(v_tilde)*qnorm(1-alpha)
  pnorm(numer/se_xs)
}

max_uncond_power =  function(g_sp_h, se_sp, l0 = 0, g_xp, f = 0.5, delta0 =0, v =1, w = 1, l1 = 0, alpha = 0.025){
  a = delta0 +(1+l0-(1-f)*(1+l1))*g_sp_h-g_xp
  if(a-w*qnorm(1-alpha)*(1+l1)*(1-f)*se_sp<0){
    return(NA)
  }else{
    pnorm(-w*qnorm(1-alpha)+a/(1-f)/(1+l1)/se_sp)
  }
}

max_up_unrest =  function(g_sp_h, se_sp, l0 = 0, g_xp, f = 0.5, delta0 =0, v =1, w = 1, l1 = 0, alpha = 0.025){
  a = delta0 +(1-(1-f)*(1+l1)/(1+l0))*(1+l0)*g_sp_h-g_xp
  pnorm(-w*qnorm(1-alpha)+a/(1-f)/(1+l1)/se_sp)
}

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


explore_max_uncond_power = function(v.list = c(1, 0, 0, 0),
                                    w.list = c(1, 1, 1, 1/(1-0.23)),
                                    l1.list = c(0, 0, -0.23, -0.23),
                                    f.preserv = 0.5,
                                    null.pe = 0.3,
                                    design.alt.pe = seq(0.65, 0.98, by = 0.01),
                                    hist.ac.pe = 0.928,
                                    hist.ac.effect.se = 0.61,
                                    lambda0.for.design = 0,
                                    sign.level = 0.025, opt.power = 0.9,
                                    legend.position = "bottom"){
  
  g_xp = log(1-design.alt.pe)
  g_sp_h = log(1-hist.ac.pe)
  se_sp = hist.ac.effect.se
  l0 = lambda0.for.design
  
  method.list = mapply(FUN = method_name, v = v.list, w = w.list, l1 = l1.list, SIMPLIFY = TRUE)
  exp_pe_expanded = rep(design.alt.pe, times = length(v.list))
  method_expanded = rep(method.list, each = length(design.alt.pe))
  
  grid = expand.grid(
    i = seq_along(g_xp),
    j = seq_along(v.list)
  )
  
  
  if(!is.null(f.preserv)){
    max_up = mapply(function(i, j) {
      max_up_unrest(
        g_sp_h = g_sp_h,
        se_sp = se_sp,
        l0 = l0,
        g_xp = g_xp[i],
        f = f.preserv,
        delta0 = 0,
        v = v.list[j],
        w = w.list[j],
        l1 = l1.list[j],
        alpha = sign.level
      )
    }, i = grid$i, j = grid$j)
    
    df.preserv.eff = data.frame(
      Exp_pe = round(100*exp_pe_expanded,2),
      Method = method_expanded,
      Max_up = max_up
    )
    
    name_f = paste0("Criterion: Preserving ", round(100*f.preserv,2),"% of active control effect")
    
    plot.preserv.eff = df.preserv.eff %>%
      ggplot2::ggplot() +
      geom_hline(yintercept = opt.power, linetype = 2) +
      geom_line(aes(x = Exp_pe, y = Max_up, color = Method, linetype = Method), size =1.2) +
      theme_bw() + scale_color_brewer(palette = "Set1")+
      theme(legend.position = "bottom") + xlab("Efficacy of experimental intervention (%)") +
      ylab("Maximum unconditional power") + ggtitle(name_f) +
      theme(legend.position = legend.position)
  }
  
  if(!is.null(null.pe)){
    d0 = log(1-null.pe)
    max_up = mapply(function(i, j) {
      max_up_unrest(
        g_sp_h = g_sp_h,
        se_sp = se_sp,
        l0 = l0,
        g_xp = g_xp[i],
        f = 0,
        delta0 = d0,
        v = v.list[j],
        w = w.list[j],
        l1 = l1.list[j],
        alpha = sign.level
      )
    }, i = grid$i, j = grid$j)
    
    df.inf.eff = data.frame(
      Exp_pe = round(100*exp_pe_expanded,2),
      Method = method_expanded,
      Max_up = max_up
    )
    
    name_d0 = paste0("Criterion: Inferred efficacy of ", round(100*null.pe,2), "%.")
    
    plot.inf.eff = df.inf.eff %>%
      ggplot2::ggplot() +
      geom_hline(yintercept = opt.power, linetype = 2) +
      geom_line(aes(x = Exp_pe, y = Max_up, color = Method, linetype = Method), size =1.2) +
      theme_bw() + scale_color_brewer(palette = "Set1")+
      theme(legend.position = "bottom") + xlab("Efficacy of experimental intervention (%)") +
      ylab("Maximum unconditional power") + ggtitle(name_d0)+
      theme(legend.position = legend.position)
  }
  
  if(!is.null(f.preserv)&!is.null(null.pe)){
    return(list(table.preserv.eff = df.preserv.eff,
                plot.preserv.eff = plot.preserv.eff,
                table.inf.eff = df.inf.eff,
                plot.inf.eff = plot.inf.eff))
  }else if(!is.null(f.preserv)){
    return(list(table.preserv.eff = df.preserv.eff,
                plot.preserv.eff = plot.preserv.eff))
  }else {
    return(list(table.inf.eff = df.inf.eff,
                plot.inf.eff = plot.inf.eff))
  }
}

##############################################################################
######################### Required sample size ###############################

### Inputs:
###  - g_sp_h: historical AC effect estimate
###  - se_sp: standard error of the historical AC effect estimate
###  - l0: fraction of the historical AC effect that has been lost.
###        l0=0 means no loss of AC effect.
###  - g_xp: True effect of experimental relative to placebo (design alternative).
###  - f: fraction of the AC effect that want to be preserved. If not provided,
###       f = 0.5 by default.
###       For inferred efficacy margins, set f=0.
###  - delta0: maximum null efficacy. delta0 = 0 for preservation of effect margins.
###  - v: tool parameter to navigate between synthesis (v=0) and fixed-margin (v=1) methods.
###  - w: tool parameter to choose a particular fixed-margin method.
###       c=1 for the 95-95 method.
###  - l1: Assumed fraction of the AC effect that has been lost.
###       l1=0 for the 95-95 method.
###  - alpha: significance level. Default is 0.025.
###  - power: desired power. Default is 0.9.
###  - k: allocation ratio of participants of experimental versus control groups


###  - v_xs: required standard error of the effect estimate of experimental relative to
###           active control.

###  Approach targeting unconditional power
solve_v_xs_uncond = function(g_sp_h, se_sp, l0=0, g_xp, f = 0.5, delta0 =0, v =1, w=1, l1=0, alpha = 0.025, power = 0.9){
  max_pow = max_uncond_power(g_sp_h = g_sp_h, se_sp=se_sp, l0 = l0, g_xp = g_xp, f = f, delta0 =delta0, v =v, w = w, l1 = l1, alpha = alpha)
  options(warn = - 1)
  if(is.na(max_pow)){
    return(NA)
  }else{
    if(max_pow<power){
      return(NA)
    }else{
      if(f<1){
        c = (1-f)*(1+l1)*se_sp
        fn_aux = function(x){
          v_tilde = sqrt(x+2*v*w*c*sqrt(x)+w^2*c^2)
          qnorm(power) + (v_tilde*qnorm(1-alpha)-delta0+(1-f)*(1+l1)*g_sp_h+g_xp-(1+l0)*g_sp_h)/sqrt(x+c^2)
        }
        result = nleqslv::nleqslv(0.001, fn_aux)
        if(result$x<0){
          return(NA)
        }else{
          return(result$x) 
        }
      }else{
        v_xs =(delta0 - g_xp+(1+l0)*g_sp_h)^2/(qnorm(1-alpha)+qnorm(power))^2
        return(v_xs)
      }
    }
  }
}

###  Approach targeting conditional power
solve_v_xs_cond = function(g_sp_h, se_sp, l0=0, g_xp, f = 0.5, delta0 =0, v =1, w=1, l1=0, alpha = 0.025, power = 0.9){
  pow_ach = 1*(delta0>(1-f)*(1+l1)*(g_sp_h+w*se_sp*qnorm(1-alpha))+g_xp-(1+l0)*g_sp_h)
  if(pow_ach==1){
    if(f<1){
      c = (1-f)*(1+l1)*se_sp
      fn_aux = function(x){
        v_tilde = sqrt(x+2*v*w*c*sqrt(x)+w^2*c^2)
        qnorm(power) + (v_tilde*qnorm(1-alpha)-delta0+(1-f)*(1+l1)*g_sp_h+g_xp-(1+l0)*g_sp_h)/sqrt(x)
      }
      result = nleqslv::nleqslv(0.001, fn_aux)
      return(result$x)
    }else{
      v_xs =(delta0 - g_xp+(1+l0)*g_sp_h)^2/(qnorm(1-alpha)+qnorm(power))^2
      return(v_xs)
    }
  }else{
    return(NA)
  }
}

### Function to calculate the required number of events targeting unconditional power
rne_uncond = function(g_sp_h, se_sp, l0=0, g_xp, f = 0.5, delta0 =0, v =1, w=1, l1=0, alpha = 0.025, power = 0.9, k =1){
  max_pow = max_uncond_power(g_sp_h = g_sp_h, se_sp=se_sp, l0 = l0, g_xp = g_xp, f = f, delta0 =delta0, v =v, w = w, l1 = l1, alpha = alpha)
  if(is.na(max_pow)){
    return(c("DNP", "-", "-"))
  }else{
    if(max_pow<power){
      return(c("PNR", "-", "-"))
    }else{
      v_xs = solve_v_xs_uncond(g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0, g_xp = g_xp, f = f, delta0 =delta0, v =v, w=w, l1 = l1, alpha = alpha, power = power)
      if(is.na(v_xs)){
        return(c("PNR", "-", "-"))
      }else{
        g_xs = g_xp-(1+l0)*g_sp_h
        Ls = (1+exp(-g_xs)/k)/v_xs 
        Lx = k*exp(g_xs)*Ls
        Lsf = round(Ls)
        Lxf = round(Lx)
        return(c(Lsf+Lxf, Lxf, Lsf)) 
      }
    }
  }
}

### Function to calculate the required number of events targeting conditional power
rne_cond = function(g_sp_h, se_sp, l0=0, g_xp, f = 0.5, delta0 =0, v =1, w=1, l1=0, alpha = 0.025, power = 0.9, k =1){
  pow_ach = 1*(delta0>(1-f)*(1+l1)*(g_sp_h+w*se_sp*qnorm(1-alpha))+g_xp-(1+l0)*g_sp_h)
  if(pow_ach==1){
    v_xs = solve_v_xs_cond(g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0, g_xp = g_xp, f = f, delta0 =delta0, v =v, w=w, l1 = l1, alpha = alpha, power = power)
    if(!is.na(v_xs)){
      g_xs = g_xp-(1+l0)*g_sp_h
      Ls = (1+exp(-g_xs)/k)/v_xs 
      Lx = k*exp(g_xs)*Ls
      Lsf = round(Ls)
      Lxf = round(Lx)
      return(c(Lsf+Lxf, Lxf, Lsf))
    }else{
      return(c(rep(NA,3)))
    }
  }else{
    return(c("DNP", "-", "-"))
  }
}

### Function to compute required sample size given required number of events
ss_formula = function(rne.obj, rate = 0.03, ltfu = 0.075, dur = 2, g_xp = 0.95, g_sp_h = 0.928, l0=0, k = 1, correction = FALSE){
  
  rne0 = as.numeric(rne.obj)
  if(!is.na(rne0[1])){
    n1 = round(rne0[2]/k,0)
    n0 = rne0[3]
    if(correction){
      n1 = n1 + 1
      n0 = n0 + 1
    }
    
    HR1 = exp(g_xp)
    HR0 = exp((1+l0)*g_sp_h)
    
    rate1 = rate*HR1*(1-ltfu)*dur # rate of events in Experimental arm 
    rate0 = rate*HR0*(1-ltfu)*dur # rate of events in Control arm 
    n1/rate1 # required N for Experimental arm
    n0/rate0 # required N for Control arm
    nper = round(max(n0/rate0,n1/rate1),0)
    nper1 = k*nper
    ntot = nper1 + nper
    
    return(c(ntot, nper1, nper))
  }else{
    return(c(NA, NA, NA))
  }
}

### Non-inferiority margin
ni_margin = function(g_sp_h, se_sp, se_xs, f = 0.5, d0 =0, v =1, w=1, l1=0, alpha = 0.025){
  v_tilde = se_xs^2+2*v*w*(1+l1)*(1-f)*se_xs*se_sp+w^2*(1+l1)^2*(1-f)^2*se_sp^2
  round(exp(d0 - (1-f)*(1+l1)*g_sp_h-qnorm(1-alpha)*(sqrt(v_tilde)-se_xs)),2)
}

### Method name
method_name = function(v, w, l1){
  if(v==1 & w==1 & l1==0){
    return("95-95 method")
  }else if(v==0 & w==1 & l1==0){
    return("Traditional SM")
  }else if(v==0 & w==1 & l1!=0){
    return(paste0("BA-SM, \u03BB1=", round(l1*100, 2), "%"))
  }else if(v==0 & w==1/(1+l1) & l1!=0){
    return(paste0("OD, \u03BB1=", round(l1*100, 2), "%"))
  }else if(v==1 & w!=1 & l1==0){
    return(paste0("FM method with w=", round(w, 2)))
  }else return(paste0("Method with (v,w,l1)= (", v, ", ", w, ", ", round(l1, 2), ")"))
}

#####################################################################################
### Main function for trial design

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

ni_designs = function(v.list = c(1, 0, 0, 0),
                      w.list = c(1, 1, 1, 1/(1-0.25)),
                      l1.list = c(0, 0, -0.25, -0.25),
                      f.preserv = 0.5,
                      null.pe = 0.3,
                      design.alternative.pe = 0.85,
                      hist.ac.pe = 0.55,
                      hist.ac.effect.se = 0.25,
                      lambda0.for.design = 0,
                      target.on.unconditional.power = TRUE,
                      allocation.ratio = 1,
                      power = 0.9, sign.level = 0.025,
                      lambda0.sens.analysis = NULL,
                      placebo.incidence.rate = 0.03,
                      loss.to.followup = 0.075,
                      trial.duration = 2,
                      correction = FALSE
){
  
  g_xp = log(1-design.alternative.pe)
  g_sp_h = log(1-hist.ac.pe)
  se_sp = hist.ac.effect.se
  l0 = lambda0.for.design
  k = allocation.ratio
  rate = placebo.incidence.rate
  ltfu = loss.to.followup
  dur = trial.duration
  ac.pe = 1-exp(g_sp_h*(1+l0))
  
  if(!is.null(lambda0.sens.analysis)){
    l0.sa = lambda0.sens.analysis
    ac.pe.sa = 1-exp(g_sp_h*(1+l0.sa))
  }
  
  
  # Check if the lengths of the input vectors are equal
  if (length(v.list) != length(w.list) || length(v.list) != length(l1.list)) {
    stop("The lists for the parameters v, w, and l1 must have the same length.")
  }
  
  ### Create list with hyperparameters for trial design
  
  ### Preservation of effect criterion
  method.list = mapply(FUN = method_name, v = v.list, w = w.list, l1 = l1.list, SIMPLIFY = TRUE)
  
  if(!is.null(f.preserv)){
    if(target.on.unconditional.power){
      se_xs = sqrt(mapply(FUN = solve_v_xs_uncond, g_sp_h = g_sp_h, se_sp = se_sp, l0= l0, g_xp = g_xp,
                          f = f.preserv, delta0 =0, v =v.list, w=w.list, l1= l1.list, alpha = sign.level, power = power))
      rne = t(mapply(FUN = rne_uncond, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0,
                     g_xp = g_xp, f = f.preserv, delta0 = 0, v= v.list, w = w.list, l1 = l1.list, k = k))
    }else{
      se_xs = sqrt(mapply(FUN = solve_v_xs_cond, g_sp_h = g_sp_h, se_sp = se_sp, l0= l0, g_xp = g_xp,
                          f = f.preserv, delta0 =0, v =v.list, w=w.list, l1= l1.list, alpha = sign.level, power = power))
      rne = t(mapply(FUN = rne_cond, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0,
                     g_xp = g_xp, f = f.preserv, delta0 = 0, v= v.list, w = w.list, l1 = l1.list, k = k))
    }
    
    ni_margin = mapply(FUN = ni_margin, g_sp_h = g_sp_h, se_sp = se_sp, se_xs = se_xs, f = f.preserv, d0 =0, v =v.list, w=w.list, l1=l1.list, alpha = 0.025)
    sample.size = t(apply(rne, 1, FUN = function(x) ss_formula(rne.obj = x, rate = rate, ltfu = ltfu, dur = dur,
                                                               g_xp = g_xp, g_sp_h = g_sp_h, l0= l0, k = k, correction = correction)))
    max_tol_loss = mapply(FUN = min_l0, se_sp = se_sp, g_sp_h = g_sp_h, se_xs = se_xs, f = f.preserv, v =v.list, w=w.list, l1=l1.list, alpha = 0.025)
    up = mapply(FUN = uncond_power, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0, g_xp = g_xp,
                se_xs = se_xs, f = f.preserv, delta0 =0, v =v.list, w = w.list,
                l1 = l1.list, alpha = sign.level)
    
    if(!is.null(lambda0.sens.analysis)&lambda0.sens.analysis!=lambda0.for.design){
      up.sa = mapply(FUN = uncond_power, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0.sa, g_xp = g_xp,
                     se_xs = se_xs, f = f.preserv, delta0 =0, v =v.list, w = w.list,
                     l1 = l1.list, alpha = sign.level)
    }
    
    df.preserv = data.frame(
      method = method.list,
      margin = ni_margin,
      rne = rne[,1],
      exp_arm = rne[,2],
      ctr_arm = rne[,3],
      ss = sample.size[,1],
      exp = sample.size[,2],
      ctr = sample.size[,3],
      min_pe_t1e_control = round(1-exp((1+max_tol_loss)*g_sp_h),3),
      uncond_power = round(up,2)
    )
    
    colnames(df.preserv) = c("Method", "NI margin", "RNE", "Exp", "Ctr",
                             "Sample size", "Exp.arm", "Ctr.arm",
                             "CNC", "U.power")
    
    if(!is.null(lambda0.sens.analysis)&lambda0.sens.analysis!=lambda0.for.design){
      df.preserv = cbind(df.preserv, lambda0.sens.analysis = round(up.sa,2))
      colnames(df.preserv)[11] = "U.power (SA)"
    }
  }
  
  if(!is.null(null.pe)){
    d0 = log(1-null.pe)
    if(target.on.unconditional.power){
      se_xs = sqrt(mapply(FUN = solve_v_xs_uncond, g_sp_h = g_sp_h, se_sp = se_sp, l0= l0, g_xp = g_xp,
                          f = 0, delta0 =d0, v =v.list, w=w.list, l1= l1.list, alpha = sign.level, power = power))
      rne = t(mapply(FUN = rne_uncond, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0,
                     g_xp = g_xp, f = 0, delta0 = d0, v= v.list, w = w.list, l1 = l1.list, k = k))
    }else{
      se_xs = sqrt(mapply(FUN = solve_v_xs_cond, g_sp_h = g_sp_h, se_sp = se_sp, l0= l0, g_xp = g_xp,
                          f = 0, delta0 =d0, v =v.list, w=w.list, l1= l1.list, alpha = sign.level, power = power))
      rne = t(mapply(FUN = rne_cond, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0,
                     g_xp = g_xp, f = 0, delta0 = d0, v= v.list, w = w.list, l1 = l1.list, k = k))
    }
    
    ni_margin = mapply(FUN = ni_margin, g_sp_h = g_sp_h, se_sp = se_sp, se_xs = se_xs, f = 0, d0 =d0, v =v.list, w=w.list, l1=l1.list, alpha = 0.025)
    sample.size = t(apply(rne, 1, FUN = function(x) ss_formula(rne.obj = x, rate = rate, ltfu = ltfu, dur = dur,
                                                               g_xp = g_xp, g_sp_h = g_sp_h, l0= l0, k = k, correction = correction)))
    max_tol_loss = mapply(FUN = min_l0, se_sp = se_sp, g_sp_h = g_sp_h, se_xs = se_xs, f = 0, v =v.list, w=w.list, l1=l1.list, alpha = 0.025)
    up = mapply(FUN = uncond_power, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0, g_xp = g_xp,
                se_xs = se_xs, f = 0, delta0 =d0, v =v.list, w = w.list,
                l1 = l1.list, alpha = sign.level)
    
    if(!is.null(lambda0.sens.analysis)&lambda0.sens.analysis!=lambda0.for.design){
      up.sa = mapply(FUN = uncond_power, g_sp_h = g_sp_h, se_sp = se_sp, l0 = l0.sa, g_xp = g_xp,
                     se_xs = se_xs, f = 0, delta0 =d0, v =v.list, w = w.list,
                     l1 = l1.list, alpha = sign.level)
    }
    
    df.inf.eff = data.frame(
      method = method.list,
      margin = ni_margin,
      rne = rne[,1],
      exp_arm = rne[,2],
      ctr_arm = rne[,3],
      ss = sample.size[,1],
      exp = sample.size[,2],
      ctr = sample.size[,3],
      min_pe_t1e_control = round(1-exp((1+max_tol_loss)*g_sp_h),3),
      uncond_power = round(up,2)
    )
    
    colnames(df.inf.eff) = c("Method", "NI margin", "RNE", "Exp", "Ctr",
                             "Sample size", "Exp.arm", "Ctr.arm",
                             "CNC", "U.power")
    
    if(!is.null(lambda0.sens.analysis)&lambda0.sens.analysis!=lambda0.for.design){
      df.inf.eff = cbind(df.inf.eff, lambda0.sens.analysis = round(up.sa,2))
      colnames(df.inf.eff)[11] = "U.power (SA)"
    }
  }
  
  ### Combine results
  type.power = ifelse(target.on.unconditional.power, "unconditional power", "conditional power")
  design.approach = paste0("Design approach targeting ", round(power*100, 2), "% ", type.power, " and assuming an active control efficacy of ", round(100*ac.pe, 1), "%")
  if(!is.null(lambda0.sens.analysis)&lambda0.sens.analysis!=lambda0.for.design){
    sens.analysis = paste0("Sensitivity analysis (SA) assumes an active control efficacy of ", round(100*ac.pe.sa, 1), "%")
    trial.spec = list(design.approach, sens.analysis)
    names(trial.spec) = c("Approach", "Sensitivity analysis")
  }else{
    trial.spec = design.approach
    names(trial.spec) = "Approach"
  }
  
  if(!is.null(f.preserv)){
    name_f = paste0("NI criterion: Preserving ", round(100*f.preserv,2),"% of active control effect")
  }
  if(!is.null(null.pe)){
    name_d0 = paste0("NI criterion: Inferred efficacy of ", round(100*null.pe,2), "% relative to hypothetical placebo")
  }
  
  if(!is.null(f.preserv) & !is.null(null.pe)){
    obj.designs = list(Specifications = trial.spec, df.preserv = df.preserv, df.inf.eff = df.inf.eff)
    names(obj.designs) = c("Specifications", name_f, name_d0)
  }
  
  if(!is.null(f.preserv) & is.null(null.pe)){
    obj.designs = list(Specifications = trial.spec, df.preserv = df.preserv)
    names(obj.designs) = c("Specifications", name_f)
  }
  
  if(is.null(f.preserv) & !is.null(null.pe)){
    obj.designs = list(Specifications = trial.spec, df.inf.eff = df.inf.eff)
    names(obj.designs) = c("Specifications", name_d0)
  }
  
  return(obj.designs)
}
