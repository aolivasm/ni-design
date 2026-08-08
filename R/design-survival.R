#' Required number of events for a log-ratio design
#'
#' Converts the required variance from the general framework into numbers of
#' events under a proportional-hazards/log-rank approximation.
#'
#' @inheritParams solve_variance
#' @param k Experimental-to-control allocation ratio.
#' @return A length-three vector containing total, experimental-arm, and
#'   control-arm event counts. `"DNP"` means the design alternative is not
#'   detectable and `"PNR"` means the requested power is not reachable.
#' @name required_events
NULL

#' @rdname required_events
#' @export
rne_uncond <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                       delta0 = 0, u = 1, l1 = 0, alpha = 0.025,
                       power = 0.9, k = 1) {
  .assert_scalar_numeric(k, "k", lower = 0, lower_open = TRUE)
  max_pow <- max_uncond_power(g_cph, se_cph, l0, g_xp, f, delta0,
                              u, l1, alpha)
  if (is.na(max_pow)) return(c("DNP", "-", "-"))
  if (max_pow < power) return(c("PNR", "-", "-"))
  variance <- solve_v_xs_uncond(g_cph, se_cph, l0, g_xp, f, delta0,
                                u, l1, alpha, power)
  if (is.na(variance)) return(c("PNR", "-", "-"))
  .events_from_variance(variance, g_xp - (1 + l0) * g_cph, k)
}

#' @rdname required_events
#' @export
rne_cond <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                     delta0 = 0, u = 1, l1 = 0, alpha = 0.025,
                     power = 0.9, k = 1) {
  .assert_scalar_numeric(k, "k", lower = 0, lower_open = TRUE)
  variance <- solve_v_xs_cond(g_cph, se_cph, l0, g_xp, f, delta0,
                              u, l1, alpha, power)
  if (is.na(variance)) return(c("DNP", "-", "-"))
  .events_from_variance(variance, g_xp - (1 + l0) * g_cph, k)
}

.events_from_variance <- function(variance, g_xs, k) {
  events_control <- (1 + exp(-g_xs) / k) / variance
  events_experimental <- k * exp(g_xs) * events_control
  events_control <- round(events_control)
  events_experimental <- round(events_experimental)
  c(events_experimental + events_control,
    events_experimental, events_control)
}

#' Convert required events to randomized sample size
#'
#' This is the manuscript's fixed-follow-up approximation. It assumes constant
#' incidence and efficacy over follow-up.
#'
#' @param rne.obj Length-three event-count vector from [rne_uncond()] or
#'   [rne_cond()].
#' @param rate Annual placebo incidence rate.
#' @param ltfu Annual loss-to-follow-up proportion.
#' @param dur Trial duration in years.
#' @param g_xp Experimental-versus-placebo log hazard ratio.
#' @param g_cph Historical active-control-versus-placebo log hazard ratio.
#' @param l0 Relative effect deviation used for planning.
#' @param k Experimental-to-control allocation ratio.
#' @param correction Whether to apply the manuscript's one-event-per-arm
#'   adjustment for interim monitoring.
#' @return Total, experimental-arm, and control-arm sample sizes.
#' @export
ss_formula <- function(rne.obj, rate = 0.03, ltfu = 0.075, dur = 2,
                       g_xp = log(0.05), g_cph = log(0.072), l0 = 0,
                       k = 1, correction = FALSE) {
  .assert_scalar_numeric(rate, "rate", lower = 0, upper = 1,
                         lower_open = TRUE)
  .assert_scalar_numeric(ltfu, "ltfu", lower = 0, upper = 1,
                         upper_open = TRUE)
  .assert_scalar_numeric(dur, "dur", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(k, "k", lower = 0, lower_open = TRUE)
  if (!is.logical(correction) || length(correction) != 1L ||
      is.na(correction)) {
    stop("`correction` must be TRUE or FALSE.", call. = FALSE)
  }
  rne0 <- suppressWarnings(as.numeric(rne.obj))
  if (length(rne0) != 3L || anyNA(rne0)) return(c(NA_real_, NA_real_, NA_real_))
  events_exp_scaled <- round(rne0[2] / k)
  events_control <- rne0[3]
  if (correction) {
    events_exp_scaled <- events_exp_scaled + 1
    events_control <- events_control + 1
  }
  rate_exp <- rate * exp(g_xp) * (1 - ltfu) * dur
  rate_control <- rate * exp((1 + l0) * g_cph) * (1 - ltfu) * dur
  n_control <- round(max(events_control / rate_control,
                         events_exp_scaled / rate_exp))
  n_exp <- k * n_control
  c(n_exp + n_control, n_exp, n_control)
}

#' Design a non-inferiority trial on the prevention-efficacy scale
#'
#' This is the primary function described in Appendix B of
#' Olivas-Martinez, Gao, and Janes (2026). It preserves the manuscript's
#' argument names and output structure.
#'
#' @param u.list,l1.list Equal-length vectors of framework parameters.
#' @param f.preserv Fraction of the active-control effect to preserve, or
#'   `NULL` to omit this criterion.
#' @param null.pe Minimum prevention efficacy for the inferred-efficacy
#'   criterion, or `NULL` to omit it.
#' @param design.alternative.pe Experimental prevention efficacy under the
#'   design alternative.
#' @param hist.ac.pe Historical active-control prevention efficacy.
#' @param hist.ac.effect.se Standard error of the historical effect estimate
#'   on the log hazard-ratio scale.
#' @param lambda0.for.design True relative effect deviation assumed at design.
#' @param target.on.unconditional.power If `TRUE`, target unconditional power;
#'   otherwise target conditional power.
#' @param allocation.ratio Experimental-to-control allocation ratio.
#' @param power Desired power.
#' @param sign.level One-sided significance level.
#' @param lambda0.sens.analysis Optional relative effect deviation for a power
#'   sensitivity analysis.
#' @param placebo.incidence.rate Annual placebo incidence.
#' @param loss.to.followup Annual loss-to-follow-up proportion.
#' @param trial.duration Fixed follow-up duration in years.
#' @param correction Apply the manuscript's one-event-per-arm adjustment for
#'   interim monitoring.
#' @return An object of class `ni.design` containing specifications and one
#'   data frame per requested success criterion.
#' @examples
#' design <- ni_design(
#'   u.list = c(1, 0),
#'   l1.list = c(0, 1.96 * 0.61 / log(1 - 0.928)),
#'   design.alternative.pe = 0.95,
#'   hist.ac.pe = 0.928,
#'   hist.ac.effect.se = 0.61
#' )
#' summary(design)
#' @export
ni_design <- function(u.list = c(1, 0, 0, 0),
                      l1.list = c(0, 0, -0.25, -0.25),
                      f.preserv = 0.5,
                      null.pe = 0.3,
                      design.alternative.pe = 0.85,
                      hist.ac.pe = 0.55,
                      hist.ac.effect.se = 0.25,
                      lambda0.for.design = 0,
                      target.on.unconditional.power = TRUE,
                      allocation.ratio = 1,
                      power = 0.9,
                      sign.level = 0.025,
                      lambda0.sens.analysis = NULL,
                      placebo.incidence.rate = 0.03,
                      loss.to.followup = 0.075,
                      trial.duration = 2,
                      correction = FALSE) {
  .assert_probability(design.alternative.pe, "design.alternative.pe",
                      open = TRUE)
  .assert_probability(hist.ac.pe, "hist.ac.pe", open = TRUE)
  if (!is.null(null.pe)) .assert_probability(null.pe, "null.pe", open = TRUE)
  g_xp <- log(1 - design.alternative.pe)
  g_cph <- log(1 - hist.ac.pe)
  display_effect <- function(l0) 1 - exp((1 + l0) * g_cph)
  criteria <- list()
  if (!is.null(f.preserv)) {
    criteria[[paste0("NI criterion: Preserving ",
                     round(100 * f.preserv, 2),
                     "% of active control effect")]] <- list(
      f = f.preserv, delta0 = 0
    )
  }
  if (!is.null(null.pe)) {
    criteria[[paste0("NI criterion: Inferred efficacy of ",
                     round(100 * null.pe, 2),
                     "% relative to hypothetical placebo")]] <- list(
      f = 0, delta0 = log(1 - null.pe)
    )
  }
  .design_ratio_engine(
    u.list = u.list, l1.list = l1.list, criteria = criteria,
    g_xp = g_xp, g_cph = g_cph, se_cph = hist.ac.effect.se,
    l0 = lambda0.for.design,
    target_unconditional = target.on.unconditional.power,
    k = allocation.ratio, power = power, alpha = sign.level,
    l0_sensitivity = lambda0.sens.analysis,
    rate = placebo.incidence.rate, ltfu = loss.to.followup,
    duration = trial.duration, correction = correction,
    display_effect = display_effect, display_name = "efficacy",
    cnc_transform = function(x) 1 - exp((1 + x) * g_cph)
  )
}

#' Design a non-inferiority trial using hazard ratios
#'
#' A hazard-ratio parameterization of [ni_design()].
#'
#' @inheritParams ni_design
#' @param null.hr Maximum experimental-versus-placebo hazard ratio allowed by
#'   the inferred-efficacy criterion, or `NULL`.
#' @param design.alternative.hr Experimental-versus-placebo hazard ratio.
#' @param hist.ac.hr Historical active-control-versus-placebo hazard ratio.
#' @return An object of class `ni.design`.
#' @export
ni_design_hr <- function(u.list = c(1, 0, 0, 0),
                         l1.list = c(0, 0, -0.25, -0.25),
                         f.preserv = 0.5,
                         null.hr = 0.7,
                         design.alternative.hr = 0.15,
                         hist.ac.hr = 0.45,
                         hist.ac.effect.se = 0.25,
                         lambda0.for.design = 0,
                         target.on.unconditional.power = TRUE,
                         allocation.ratio = 1,
                         power = 0.9,
                         sign.level = 0.025,
                         lambda0.sens.analysis = NULL,
                         placebo.incidence.rate = 0.03,
                         loss.to.followup = 0.075,
                         trial.duration = 2,
                         correction = FALSE) {
  .assert_ratio(design.alternative.hr, "design.alternative.hr")
  .assert_ratio(hist.ac.hr, "hist.ac.hr")
  if (!is.null(null.hr)) .assert_ratio(null.hr, "null.hr")
  g_xp <- log(design.alternative.hr)
  g_cph <- log(hist.ac.hr)
  display_effect <- function(l0) exp((1 + l0) * g_cph)
  criteria <- list()
  if (!is.null(f.preserv)) {
    criteria[[paste0("NI criterion: Preserving ",
                     round(100 * f.preserv, 2),
                     "% of active control effect")]] <- list(
      f = f.preserv, delta0 = 0
    )
  }
  if (!is.null(null.hr)) {
    criteria[[paste0("NI criterion: Inferred hazards ratio of ",
                     round(null.hr, 2),
                     " relative to hypothetical placebo")]] <- list(
      f = 0, delta0 = log(null.hr)
    )
  }
  .design_ratio_engine(
    u.list = u.list, l1.list = l1.list, criteria = criteria,
    g_xp = g_xp, g_cph = g_cph, se_cph = hist.ac.effect.se,
    l0 = lambda0.for.design,
    target_unconditional = target.on.unconditional.power,
    k = allocation.ratio, power = power, alpha = sign.level,
    l0_sensitivity = lambda0.sens.analysis,
    rate = placebo.incidence.rate, ltfu = loss.to.followup,
    duration = trial.duration, correction = correction,
    display_effect = display_effect, display_name = "effect",
    cnc_transform = function(x) exp((1 + x) * g_cph)
  )
}

.design_ratio_engine <- function(u.list, l1.list, criteria, g_xp, g_cph,
                                 se_cph, l0, target_unconditional, k,
                                 power, alpha, l0_sensitivity, rate, ltfu,
                                 duration, correction, display_effect,
                                 display_name, cnc_transform) {
  .validate_design_vectors(u.list, l1.list)
  .assert_scalar_numeric(se_cph, "hist.ac.effect.se", lower = 0,
                         lower_open = TRUE)
  .assert_scalar_numeric(k, "allocation.ratio", lower = 0,
                         lower_open = TRUE)
  .assert_scalar_numeric(power, "power", lower = 0.5, upper = 1,
                         lower_open = TRUE, upper_open = TRUE)
  .assert_scalar_numeric(alpha, "sign.level", lower = 0, upper = 0.5,
                         lower_open = TRUE, upper_open = TRUE)
  if (length(criteria) == 0L) {
    stop("At least one non-inferiority criterion must be requested.",
         call. = FALSE)
  }
  method.list <- mapply(method_name, u.list, l1.list,
                        MoreArgs = list(g_cph = g_cph, se_cph = se_cph),
                        USE.NAMES = FALSE)
  tables <- lapply(criteria, function(criterion) {
    .build_survival_table(
      method.list, u.list, l1.list, g_xp, g_cph, se_cph, l0,
      criterion$f, criterion$delta0, target_unconditional, k, power,
      alpha, l0_sensitivity, rate, ltfu, duration, correction,
      cnc_transform
    )
  })
  type_power <- if (target_unconditional) "unconditional power" else
    "conditional power"
  approach <- paste0("Design approach targeting ", round(power * 100, 2),
                     "% ", type_power, " and assuming an active control ",
                     display_name, " of ",
                     round(100 * display_effect(l0), 1),
                     if (display_name == "efficacy") "%" else "")
  if (!is.null(l0_sensitivity) && l0_sensitivity != l0) {
    sensitivity <- paste0(
      "Sensitivity analysis (SA) assumes an active control ", display_name,
      " of ", round(100 * display_effect(l0_sensitivity), 1),
      if (display_name == "efficacy") "%" else ""
    )
    specifications <- list(
      Approach = approach,
      `Sensitivity analysis` = sensitivity
    )
  } else {
    specifications <- approach
    names(specifications) <- "Approach"
  }
  out <- c(list(Specifications = specifications), tables)
  class(out) <- "ni.design"
  invisible(out)
}

.build_survival_table <- function(methods, u, l1, g_xp, g_cph, se_cph,
                                  l0, f, delta0, target_unconditional, k,
                                  power, alpha, l0_sensitivity, rate, ltfu,
                                  duration, correction, cnc_transform) {
  solver <- if (target_unconditional) solve_v_xs_uncond else solve_v_xs_cond
  event_solver <- if (target_unconditional) rne_uncond else rne_cond
  variance <- mapply(
    solver, u = u, l1 = l1,
    MoreArgs = list(g_cph = g_cph, se_cph = se_cph, l0 = l0,
                    g_xp = g_xp, f = f, delta0 = delta0,
                    alpha = alpha, power = power),
    USE.NAMES = FALSE
  )
  se_xs <- sqrt(variance)
  rne <- t(mapply(
    event_solver, u = u, l1 = l1,
    MoreArgs = list(g_cph = g_cph, se_cph = se_cph, l0 = l0,
                    g_xp = g_xp, f = f, delta0 = delta0,
                    alpha = alpha, power = power, k = k),
    USE.NAMES = FALSE
  ))
  margins <- mapply(
    ni_margin, se_xs = se_xs, u = u, l1 = l1,
    MoreArgs = list(g_cph = g_cph, se_cph = se_cph, f = f,
                    d0 = delta0, alpha = alpha),
    USE.NAMES = FALSE
  )
  sample_size <- t(apply(
    rne, 1,
    function(x) ss_formula(x, rate, ltfu, duration, g_xp, g_cph, l0,
                           k, correction)
  ))
  tolerance <- vapply(seq_along(se_xs), function(i) {
    if (!is.finite(se_xs[i]) || se_xs[i] <= 0) return(NA_real_)
    min_l0(g_cph, se_cph, se_xs[i], f, u[i], l1[i], alpha)
  }, numeric(1))
  unconditional <- vapply(seq_along(se_xs), function(i) {
    if (!is.finite(se_xs[i]) || se_xs[i] <= 0) return(NA_real_)
    uncond_power(g_cph, se_cph, l0, g_xp, se_xs[i], f, delta0,
                 u[i], l1[i], alpha)
  }, numeric(1))
  out <- data.frame(
    Method = methods,
    `NI margin` = margins,
    RNE = rne[, 1],
    Exp = rne[, 2],
    Ctr = rne[, 3],
    `Sample size` = sample_size[, 1],
    Exp.arm = sample_size[, 2],
    Ctr.arm = sample_size[, 3],
    CNC = round(cnc_transform(tolerance), 3),
    U.power = round(unconditional, 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  for (j in 3:5) out[[j]] <- utils::type.convert(out[[j]], as.is = TRUE)
  if (!is.null(l0_sensitivity) && l0_sensitivity != l0) {
    sensitivity_power <- vapply(seq_along(se_xs), function(i) {
      if (!is.finite(se_xs[i]) || se_xs[i] <= 0) return(NA_real_)
      uncond_power(g_cph, se_cph, l0_sensitivity, g_xp, se_xs[i], f,
                   delta0, u[i], l1[i], alpha)
    }, numeric(1))
    out[["U.power (SA)"]] <- round(sensitivity_power, 2)
  }
  out
}

#' @export
print.ni.design <- function(x, ...) {
  cat("\n=== Non-Inferiority Design Results ===\n")
  tables <- x[setdiff(names(x), "Specifications")]
  for (nm in names(tables)) {
    cat("\n---", nm, "---\n")
    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::kable(tables[[nm]], align = "c"))
    } else {
      print(tables[[nm]], row.names = FALSE)
    }
  }
  invisible(x)
}

#' @export
summary.ni.design <- function(object, ...) {
  cat("\n=== Summary of Non-Inferiority Trial Design ===\n")
  cat("\nDesign Specifications:\n")
  specs <- object$Specifications
  if (is.list(specs)) {
    for (nm in names(specs)) cat(" *", nm, ":", specs[[nm]], "\n")
  } else {
    cat(" *", names(specs), ":", specs, "\n")
  }
  tables <- object[setdiff(names(object), "Specifications")]
  for (nm in names(tables)) {
    cat("\n---", nm, "---\n")
    preview <- utils::head(tables[[nm]], 5)
    if (requireNamespace("knitr", quietly = TRUE)) {
      print(knitr::kable(preview, align = "c"))
    } else {
      print(preview, row.names = FALSE)
    }
  }
  invisible(object)
}

.validate_design_vectors <- function(u, l1) {
  if (!is.numeric(u) || !is.numeric(l1) || length(u) == 0L ||
      length(u) != length(l1) || anyNA(u) || anyNA(l1) ||
      any(!is.finite(u)) || any(!is.finite(l1))) {
    stop("`u.list` and `l1.list` must be equal-length finite numeric vectors.",
         call. = FALSE)
  }
  if (any(u < 0)) stop("All `u.list` values must be non-negative.",
                       call. = FALSE)
  if (any(l1 <= -1)) stop("All `l1.list` values must be greater than -1.",
                          call. = FALSE)
  invisible(NULL)
}

.assert_probability <- function(x, name, open = FALSE) {
  .assert_scalar_numeric(x, name, lower = 0, upper = 1,
                         lower_open = open, upper_open = open)
}

.assert_ratio <- function(x, name) {
  .assert_scalar_numeric(x, name, lower = 0, lower_open = TRUE)
}
