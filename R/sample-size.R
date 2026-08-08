#' Convert a target effect-estimator variance to sample size
#'
#' Computes two-arm sample sizes using large-sample variances for the log
#' relative risk, risk difference, or mean difference. The experimental-to-
#' control allocation ratio is `n_experimental / n_control`.
#'
#' @param required.variance Maximum variance allowed for the estimated
#'   experimental-versus-control effect.
#' @param effect.measure One of `"relative_risk"`, `"risk_difference"`, or
#'   `"mean_difference"`.
#' @param allocation.ratio Experimental-to-control allocation ratio.
#' @param dropout Expected dropout proportion. Supply one value for both arms
#'   or two values ordered experimental, control.
#' @param risk.experimental,risk.control Event risks at the analysis time;
#'   required for binary effect measures.
#' @param sd.experimental,sd.control Outcome standard deviations; required for
#'   a mean difference.
#' @return A one-row data frame with analyzable and randomized sample sizes and
#'   the variance achieved by the analyzable sample.
#' @details
#' The formulas are
#' \deqn{Var(log(RR)) = (1-p_X)/(n_X p_X) + (1-p_C)/(n_C p_C),}
#' \deqn{Var(RD) = p_X(1-p_X)/n_X + p_C(1-p_C)/n_C,}
#' and
#' \deqn{Var(MD) = sigma_X^2/n_X + sigma_C^2/n_C.}
#' Dropout inflates randomized counts and is assumed non-informative; it does
#' not alter the analyzable-sample variance formula.
#' @examples
#' ni_sample_size(
#'   required.variance = 0.01,
#'   effect.measure = "relative_risk",
#'   risk.experimental = 0.10,
#'   risk.control = 0.15
#' )
#' ni_sample_size(
#'   required.variance = 0.04,
#'   effect.measure = "mean_difference",
#'   sd.experimental = 2,
#'   sd.control = 2
#' )
#' @export
ni_sample_size <- function(required.variance,
                           effect.measure = c("relative_risk",
                                              "risk_difference",
                                              "mean_difference"),
                           allocation.ratio = 1,
                           dropout = 0,
                           risk.experimental = NULL,
                           risk.control = NULL,
                           sd.experimental = NULL,
                           sd.control = NULL) {
  effect.measure <- match.arg(effect.measure)
  .assert_scalar_numeric(required.variance, "required.variance", lower = 0,
                         lower_open = TRUE)
  .assert_scalar_numeric(allocation.ratio, "allocation.ratio", lower = 0,
                         lower_open = TRUE)
  dropout <- .validate_dropout(dropout)
  coefficient <- switch(
    effect.measure,
    relative_risk = {
      .assert_probability(risk.experimental, "risk.experimental", open = TRUE)
      .assert_probability(risk.control, "risk.control", open = TRUE)
      (1 - risk.experimental) / (allocation.ratio * risk.experimental) +
        (1 - risk.control) / risk.control
    },
    risk_difference = {
      .assert_probability(risk.experimental, "risk.experimental", open = TRUE)
      .assert_probability(risk.control, "risk.control", open = TRUE)
      risk.experimental * (1 - risk.experimental) / allocation.ratio +
        risk.control * (1 - risk.control)
    },
    mean_difference = {
      .assert_scalar_numeric(sd.experimental, "sd.experimental", lower = 0,
                             lower_open = TRUE)
      .assert_scalar_numeric(sd.control, "sd.control", lower = 0,
                             lower_open = TRUE)
      sd.experimental^2 / allocation.ratio + sd.control^2
    }
  )
  n_control <- max(1L, ceiling(coefficient / required.variance))
  n_experimental <- max(1L, ceiling(allocation.ratio * n_control))
  achieved <- switch(
    effect.measure,
    relative_risk =
      (1 - risk.experimental) / (n_experimental * risk.experimental) +
        (1 - risk.control) / (n_control * risk.control),
    risk_difference =
      risk.experimental * (1 - risk.experimental) / n_experimental +
        risk.control * (1 - risk.control) / n_control,
    mean_difference =
      sd.experimental^2 / n_experimental + sd.control^2 / n_control
  )
  randomized_experimental <- ceiling(n_experimental / (1 - dropout[1]))
  randomized_control <- ceiling(n_control / (1 - dropout[2]))
  data.frame(
    Effect.measure = effect.measure,
    Target.variance = required.variance,
    Achieved.variance = achieved,
    Analyzable.total = n_experimental + n_control,
    Analyzable.experimental = n_experimental,
    Analyzable.control = n_control,
    Randomized.total = randomized_experimental + randomized_control,
    Randomized.experimental = randomized_experimental,
    Randomized.control = randomized_control,
    check.names = FALSE
  )
}

#' Design a binary-outcome non-inferiority trial using relative risks
#'
#' The general framework is applied on the log-relative-risk scale; sample
#' size is then obtained from the two-sample log-relative-risk variance.
#'
#' @param u.list,l1.list Equal-length vectors of framework parameters.
#' @param f.preserv Fraction of the active-control log relative risk to
#'   preserve, or `NULL`.
#' @param null.rr Maximum experimental-versus-placebo relative risk under the
#'   inferred-effect criterion, or `NULL`.
#' @param design.alternative.rr Experimental-versus-placebo relative risk under
#'   the design alternative.
#' @param hist.ac.rr Historical active-control-versus-placebo relative risk.
#' @param hist.ac.effect.se Historical-effect standard error on the log-RR
#'   scale.
#' @param placebo.risk Event risk under hypothetical placebo at the analysis
#'   time.
#' @param lambda0.for.design,target.on.unconditional.power,allocation.ratio
#'   Framework design inputs as in [ni_design()].
#' @param power Desired power.
#' @param sign.level One-sided significance level.
#' @param lambda0.sens.analysis Optional relative effect deviation for a power
#'   sensitivity analysis.
#' @param dropout Expected dropout proportion, one value or experimental and
#'   control values.
#' @return An object of class `c("ni.effect.design", "ni.design")`.
#' @examples
#' ni_design_rr(
#'   u.list = c(1, 0), l1.list = c(0, 0),
#'   design.alternative.rr = 0.50,
#'   hist.ac.rr = 0.65,
#'   hist.ac.effect.se = 0.12,
#'   placebo.risk = 0.20
#' )
#' @export
ni_design_rr <- function(u.list = c(1, 0), l1.list = c(0, 0),
                         f.preserv = 0.5, null.rr = NULL,
                         design.alternative.rr,
                         hist.ac.rr,
                         hist.ac.effect.se,
                         placebo.risk,
                         lambda0.for.design = 0,
                         target.on.unconditional.power = TRUE,
                         allocation.ratio = 1,
                         power = 0.9,
                         sign.level = 0.025,
                         lambda0.sens.analysis = NULL,
                         dropout = 0) {
  .assert_probability(design.alternative.rr, "design.alternative.rr",
                      open = TRUE)
  .assert_probability(hist.ac.rr, "hist.ac.rr", open = TRUE)
  if (!is.null(null.rr)) .assert_probability(null.rr, "null.rr", open = TRUE)
  .assert_probability(placebo.risk, "placebo.risk", open = TRUE)
  g_xp <- log(design.alternative.rr)
  g_cph <- log(hist.ac.rr)
  arm_parameters <- function(l0) {
    list(
      risk.experimental = placebo.risk * exp(g_xp),
      risk.control = placebo.risk * exp((1 + l0) * g_cph)
    )
  }
  criteria <- .effect_criteria(f.preserv, null.rr,
                               transform_null = log,
                               natural_label = "relative risk")
  .design_effect_engine(
    effect_measure = "relative_risk", u.list = u.list,
    l1.list = l1.list, criteria = criteria, g_xp = g_xp,
    g_cph = g_cph, se_cph = hist.ac.effect.se,
    l0 = lambda0.for.design,
    target_unconditional = target.on.unconditional.power,
    k = allocation.ratio, power = power, alpha = sign.level,
    l0_sensitivity = lambda0.sens.analysis, dropout = dropout,
    arm_parameters = arm_parameters,
    margin_transform = exp,
    active_effect_transform = function(x) exp((1 + x) * g_cph),
    effect_label = "log relative risk"
  )
}

#' Design a binary-outcome non-inferiority trial using risk differences
#'
#' The risk difference is already additive, so the general framework and the
#' two-sample risk-difference variance are applied on the same scale.
#'
#' @inheritParams ni_design_rr
#' @param null.rd Maximum experimental-versus-placebo risk difference for the
#'   inferred-effect criterion, or `NULL`.
#' @param design.alternative.rd Experimental-minus-placebo risk difference.
#' @param hist.ac.rd Historical active-control-minus-placebo risk difference.
#' @param hist.ac.effect.se Historical-effect standard error on the risk-
#'   difference scale.
#' @export
ni_design_rd <- function(u.list = c(1, 0), l1.list = c(0, 0),
                         f.preserv = 0.5, null.rd = NULL,
                         design.alternative.rd,
                         hist.ac.rd,
                         hist.ac.effect.se,
                         placebo.risk,
                         lambda0.for.design = 0,
                         target.on.unconditional.power = TRUE,
                         allocation.ratio = 1,
                         power = 0.9,
                         sign.level = 0.025,
                         lambda0.sens.analysis = NULL,
                         dropout = 0) {
  .assert_negative(design.alternative.rd, "design.alternative.rd")
  .assert_negative(hist.ac.rd, "hist.ac.rd")
  if (!is.null(null.rd)) .assert_nonpositive(null.rd, "null.rd")
  .assert_probability(placebo.risk, "placebo.risk", open = TRUE)
  arm_parameters <- function(l0) {
    list(
      risk.experimental = placebo.risk + design.alternative.rd,
      risk.control = placebo.risk + (1 + l0) * hist.ac.rd
    )
  }
  .validate_binary_arms(arm_parameters(lambda0.for.design))
  if (!is.null(lambda0.sens.analysis)) {
    .validate_binary_arms(arm_parameters(lambda0.sens.analysis))
  }
  criteria <- .effect_criteria(f.preserv, null.rd,
                               transform_null = identity,
                               natural_label = "risk difference")
  .design_effect_engine(
    effect_measure = "risk_difference", u.list = u.list,
    l1.list = l1.list, criteria = criteria,
    g_xp = design.alternative.rd, g_cph = hist.ac.rd,
    se_cph = hist.ac.effect.se, l0 = lambda0.for.design,
    target_unconditional = target.on.unconditional.power,
    k = allocation.ratio, power = power, alpha = sign.level,
    l0_sensitivity = lambda0.sens.analysis, dropout = dropout,
    arm_parameters = arm_parameters,
    margin_transform = identity,
    active_effect_transform = function(x) (1 + x) * hist.ac.rd,
    effect_label = "risk difference"
  )
}

#' Design a continuous-outcome non-inferiority trial using mean differences
#'
#' @inheritParams ni_design_rr
#' @param null.md Maximum experimental-versus-placebo mean difference for the
#'   inferred-effect criterion, or `NULL`.
#' @param design.alternative.md Experimental-minus-placebo mean difference.
#' @param hist.ac.md Historical active-control-minus-placebo mean difference.
#' @param hist.ac.effect.se Historical-effect standard error on the mean-
#'   difference scale.
#' @param sd.experimental,sd.control Planning standard deviations in the two
#'   randomized arms.
#' @export
ni_design_md <- function(u.list = c(1, 0), l1.list = c(0, 0),
                         f.preserv = 0.5, null.md = NULL,
                         design.alternative.md,
                         hist.ac.md,
                         hist.ac.effect.se,
                         sd.experimental,
                         sd.control,
                         lambda0.for.design = 0,
                         target.on.unconditional.power = TRUE,
                         allocation.ratio = 1,
                         power = 0.9,
                         sign.level = 0.025,
                         lambda0.sens.analysis = NULL,
                         dropout = 0) {
  .assert_negative(design.alternative.md, "design.alternative.md")
  .assert_negative(hist.ac.md, "hist.ac.md")
  if (!is.null(null.md)) .assert_nonpositive(null.md, "null.md")
  .assert_scalar_numeric(sd.experimental, "sd.experimental", lower = 0,
                         lower_open = TRUE)
  .assert_scalar_numeric(sd.control, "sd.control", lower = 0,
                         lower_open = TRUE)
  arm_parameters <- function(l0) {
    list(sd.experimental = sd.experimental, sd.control = sd.control)
  }
  criteria <- .effect_criteria(f.preserv, null.md,
                               transform_null = identity,
                               natural_label = "mean difference")
  .design_effect_engine(
    effect_measure = "mean_difference", u.list = u.list,
    l1.list = l1.list, criteria = criteria,
    g_xp = design.alternative.md, g_cph = hist.ac.md,
    se_cph = hist.ac.effect.se, l0 = lambda0.for.design,
    target_unconditional = target.on.unconditional.power,
    k = allocation.ratio, power = power, alpha = sign.level,
    l0_sensitivity = lambda0.sens.analysis, dropout = dropout,
    arm_parameters = arm_parameters,
    margin_transform = identity,
    active_effect_transform = function(x) (1 + x) * hist.ac.md,
    effect_label = "mean difference"
  )
}

.design_effect_engine <- function(effect_measure, u.list, l1.list, criteria,
                                  g_xp, g_cph, se_cph, l0,
                                  target_unconditional, k, power, alpha,
                                  l0_sensitivity, dropout, arm_parameters,
                                  margin_transform, active_effect_transform,
                                  effect_label) {
  .validate_design_vectors(u.list, l1.list)
  .validate_design_inputs(g_cph, se_cph, l0, g_xp, 0, 0, u.list[1],
                          l1.list[1], alpha, power)
  .assert_scalar_numeric(k, "allocation.ratio", lower = 0,
                         lower_open = TRUE)
  dropout <- .validate_dropout(dropout)
  if (!is.null(l0_sensitivity)) {
    .assert_scalar_numeric(l0_sensitivity, "lambda0.sens.analysis",
                           lower = -1, lower_open = TRUE)
  }
  if (length(criteria) == 0L) {
    stop("At least one non-inferiority criterion must be requested.",
         call. = FALSE)
  }
  arms <- arm_parameters(l0)
  if (effect_measure != "mean_difference") .validate_binary_arms(arms)
  methods <- mapply(method_name, u.list, l1.list,
                    MoreArgs = list(g_cph = g_cph, se_cph = se_cph),
                    USE.NAMES = FALSE)
  tables <- lapply(criteria, function(criterion) {
    .build_effect_table(
      effect_measure, methods, u.list, l1.list, g_xp, g_cph, se_cph,
      l0, criterion$f, criterion$delta0, target_unconditional, k,
      power, alpha, l0_sensitivity, dropout, arms, margin_transform,
      active_effect_transform
    )
  })
  type_power <- if (target_unconditional) "unconditional" else "conditional"
  specifications <- list(
    Approach = paste0("Design targeting ", round(100 * power, 2), "% ",
                      type_power, " power"),
    `Effect measure` = effect_label,
    `Planning active-control effect` =
      format(active_effect_transform(l0), digits = 5),
    `Allocation ratio` = paste0(k, " experimental:1 control")
  )
  if (!is.null(l0_sensitivity) && l0_sensitivity != l0) {
    specifications[["Sensitivity active-control effect"]] <-
      format(active_effect_transform(l0_sensitivity), digits = 5)
  }
  out <- c(list(Specifications = specifications), tables)
  class(out) <- c("ni.effect.design", "ni.design")
  invisible(out)
}

.build_effect_table <- function(effect_measure, methods, u, l1, g_xp,
                                g_cph, se_cph, l0, f, delta0,
                                target_unconditional, k, power, alpha,
                                l0_sensitivity, dropout, arms,
                                margin_transform, active_effect_transform) {
  solver <- if (target_unconditional) solve_v_xs_uncond else solve_v_xs_cond
  variance <- mapply(
    solver, u = u, l1 = l1,
    MoreArgs = list(g_cph = g_cph, se_cph = se_cph, l0 = l0,
                    g_xp = g_xp, f = f, delta0 = delta0,
                    alpha = alpha, power = power),
    USE.NAMES = FALSE
  )
  make_size <- function(target) {
    if (!is.finite(target) || target <= 0) {
      return(data.frame(
        Target.variance = NA_real_, Achieved.variance = NA_real_,
        Analyzable.total = NA_integer_, Analyzable.experimental = NA_integer_,
        Analyzable.control = NA_integer_, Randomized.total = NA_integer_,
        Randomized.experimental = NA_integer_, Randomized.control = NA_integer_
      ))
    }
    args <- c(
      list(required.variance = target, effect.measure = effect_measure,
           allocation.ratio = k, dropout = dropout),
      arms
    )
    size <- do.call(ni_sample_size, args)
    size[names(size) != "Effect.measure"]
  }
  sizes <- do.call(rbind, lapply(variance, make_size))
  margin <- vapply(seq_along(variance), function(i) {
    if (!is.finite(variance[i]) || variance[i] <= 0) return(NA_real_)
    margin_transform(.margin_linear(
      g_cph, se_cph, sqrt(variance[i]), f, delta0, u[i], l1[i], alpha
    ))
  }, numeric(1))
  tolerance <- vapply(seq_along(variance), function(i) {
    if (!is.finite(variance[i]) || variance[i] <= 0) return(NA_real_)
    min_l0(g_cph, se_cph, sqrt(variance[i]), f, u[i], l1[i], alpha)
  }, numeric(1))
  unconditional <- vapply(seq_along(variance), function(i) {
    if (!is.finite(variance[i]) || variance[i] <= 0) return(NA_real_)
    uncond_power(g_cph, se_cph, l0, g_xp, sqrt(variance[i]), f,
                 delta0, u[i], l1[i], alpha)
  }, numeric(1))
  out <- data.frame(
    Method = methods,
    `NI margin` = signif(margin, 5),
    `Target variance` = signif(variance, 6),
    `Achieved variance` = signif(sizes$Achieved.variance, 6),
    `Sample size` = sizes$Randomized.total,
    Exp.arm = sizes$Randomized.experimental,
    Ctr.arm = sizes$Randomized.control,
    `Analyzable sample size` = sizes$Analyzable.total,
    `AC effect at T1E limit` = signif(active_effect_transform(tolerance), 5),
    U.power = round(unconditional, 3),
    Status = ifelse(is.finite(variance), "OK", "Not attainable"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!is.null(l0_sensitivity) && l0_sensitivity != l0) {
    sensitivity <- vapply(seq_along(variance), function(i) {
      if (!is.finite(variance[i]) || variance[i] <= 0) return(NA_real_)
      uncond_power(g_cph, se_cph, l0_sensitivity, g_xp,
                   sqrt(variance[i]), f, delta0, u[i], l1[i], alpha)
    }, numeric(1))
    out[["U.power (SA)"]] <- round(sensitivity, 3)
  }
  out
}

.effect_criteria <- function(f, null, transform_null, natural_label) {
  out <- list()
  if (!is.null(f)) {
    .assert_scalar_numeric(f, "f.preserv", lower = 0, upper = 1)
    out[[paste0("NI criterion: Preserving ", round(100 * f, 2),
                "% of active control effect")]] <- list(f = f, delta0 = 0)
  }
  if (!is.null(null)) {
    out[[paste0("NI criterion: Inferred ", natural_label, " of ",
                format(null, trim = TRUE))]] <-
      list(f = 0, delta0 = transform_null(null))
  }
  out
}

.validate_dropout <- function(dropout) {
  if (!is.numeric(dropout) || length(dropout) < 1L || length(dropout) > 2L ||
      anyNA(dropout) || any(!is.finite(dropout)) ||
      any(dropout < 0 | dropout >= 1)) {
    stop("`dropout` must contain one or two proportions in [0, 1).",
         call. = FALSE)
  }
  if (length(dropout) == 1L) rep(dropout, 2L) else dropout
}

.validate_binary_arms <- function(arms) {
  .assert_probability(arms$risk.experimental, "implied experimental risk",
                      open = TRUE)
  .assert_probability(arms$risk.control, "implied control risk", open = TRUE)
  invisible(NULL)
}

.assert_negative <- function(x, name) {
  .assert_scalar_numeric(x, name, upper = 0, upper_open = TRUE)
}

.assert_nonpositive <- function(x, name) {
  .assert_scalar_numeric(x, name, upper = 0)
}
