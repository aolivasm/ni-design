# Core normal-approximation formulas from Olivas-Martinez et al. (2026).

#' Operating characteristics for a non-inferiority test
#'
#' These functions implement equations (7)--(10) and related results in
#' Olivas-Martinez, Gao, and Janes (2026). Effects must be expressed on an
#' additive scale for which smaller values indicate benefit.
#'
#' @param g_cph Historical active-control effect estimate on the additive
#'   analysis scale.
#' @param se_cph Standard error of `g_cph`.
#' @param l0 True relative effect deviation from the historical effect.
#' @param se_xs Standard error of the experimental-versus-control estimate.
#' @param f Fraction of the active-control effect to preserve. Use `0` for an
#'   inferred-efficacy criterion.
#' @param u Unifying parameter: `0` for fixed-margin, `1` for synthesis, and
#'   `1 / (1 + l1)` for Odem-Davis.
#' @param l1 Assumed relative effect deviation.
#' @param alpha One-sided significance level.
#' @param g_xp Experimental-versus-placebo effect under the design alternative.
#' @param delta0 Fixed null efficacy on the additive scale; use `0` for a
#'   preservation-of-effect criterion.
#' @return A numeric probability, except `min_l0()`, which returns the smallest
#'   `l0` for which unconditional type-I error is controlled.
#' @details
#' In the maximum-unconditional-power result, the historical variance is
#' `se_cph^2` for fixed-margin methods (`u = 0`) and
#' `(1 + l1)^2 * se_cph^2` otherwise, as specified by equation (4) of the
#' manuscript.
#' @name operating_characteristics
NULL

#' @rdname operating_characteristics
#' @export
uncond_t1e <- function(g_cph, se_cph, l0 = 0, se_xs, f = 0.5, u = 1,
                       l1 = 0, alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(se_xs, "se_xs", lower = 0, lower_open = TRUE)
  bias <- (1 - f) * (l0 - l1) * g_cph
  v_aux <- se_xs^2 + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2
  numer <- bias - sqrt(v_aux) * stats::qnorm(1 - alpha)
  denom <- sqrt(se_xs^2 + (1 - f)^2 * (1 + (u > 0) * l1)^2 * se_cph^2)
  stats::pnorm(numer / denom)
}

#' @rdname operating_characteristics
#' @export
cond_t1e <- function(g_cph, se_cph, l0 = 0, se_xs, f = 0.5, u = 1,
                     l1 = 0, alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(se_xs, "se_xs", lower = 0, lower_open = TRUE)
  bias <- (1 - f) * (l0 - l1) * g_cph
  v_aux <- se_xs^2 + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2
  numer <- bias - sqrt(v_aux) * stats::qnorm(1 - alpha)
  stats::pnorm(numer / se_xs)
}

#' @rdname operating_characteristics
#' @export
min_l0 <- function(g_cph, se_cph, se_xs, f = 0.5, u = 1, l1 = 0,
                   alpha = 0.025) {
  .validate_core(g_cph, se_cph, 0, f, u, l1, alpha)
  .assert_scalar_numeric(se_xs, "se_xs", lower = 0, lower_open = TRUE)
  if (f == 1) {
    return(0)
  }
  v_aux <- se_xs^2 + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2
  numer <- sqrt(v_aux) -
    sqrt(se_xs^2 + (1 - f)^2 * (1 + (u > 0) * l1)^2 * se_cph^2)
  l1 + numer * stats::qnorm(1 - alpha) / (1 - f) / g_cph
}

#' @rdname operating_characteristics
#' @export
uncond_power <- function(g_cph, se_cph, l0 = 0, g_xp, se_xs, f = 0.5,
                         delta0 = 0, u = 1, l1 = 0, alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(g_xp, "g_xp")
  .assert_scalar_numeric(se_xs, "se_xs", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(delta0, "delta0")
  v_aux <- se_xs^2 + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2
  numer <- delta0 + ((1 + l0) - (1 - f) * (1 + l1)) * g_cph -
    g_xp - sqrt(v_aux) * stats::qnorm(1 - alpha)
  denom <- sqrt(se_xs^2 + (1 - f)^2 *
                  (1 + (u > 0) * l1)^2 * se_cph^2)
  stats::pnorm(numer / denom)
}

#' @rdname operating_characteristics
#' @export
cond_power <- function(g_cph, se_cph, l0 = 0, g_xp, se_xs, f = 0.5,
                       delta0 = 0, u = 1, l1 = 0, alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(g_xp, "g_xp")
  .assert_scalar_numeric(se_xs, "se_xs", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(delta0, "delta0")
  v_aux <- se_xs^2 + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2
  numer <- delta0 + ((1 + l0) - (1 - f) * (1 + l1)) * g_cph -
    g_xp - sqrt(v_aux) * stats::qnorm(1 - alpha)
  stats::pnorm(numer / se_xs)
}

#' @rdname operating_characteristics
#' @export
max_uncond_power <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                             delta0 = 0, u = 1, l1 = 0,
                             alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(g_xp, "g_xp")
  .assert_scalar_numeric(delta0, "delta0")
  a <- delta0 + (1 + l0 - (1 - f) * (1 + l1)) * g_cph - g_xp
  if (f == 1) {
    if (a > 0) return(1)
    return(if (a == 0) alpha else NA_real_)
  }
  if (a - u * stats::qnorm(1 - alpha) * (1 + l1) *
      (1 - f) * se_cph < 0) {
    return(NA_real_)
  }
  historical_scale <- 1 + (u > 0) * l1
  stats::pnorm(-u * stats::qnorm(1 - alpha) +
                 a / (1 - f) / historical_scale / se_cph)
}

#' @rdname operating_characteristics
#' @export
max_up_unrest <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                          delta0 = 0, u = 1, l1 = 0, alpha = 0.025) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(g_xp, "g_xp")
  .assert_scalar_numeric(delta0, "delta0")
  a <- delta0 + ((1 + l0) - (1 - f) * (1 + l1)) * g_cph - g_xp
  if (f == 1) {
    if (a > 0) return(1)
    return(if (a == 0) alpha else 0)
  }
  stats::pnorm(-u * stats::qnorm(1 - alpha) +
                 a / (1 - f) / (1 + (u > 0) * l1) / se_cph)
}

#' Label a framework method
#'
#' @param u,l1 Framework parameters.
#' @param g_cph,se_cph Historical effect estimate and standard error.
#' @return A character label.
#' @export
method_name <- function(u, l1, g_cph, se_cph) {
  # The published implementation and Appendix B use 1.96 exactly.
  l1_95 <- 1.96 * se_cph / g_cph
  close <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-8))
  if (close(u, 0)) {
    if (close(l1, l1_95)) {
      return("95-95 method")
    }
    if (close(l1, 0)) {
      return("0-95 method")
    }
    return(paste0("FM method with lambda1=", round(l1, 4)))
  }
  if (close(u, 1)) {
    if (close(l1, 0)) {
      return("Traditional SM")
    }
    return(paste0("BA-SM, \u03bb1=", round(l1 * 100, 2), "%"))
  }
  if (l1 != 0 && close(u, 1 / (1 + l1))) {
    return(paste0("OD, \u03bb1=", round(l1 * 100, 2), "%"))
  }
  paste0("Method with (u, lambda1)= (", round(u, 4), ", ",
         round(l1, 4), ")")
}

#' Solve the required active-trial variance
#'
#' Solves equation (7) or (9) for the variance of the
#' experimental-versus-control effect estimator.
#'
#' @inheritParams operating_characteristics
#' @param power Desired power.
#' @return The required variance, or `NA_real_` when the design alternative
#'   cannot attain the requested power.
#' @name solve_variance
NULL

#' @rdname solve_variance
#' @export
solve_v_xs_uncond <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                              delta0 = 0, u = 1, l1 = 0,
                              alpha = 0.025, power = 0.9) {
  .validate_design_inputs(g_cph, se_cph, l0, g_xp, f, delta0, u, l1,
                          alpha, power)
  max_pow <- max_uncond_power(g_cph, se_cph, l0, g_xp, f, delta0,
                              u, l1, alpha)
  if (is.na(max_pow) || max_pow < power) {
    return(NA_real_)
  }
  if (f == 1) {
    return((delta0 - g_xp + (1 + l0) * g_cph)^2 /
             (stats::qnorm(1 - alpha) + stats::qnorm(power))^2)
  }
  c_hist <- (1 - f) * (1 + (u > 0) * l1) * se_cph
  objective <- function(x) {
    v_aux <- sqrt(x + u^2 * (1 - f)^2 * (1 + l1)^2 * se_cph^2)
    stats::qnorm(power) +
      (v_aux * stats::qnorm(1 - alpha) - delta0 +
         (1 - f) * (1 + l1) * g_cph + g_xp -
         (1 + l0) * g_cph) / sqrt(x + c_hist^2)
  }
  .solve_positive_variance(objective)
}

#' @rdname solve_variance
#' @export
solve_v_xs_cond <- function(g_cph, se_cph, l0 = 0, g_xp, f = 0.5,
                            delta0 = 0, u = 1, l1 = 0,
                            alpha = 0.025, power = 0.9) {
  .validate_design_inputs(g_cph, se_cph, l0, g_xp, f, delta0, u, l1,
                          alpha, power)
  detectable <- delta0 > (1 - f) * (1 + l1) *
    (g_cph + u * se_cph * stats::qnorm(1 - alpha)) +
    g_xp - (1 + l0) * g_cph
  if (!detectable) {
    return(NA_real_)
  }
  if (f == 1) {
    return((delta0 - g_xp + (1 + l0) * g_cph)^2 /
             (stats::qnorm(1 - alpha) + stats::qnorm(power))^2)
  }
  c_hist <- (1 - f) * (1 + l1) * se_cph
  objective <- function(x) {
    v_aux <- sqrt(x + u^2 * c_hist^2)
    stats::qnorm(power) +
      (v_aux * stats::qnorm(1 - alpha) - delta0 +
         (1 - f) * (1 + l1) * g_cph + g_xp -
         (1 + l0) * g_cph) / sqrt(x)
  }
  .solve_positive_variance(objective)
}

.solve_positive_variance <- function(objective) {
  result <- suppressWarnings(nleqslv::nleqslv(0.001, objective))
  value <- as.numeric(result$x)
  if (length(value) != 1L || !is.finite(value) || value <= 0 ||
      abs(objective(value)) > 1e-5) {
    return(NA_real_)
  }
  value
}

.margin_linear <- function(g_cph, se_cph, se_xs, f = 0.5, delta0 = 0,
                           u = 1, l1 = 0, alpha = 0.025) {
  v_aux <- se_xs^2 + u^2 * (1 + l1)^2 * (1 - f)^2 * se_cph^2
  delta0 - (1 - f) * (1 + l1) * g_cph -
    stats::qnorm(1 - alpha) * (sqrt(v_aux) - se_xs)
}

#' Compute the non-inferiority margin on a ratio scale
#'
#' @inheritParams operating_characteristics
#' @param d0 Fixed null efficacy on the log scale.
#' @return The multiplicative non-inferiority margin, rounded to two decimals
#'   for compatibility with the manuscript output.
#' @export
ni_margin <- function(g_cph, se_cph, se_xs, f = 0.5, d0 = 0, u = 1,
                      l1 = 0, alpha = 0.025) {
  round(exp(.margin_linear(g_cph, se_cph, se_xs, f, d0, u, l1,
                           alpha)), 2)
}

.validate_core <- function(g_cph, se_cph, l0, f, u, l1, alpha) {
  .assert_scalar_numeric(g_cph, "g_cph")
  if (g_cph == 0) stop("`g_cph` must not be zero.", call. = FALSE)
  .assert_scalar_numeric(se_cph, "se_cph", lower = 0, lower_open = TRUE)
  .assert_scalar_numeric(l0, "l0")
  .assert_scalar_numeric(f, "f", lower = 0, upper = 1)
  .assert_scalar_numeric(u, "u", lower = 0)
  .assert_scalar_numeric(l1, "l1", lower = -1, lower_open = TRUE)
  .assert_scalar_numeric(alpha, "alpha", lower = 0, upper = 0.5,
                         lower_open = TRUE, upper_open = TRUE)
  invisible(NULL)
}

.validate_design_inputs <- function(g_cph, se_cph, l0, g_xp, f, delta0,
                                    u, l1, alpha, power) {
  .validate_core(g_cph, se_cph, l0, f, u, l1, alpha)
  .assert_scalar_numeric(g_xp, "g_xp")
  .assert_scalar_numeric(delta0, "delta0")
  .assert_scalar_numeric(power, "power", lower = 0.5, upper = 1,
                         lower_open = TRUE, upper_open = TRUE)
  invisible(NULL)
}

.assert_scalar_numeric <- function(x, name, lower = -Inf, upper = Inf,
                                   lower_open = FALSE, upper_open = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be one finite numeric value.", name), call. = FALSE)
  }
  lower_bad <- if (lower_open) x <= lower else x < lower
  upper_bad <- if (upper_open) x >= upper else x > upper
  if (lower_bad || upper_bad) {
    interval <- paste0(if (lower_open) "(" else "[", lower, ", ", upper,
                       if (upper_open) ")" else "]")
    stop(sprintf("`%s` must lie in %s.", name, interval), call. = FALSE)
  }
  invisible(NULL)
}
