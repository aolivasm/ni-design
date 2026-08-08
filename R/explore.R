#' Explore maximum unconditional power
#'
#' Evaluates maximum attainable unconditional power across experimental
#' prevention-efficacy alternatives for one or both manuscript criteria.
#'
#' @param u.list,l1.list Equal-length vectors of framework parameters.
#' @param f.preserv Fraction of active-control effect to preserve, or `NULL`.
#' @param null.pe Minimum inferred prevention efficacy, or `NULL`.
#' @param design.alt.pe Vector of experimental prevention efficacies.
#' @param hist.ac.pe Historical active-control prevention efficacy.
#' @param hist.ac.effect.se Historical log hazard-ratio standard error.
#' @param lambda0.for.design Relative effect deviation assumed for planning.
#' @param sign.level One-sided significance level.
#' @param opt.power Reference power shown as a horizontal line.
#' @param legend.position A ggplot2 legend position.
#' @return A list containing a data frame and plot for each requested
#'   criterion.
#' @export
explore_max_uncond_power <- function(
    u.list = c(1, 1, 0, 0),
    l1.list = c(0, -0.23, -0.45, 0),
    f.preserv = 0.5,
    null.pe = 0.3,
    design.alt.pe = seq(0.65, 0.98, by = 0.01),
    hist.ac.pe = 0.928,
    hist.ac.effect.se = 0.61,
    lambda0.for.design = 0,
    sign.level = 0.025,
    opt.power = 0.9,
    legend.position = "bottom") {
  .validate_design_vectors(u.list, l1.list)
  if (!is.numeric(design.alt.pe) || length(design.alt.pe) == 0L ||
      anyNA(design.alt.pe) || any(design.alt.pe <= 0 | design.alt.pe >= 1)) {
    stop("`design.alt.pe` must contain prevention efficacies in (0, 1).",
         call. = FALSE)
  }
  .assert_probability(hist.ac.pe, "hist.ac.pe", open = TRUE)
  .assert_scalar_numeric(hist.ac.effect.se, "hist.ac.effect.se", lower = 0,
                         lower_open = TRUE)
  .assert_scalar_numeric(opt.power, "opt.power", lower = 0, upper = 1)
  if (is.null(f.preserv) && is.null(null.pe)) {
    stop("At least one criterion must be requested.", call. = FALSE)
  }
  g_xp <- log(1 - design.alt.pe)
  g_cph <- log(1 - hist.ac.pe)
  methods <- mapply(method_name, u.list, l1.list,
                    MoreArgs = list(g_cph = g_cph,
                                    se_cph = hist.ac.effect.se),
                    USE.NAMES = FALSE)
  build <- function(f, delta0, title) {
    grid <- expand.grid(i = seq_along(g_xp), j = seq_along(u.list))
    power_values <- mapply(
      function(i, j) max_up_unrest(
        g_cph = g_cph, se_cph = hist.ac.effect.se,
        l0 = lambda0.for.design, g_xp = g_xp[i], f = f,
        delta0 = delta0, u = u.list[j], l1 = l1.list[j],
        alpha = sign.level
      ),
      grid$i, grid$j
    )
    table <- data.frame(
      Exp_pe = round(100 * rep(design.alt.pe, times = length(u.list)), 2),
      Method = rep(methods, each = length(design.alt.pe)),
      Max_up = power_values
    )
    plot <- ggplot2::ggplot(
      table,
      ggplot2::aes(x = Exp_pe, y = Max_up, color = Method,
                   linetype = Method)
    ) +
      ggplot2::geom_hline(yintercept = opt.power, linetype = 2) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::labs(
        x = "Efficacy of experimental intervention (%)",
        y = "Maximum unconditional power",
        title = title
      ) +
      ggplot2::theme(legend.position = legend.position)
    list(table = table, plot = plot)
  }
  out <- list()
  if (!is.null(f.preserv)) {
    item <- build(
      f.preserv, 0,
      paste0("Criterion: Preserving ", round(100 * f.preserv, 2),
             "% of active control effect")
    )
    out$table.preserv.eff <- item$table
    out$plot.preserv.eff <- item$plot
  }
  if (!is.null(null.pe)) {
    .assert_probability(null.pe, "null.pe", open = TRUE)
    item <- build(
      0, log(1 - null.pe),
      paste0("Criterion: Inferred efficacy of ",
             round(100 * null.pe, 2), "%.")
    )
    out$table.inf.eff <- item$table
    out$plot.inf.eff <- item$plot
  }
  out
}

utils::globalVariables(c("Exp_pe", "Max_up", "Method"))
