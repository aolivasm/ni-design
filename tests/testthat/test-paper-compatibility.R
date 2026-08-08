paper_design <- function() {
  ni_design(
    u.list = c(1, 1, 1 / (1 - 0.23), 0, 0),
    l1.list = c(0, -0.23, -0.23, 1.96 * 0.61 / log(1 - 0.928), 0),
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
}

test_that("Appendix B prevention-efficacy output is preserved", {
  design <- paper_design()
  expect_s3_class(design, "ni.design")
  expect_equal(
    design$Specifications,
    list(
      Approach = paste0(
        "Design approach targeting 90% conditional power and assuming ",
        "an active control efficacy of 92.8%"
      ),
      `Sensitivity analysis` = paste0(
        "Sensitivity analysis (SA) assumes an active control efficacy ",
        "of 94.7%"
      )
    )
  )

  preservation <- design[[2]]
  expect_equal(
    preservation$Method,
    c("Traditional SM", "BA-SM, \u03bb1=-23%", "OD, \u03bb1=-23%",
      "95-95 method", "0-95 method")
  )
  expect_equal(preservation$`NI margin`, c(3.12, 2.42, 2.21, 2.05, 3.73))
  expect_equal(preservation$RNE, c(19, 27, 32, 37, 15))
  expect_equal(preservation$Exp, c(8, 11, 13, 15, 6))
  expect_equal(preservation$Ctr, c(11, 16, 19, 22, 9))
  expect_equal(preservation$`Sample size`,
               c(5766, 8008, 9510, 11012, 4504))
  expect_equal(preservation$CNC, c(0.928, 0.868, 0.844, 0.850, 0.948))
  expect_equal(preservation$U.power, c(0.86, 0.86, 0.86, 0.83, 0.87))
  expect_equal(preservation$`U.power (SA)`,
               c(0.69, 0.65, 0.63, 0.60, 0.72))

  inferred <- design[[3]]
  expect_equal(inferred$`NI margin`, c(6.13, 3.72, 2.88, 2.94, 9.72))
  expect_equal(inferred$RNE, c(9, 15, 22, 21, 7))
  expect_equal(inferred$`Sample size`, c(2882, 4504, 6506, 6486, 2162))
  expect_equal(inferred$CNC, c(0.928, 0.868, 0.837, 0.870, 0.952))
  expect_equal(inferred$U.power, c(0.83, 0.83, 0.81, 0.78, 0.85))
  expect_equal(inferred$`U.power (SA)`,
               c(0.73, 0.69, 0.65, 0.63, 0.76))
})

test_that("non-default alpha and power are propagated consistently", {
  design <- ni_design(
    u.list = 0, l1.list = 0,
    f.preserv = 0.5, null.pe = NULL,
    design.alternative.pe = 0.85,
    hist.ac.pe = 0.55,
    hist.ac.effect.se = 0.25,
    power = 0.8, sign.level = 0.05
  )
  table <- design[[2]]
  variance <- solve_v_xs_uncond(
    g_cph = log(0.45), se_cph = 0.25,
    g_xp = log(0.15), f = 0.5,
    u = 0, l1 = 0, alpha = 0.05, power = 0.8
  )
  expected_events <- rne_uncond(
    g_cph = log(0.45), se_cph = 0.25,
    g_xp = log(0.15), f = 0.5,
    u = 0, l1 = 0, alpha = 0.05, power = 0.8
  )
  expect_equal(table$RNE, as.numeric(expected_events[1]))
  expect_equal(table$`NI margin`,
               ni_margin(log(0.45), 0.25, sqrt(variance),
                         f = 0.5, u = 0, l1 = 0, alpha = 0.05))
})

test_that("unattainable paper-scale designs return missing results cleanly", {
  expect_no_error(
    design <- ni_design(
      u.list = 1, l1.list = 0,
      f.preserv = 0.5, null.pe = NULL,
      design.alternative.pe = 0.20,
      hist.ac.pe = 0.80,
      hist.ac.effect.se = 0.50,
      power = 0.90
    )
  )
  expect_true(is.na(design[[2]]$`Sample size`))
  expect_true(is.na(design[[2]]$U.power))
})
