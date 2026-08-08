test_that("mean-difference formula is exact before integer rounding", {
  result <- ni_sample_size(
    required.variance = 0.04,
    effect.measure = "mean_difference",
    sd.experimental = 2,
    sd.control = 2
  )
  expect_equal(result$Analyzable.experimental, 200)
  expect_equal(result$Analyzable.control, 200)
  expect_equal(result$Achieved.variance, 0.04)
})

test_that("binary formulas achieve the requested variance", {
  for (measure in c("relative_risk", "risk_difference")) {
    result <- ni_sample_size(
      required.variance = 0.01,
      effect.measure = measure,
      allocation.ratio = 2,
      dropout = c(0.10, 0.20),
      risk.experimental = 0.10,
      risk.control = 0.15
    )
    expect_lte(result$Achieved.variance, result$Target.variance)
    expect_gte(result$Randomized.experimental,
               result$Analyzable.experimental)
    expect_gte(result$Randomized.control, result$Analyzable.control)
    expect_equal(
      result$Randomized.experimental,
      ceiling(result$Analyzable.experimental / 0.9)
    )
  }
})

test_that("new design wrappers return finite, variance-valid sample sizes", {
  designs <- list(
    ni_design_rr(
      design.alternative.rr = 0.50, hist.ac.rr = 0.65,
      hist.ac.effect.se = 0.12, placebo.risk = 0.20
    ),
    ni_design_rd(
      design.alternative.rd = -0.10, hist.ac.rd = -0.07,
      hist.ac.effect.se = 0.02, placebo.risk = 0.20
    ),
    ni_design_md(
      design.alternative.md = -2, hist.ac.md = -1.5,
      hist.ac.effect.se = 0.4,
      sd.experimental = 5, sd.control = 5
    )
  )
  for (design in designs) {
    expect_s3_class(design, "ni.effect.design")
    table <- design[[2]]
    expect_true(all(table$Status == "OK"))
    expect_true(all(table$`Achieved variance` <=
                      table$`Target variance` + 1e-6))
    expect_true(all(table$`Sample size` > 0))
    expect_equal(table$U.power, rep(0.9, nrow(table)), tolerance = 0.001)
  }
})

test_that("invalid planning values are rejected", {
  expect_error(
    ni_sample_size(0.01, "relative_risk",
                   risk.experimental = 0, risk.control = 0.2),
    "risk.experimental"
  )
  expect_error(
    ni_design_rd(
      design.alternative.rd = -0.3, hist.ac.rd = -0.1,
      hist.ac.effect.se = 0.02, placebo.risk = 0.2
    ),
    "implied experimental risk"
  )
  expect_error(
    ni_sample_size(0.01, "mean_difference", dropout = 1,
                   sd.experimental = 1, sd.control = 1),
    "dropout"
  )
})
