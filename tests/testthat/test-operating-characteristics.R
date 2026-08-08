test_that("conditional operating characteristics are more optimistic", {
  g_cph <- log(0.072)
  se_cph <- 0.61
  g_xp <- log(0.05)
  se_xs <- 0.40

  expect_gt(
    cond_power(g_cph, se_cph, l0 = 0, g_xp, se_xs,
               f = 0.5, u = 1, l1 = 0),
    uncond_power(g_cph, se_cph, l0 = 0, g_xp, se_xs,
                 f = 0.5, u = 1, l1 = 0)
  )
  expect_lt(
    cond_t1e(g_cph, se_cph, l0 = 0, se_xs,
             f = 0.5, u = 1, l1 = 0),
    uncond_t1e(g_cph, se_cph, l0 = 0, se_xs,
               f = 0.5, u = 1, l1 = 0)
  )
})

test_that("maximum unconditional power uses the manuscript variance", {
  g_cph <- log(0.072)
  se_cph <- 0.61
  g_xp <- log(0.05)
  f <- 0.5
  l1 <- 1.96 * se_cph / g_cph
  a <- (1 - (1 - f) * (1 + l1)) * g_cph - g_xp

  expected <- stats::pnorm(a / ((1 - f) * se_cph))
  actual <- max_uncond_power(
    g_cph, se_cph, l0 = 0, g_xp, f = f,
    delta0 = 0, u = 0, l1 = l1
  )

  expect_equal(actual, expected)
  expect_equal(
    actual,
    max_up_unrest(g_cph, se_cph, l0 = 0, g_xp, f = f,
                  delta0 = 0, u = 0, l1 = l1)
  )
})
