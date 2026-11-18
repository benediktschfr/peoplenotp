test_that("joint_individual_certainty() computes joint summaries, Makowski tiers, and inverted evidence matrix", {
  # 6 individuals across 200 draws each for 2 outcomes:
  # outcome 1 (strength): better_is = "higher", MID = 4.0
  # outcome 2 (anxiety): better_is = "lower", MID = 2.0
  df_draws <- tibble::tibble(
    id = rep(paste0("P", 1:6), each = 200),
    .draw = rep(1:200, times = 6),
    strength = c(
      rep(-6, 200), # P1: 100% det on both -> Definite Deterioration
      c(rep(-6, 197), rep(0, 3)), # P2: 98.5% det on both -> Probable Deterioration
      c(rep(6, 100), rep(0, 100)), # P3: 50% imp -> Undecided
      rep(0, 200), # P4: 100% in joint ROPE -> Practically Equivalent
      c(rep(6, 197), rep(0, 3)), # P5: 98.5% imp on both -> Probable Improvement
      rep(6, 200) # P6: 100% imp on both -> Definite Improvement
    ),
    anxiety = c(
      rep(4, 200), # P1: det (> 2)
      c(rep(4, 197), rep(0, 3)), # P2: det (> 2)
      c(rep(-4, 100), rep(0, 100)), # P3: 50% imp
      rep(0, 200), # P4: in ROPE
      c(rep(-4, 197), rep(0, 3)), # P5: imp (< -2)
      rep(-4, 200) # P6: imp (< -2)
    )
  )

  res <- joint_individual_certainty(
    data = df_draws,
    id = id,
    values = c(strength, anxiety),
    mid = list(strength = c(3.0, 4.0), anxiety = c(1.5, 2.0)),
    better_is = c("higher", "lower"),
    include_marginals = TRUE
  )

  expect_s3_class(res, "pnp_joint_individual_certainty")
  expect_named(res, c("individual_summary", "matrix_data", "meta"))

  ind_sum <- res$individual_summary
  expect_s3_class(ind_sum, "tbl_df")
  expect_identical(nrow(ind_sum), 12L) # 6 individuals * 2 MID steps

  expected_levels <- c(
    "Definite Deterioration (>99%)",
    "Probable Deterioration (>97.5%)",
    "Undecided",
    "Practically Equivalent (>97.5% ROPE)",
    "Probable Improvement (>97.5%)",
    "Definite Improvement (>99%)"
  )
  sub_step2 <- ind_sum[ind_sum$mid == "4.0 / 2.0", ]
  expect_identical(as.character(sub_step2$makowski_tier), expected_levels)
  expect_true("improvement_pct_strength" %in% names(ind_sum))
  expect_true("improvement_pct_anxiety" %in% names(ind_sum))

  mat <- res$matrix_data
  expect_identical(nrow(mat), 6L)
  expect_named(mat, c("makowski_tier", "MID = 3.0 / 1.5", "MID = 4.0 / 2.0"))
  expect_identical(mat$`MID = 4.0 / 2.0`, rep("1 (16.7%)", 6))
})

test_that("joint_individual_certainty() supports >= 3 outcomes, .by grouping, and single-case data", {
  df_three <- tibble::tibble(
    id = rep(c("P1", "P2"), each = 100),
    group = rep(c("Control", "Treatment"), each = 100),
    .draw = rep(1:100, times = 2),
    y1 = c(rep(0, 100), rep(-5, 100)),
    y2 = c(rep(0, 100), rep(-5, 100)),
    y3 = c(rep(0, 100), rep(5, 100))
  )

  res3 <- joint_individual_certainty(
    data = df_three,
    id = id,
    values = c(y1, y2, y3),
    mid = c(2, 2, 2),
    .by = group,
    better_is = c("lower", "lower", "higher")
  )

  expect_s3_class(res3, "pnp_joint_individual_certainty")
  expect_identical(res3$meta$value_cols, c("y1", "y2", "y3"))
  expect_equal(res3$individual_summary$joint_improvement_pct, c(0, 1))
  expect_equal(res3$individual_summary$joint_rope_pct, c(1, 0))

  out_print <- paste(capture.output(print(res3)), collapse = "\n")
  expect_match(out_print, "# Joint Clinical Significance Evidence Matrix")
  expect_match(out_print, "group = Control")
  expect_match(out_print, "group = Treatment")

  # Plot with outcomes_select for 3 outcomes
  p_pair <- plot(res3, outcomes_select = c("y1", "y3"))
  expect_s3_class(p_pair, "ggplot")
  expect_match(p_pair$labels$x, "^y1")
  expect_match(p_pair$labels$y, "^y3")

  # Single-case without id column
  df_single <- tibble::tibble(
    .draw = 1:50,
    y1 = rep(-4, 50),
    y2 = rep(-4, 50)
  )
  res_single <- joint_individual_certainty(
    df_single,
    values = c(y1, y2),
    mid = 2,
    better_is = "lower"
  )
  expect_null(res_single$meta$id_col)
  expect_identical(nrow(res_single$individual_summary), 1L)
  expect_equal(res_single$individual_summary$joint_improvement_pct, 1)
})

test_that("plot.pnp_joint_individual_certainty() builds 2D joint scatter and sensitivity heatmaps", {
  df_draws <- tibble::tibble(
    id = rep(paste0("P", 1:4), each = 100),
    condition = rep(c("A", "A", "B", "B"), each = 100),
    .draw = rep(1:100, times = 4),
    y1 = c(rep(-5, 100), rep(-3, 100), rep(5, 100), rep(0, 100)),
    y2 = c(rep(-5, 100), rep(-3, 100), rep(5, 100), rep(0, 100))
  )

  res <- joint_individual_certainty(
    df_draws,
    id = id,
    values = c(y1, y2),
    mid = list(c(2, 4), c(2, 4)),
    .by = condition,
    better_is = "lower"
  )

  p_joint <- plot(res, mid_select = 1)
  expect_s3_class(p_joint, "ggplot")
  expect_s3_class(p_joint$facet, "FacetWrap")
  expect_identical(p_joint$labels$colour, "Signed Certainty")

  p_no_ci <- plot(
    res,
    mid_select = c(2, 2),
    show_ci = FALSE,
    show_rops = FALSE,
    color_by = "makowski_tier",
    facet = FALSE
  )
  expect_s3_class(p_no_ci, "ggplot")
  expect_s3_class(p_no_ci$facet, "FacetNull")
  expect_identical(p_no_ci$labels$colour, "Evidence Tier")

  p_sens <- plot(res, type = "sensitivity", thresholds = c(0.5, 0.9))
  expect_s3_class(p_sens, "ggplot")
  expect_s3_class(p_sens$layers[[1]]$geom, "GeomTile")

  p_sens_imp <- plot(
    res,
    type = "sensitivity",
    outcome = "improvement",
    thresholds = c(0.5, 0.9)
  )
  expect_s3_class(p_sens_imp, "ggplot")
  expect_identical(unique(p_sens_imp$data$direction), "Improvement")
})

test_that("joint_individual_certainty() validates arguments with informative cli errors", {
  df_valid <- tibble::tibble(
    id = c(1, 1),
    .draw = 1:2,
    y1 = c(-1, 2),
    y2 = c(3, -1)
  )

  expect_error(
    joint_individual_certainty("not_a_df", values = c(y1, y2), mid = 2),
    class = "peoplenotp_error_bad_data"
  )

  expect_error(
    joint_individual_certainty(df_valid, values = y1, mid = 2),
    class = "peoplenotp_error_bad_values"
  )

  expect_error(
    joint_individual_certainty(df_valid, values = c(y1, y2), mid = c(1, 2, 3)),
    class = "peoplenotp_error_bad_mid"
  )

  expect_error(
    joint_individual_certainty(
      df_valid,
      values = c(y1, y2),
      mid = 2,
      better_is = "invalid"
    ),
    class = "peoplenotp_error_bad_better_is"
  )

  expect_error(
    joint_individual_certainty(
      df_valid,
      values = c(y1, y2),
      mid = 2,
      include_marginals = NA
    ),
    class = "peoplenotp_error_bad_include_marginals"
  )
})
