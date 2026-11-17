test_that("individual_certainty() computes summaries, Makowski tiers, and inverted evidence matrix", {
  # Construct deterministic draws for 5 individuals across 200 draws each
  # so each individual lands in a specific Makowski tier at MID = 2.5
  draws_p1 <- rep(0, 200) # 100% in ROPE -> Negligible (>99% ROPE)
  draws_p2 <- c(rep(0, 197), rep(-4, 3)) # 98.5% in ROPE -> Prob. Negligible (>97.5%)
  draws_p3 <- c(rep(0, 100), rep(-4, 100)) # 50% in ROPE -> Undecided
  draws_p4 <- c(rep(0, 3), rep(-4, 197)) # 1.5% in ROPE -> Prob. Significant (<2.5%)
  draws_p5 <- rep(-5, 200) # 0% in ROPE -> Significant (<1% ROPE)

  df_draws <- tibble::tibble(
    participant = rep(paste0("P", 1:5), each = 200),
    .draw = rep(1:200, times = 5),
    net_change = c(draws_p1, draws_p2, draws_p3, draws_p4, draws_p5)
  )

  res <- individual_certainty(
    data = df_draws,
    id = participant,
    value = net_change,
    mid = c(2.0, 2.5, 3.0),
    better_is = "lower",
    draw_col = .draw,
    ci = 0.95
  )

  expect_s3_class(res, "pnp_individual_certainty")
  expect_named(res, c("individual_summary", "matrix_data", "meta"))

  # Verify $individual_summary structure and calculations
  ind_sum <- res$individual_summary
  expect_s3_class(ind_sum, "tbl_df")
  expect_identical(nrow(ind_sum), 15L) # 5 individuals * 3 MIDs
  expect_named(
    ind_sum,
    c(
      "participant",
      "mid",
      "median",
      "ci_low",
      "ci_high",
      "pd",
      "rope_pct",
      "improvement_pct",
      "deterioration_pct",
      "directed_certainty",
      "makowski_tier"
    )
  )

  # Full posterior partitioning sums to 1
  expect_equal(
    ind_sum$rope_pct + ind_sum$improvement_pct + ind_sum$deterioration_pct,
    rep(1, nrow(ind_sum))
  )

  # Directed certainty equals max(improvement_pct, deterioration_pct)
  expect_equal(
    ind_sum$directed_certainty,
    pmax(ind_sum$improvement_pct, ind_sum$deterioration_pct)
  )

  # Check ordered factor levels for makowski_tier at MID = 2.5
  sub_25 <- ind_sum[ind_sum$mid == 2.5, ]
  expect_true(is.ordered(sub_25$makowski_tier))
  expect_identical(
    levels(sub_25$makowski_tier),
    c(
      "Negligible (>99% ROPE)",
      "Prob. Negligible (>97.5%)",
      "Undecided",
      "Prob. Significant (<2.5%)",
      "Significant (<1% ROPE)"
    )
  )
  expect_identical(
    as.character(sub_25$makowski_tier),
    c(
      "Negligible (>99% ROPE)",
      "Prob. Negligible (>97.5%)",
      "Undecided",
      "Prob. Significant (<2.5%)",
      "Significant (<1% ROPE)"
    )
  )

  # Verify $matrix_data structure (5 rows x MID columns)
  mat <- res$matrix_data
  expect_s3_class(mat, "tbl_df")
  expect_identical(nrow(mat), 5L)
  expect_named(mat, c("makowski_tier", "MID = 2.0", "MID = 2.5", "MID = 3.0"))
  expect_identical(mat$`MID = 2.5`, rep("1 (20.0%)", 5))
})

test_that("print.pnp_individual_certainty() outputs the inverted Evidence Matrix via insight::export_table()", {
  df_draws <- tibble::tibble(
    id = rep(c("P1", "P2"), each = 100),
    .draw = rep(1:100, times = 2),
    val = c(rep(0, 100), rep(5, 100))
  )

  res <- individual_certainty(
    data = df_draws,
    id = id,
    value = val,
    mid = c(2.0, 3.0),
    better_is = "higher"
  )

  out_both <- capture.output(ret <- print(res, format = "both"))
  expect_identical(ret, res)
  text_both <- paste(out_both, collapse = "\n")
  expect_match(text_both, "# Clinical Significance Evidence Matrix")
  expect_match(text_both, "Evidence Tier")
  expect_match(text_both, "MID = 2.0")
  expect_match(text_both, "MID = 3.0")
  expect_match(text_both, "1 \\(50\\.0%\\)")

  out_count <- capture.output(print(res, format = "count"))
  expect_false(any(grepl("50\\.0%", out_count)))

  out_prop <- capture.output(print(res, format = "prop"))
  expect_match(paste(out_prop, collapse = "\n"), "50\\.0%")
})

test_that("plot.pnp_individual_certainty() builds a volcano plot with ROPE and directional annotations", {
  set.seed(99)
  df_draws <- tibble::tibble(
    id = rep(paste0("P", 1:6), each = 100),
    .draw = rep(1:100, times = 6),
    change = stats::rnorm(
      600,
      mean = rep(c(-5, -3, 0, 1, 3, 5), each = 100),
      sd = 0.8
    )
  )

  res_lower <- individual_certainty(
    df_draws,
    id = id,
    value = change,
    mid = c(2.0, 2.5),
    better_is = "lower"
  )

  p_lower <- plot(
    res_lower,
    mid_select = 2.5,
    point_alpha = 0.8,
    interval_alpha = 0.4
  )
  expect_s3_class(p_lower, "ggplot")

  layer_labels <- unlist(lapply(p_lower$layers, function(l) {
    l$aes_params$label
  }))
  expect_in("\u2190 Favors Improvement", layer_labels)
  expect_in("Favors Deterioration \u2192", layer_labels)

  res_higher <- individual_certainty(
    df_draws,
    id = id,
    value = change,
    mid = 2.5,
    better_is = "higher"
  )
  p_higher <- plot(res_higher)
  layer_labels_higher <- unlist(lapply(p_higher$layers, function(l) {
    l$aes_params$label
  }))
  expect_in("\u2190 Favors Deterioration", layer_labels_higher)
  expect_in("Favors Improvement \u2192", layer_labels_higher)

  expect_error(
    plot(res_lower, mid_select = 99),
    class = "peoplenotp_error_bad_mid_select"
  )

  expect_error(
    plot(res_lower, point_alpha = 2),
    class = "peoplenotp_error_bad_alpha"
  )
})

test_that("individual_certainty() validates arguments with informative cli errors", {
  df_valid <- tibble::tibble(id = c(1, 1), .draw = 1:2, val = c(-1, 2))

  expect_error(
    individual_certainty("not_a_df", id = id, value = val, mid = 2),
    class = "peoplenotp_error_bad_data"
  )

  expect_error(
    individual_certainty(df_valid, id = missing_id, value = val, mid = 2),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    individual_certainty(df_valid, id = id, value = missing_val, mid = 2),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    individual_certainty(
      df_valid,
      id = id,
      value = val,
      mid = 2,
      draw_col = missing_draw
    ),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    individual_certainty(
      tibble::tibble(id = 1, .draw = 1, val = "a"),
      id = id,
      value = val,
      mid = 2
    ),
    class = "peoplenotp_error_bad_column_type"
  )

  expect_error(
    individual_certainty(df_valid, id = id, value = val, mid = -1),
    class = "peoplenotp_error_bad_mid"
  )

  expect_error(
    individual_certainty(
      df_valid,
      id = id,
      value = val,
      mid = 2,
      better_is = "invalid"
    ),
    class = "peoplenotp_error_bad_better_is"
  )

  expect_error(
    individual_certainty(df_valid, id = id, value = val, mid = 2, ci = 0),
    class = "peoplenotp_error_bad_ci"
  )
})
