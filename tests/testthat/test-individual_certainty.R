test_that("individual_certainty() computes summaries, directed Makowski tiers, and inverted evidence matrix", {
  # Construct deterministic draws for 6 individuals across 200 draws each
  # so each individual lands in one of the 6 directed Makowski tiers at MID = 2.5 (better_is = "lower")
  draws_p1 <- rep(5, 200) # 100% det -> Definite Deterioration (>99%)
  draws_p2 <- c(rep(5, 197), rep(0, 3)) # 98.5% det -> Probable Deterioration (>97.5%)
  draws_p3 <- c(rep(0, 100), rep(-4, 100)) # 50% ROPE, 50% imp -> Undecided
  draws_p4 <- rep(0, 200) # 100% in ROPE -> Practically Equivalent (>97.5% ROPE)
  draws_p5 <- c(rep(0, 3), rep(-4, 197)) # 98.5% imp -> Probable Improvement (>97.5%)
  draws_p6 <- rep(-5, 200) # 100% imp -> Definite Improvement (>99%)

  df_draws <- tibble::tibble(
    participant = rep(paste0("P", 1:6), each = 200),
    .draw = rep(1:200, times = 6),
    net_change = c(
      draws_p1,
      draws_p2,
      draws_p3,
      draws_p4,
      draws_p5,
      draws_p6
    )
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
  expect_identical(nrow(ind_sum), 18L) # 6 individuals * 3 MIDs
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
      "signed_certainty",
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

  # Signed certainty equals improvement_pct - deterioration_pct in [-1, 1]
  expect_equal(
    ind_sum$signed_certainty,
    ind_sum$improvement_pct - ind_sum$deterioration_pct
  )

  # Check ordered factor levels for makowski_tier at MID = 2.5
  expected_levels <- c(
    "Definite Deterioration (>99%)",
    "Probable Deterioration (>97.5%)",
    "Undecided",
    "Practically Equivalent (>97.5% ROPE)",
    "Probable Improvement (>97.5%)",
    "Definite Improvement (>99%)"
  )
  sub_25 <- ind_sum[ind_sum$mid == 2.5, ]
  expect_true(is.ordered(sub_25$makowski_tier))
  expect_identical(levels(sub_25$makowski_tier), expected_levels)
  expect_identical(as.character(sub_25$makowski_tier), expected_levels)

  # Verify $matrix_data structure (6 rows x MID columns)
  mat <- res$matrix_data
  expect_s3_class(mat, "tbl_df")
  expect_identical(nrow(mat), 6L)
  expect_named(mat, c("makowski_tier", "MID = 2.0", "MID = 2.5", "MID = 3.0"))
  expect_identical(mat$`MID = 2.5`, rep("1 (16.7%)", 6))
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

  # Auto min/middle/max selection when > 3 MIDs are supplied and mid_select = NULL
  res_many <- individual_certainty(
    data = df_draws,
    id = id,
    value = val,
    mid = c(1.0, 2.0, 3.0, 4.0, 5.0),
    better_is = "higher"
  )
  out_many <- paste(capture.output(print(res_many)), collapse = "\n")
  expect_match(out_many, "MID = 1\\.0")
  expect_match(out_many, "MID = 3\\.0")
  expect_match(out_many, "MID = 5\\.0")
  expect_false(grepl("MID = 2\\.0", out_many))
  expect_false(grepl("MID = 4\\.0", out_many))

  # Explicit mid_select vector in print()
  out_explicit <- paste(
    capture.output(print(res_many, mid_select = c(2.0, 4.0))),
    collapse = "\n"
  )
  expect_match(out_explicit, "MID = 2\\.0")
  expect_match(out_explicit, "MID = 4\\.0")
  expect_false(grepl("MID = 1\\.0", out_explicit))

  expect_error(
    print(res_many, mid_select = 99),
    class = "peoplenotp_error_bad_mid_select"
  )
})

test_that("individual_certainty() supports optional .by grouping and dplyr::group_by()", {
  df_grouped <- tibble::tibble(
    id = rep(paste0("P", 1:4), each = 100),
    condition = rep(
      c("Waitlist", "Waitlist", "Intervention", "Intervention"),
      each = 100
    ),
    .draw = rep(1:100, times = 4),
    .epred = c(
      rep(0, 100), # P1 Waitlist: Negligible
      rep(0, 100), # P2 Waitlist: Negligible
      rep(-5, 100), # P3 Intervention: Significant
      rep(-5, 100) # P4 Intervention: Significant
    )
  )

  res_by <- individual_certainty(
    df_grouped,
    mid = c(2.0, 3.0),
    .by = condition,
    better_is = "lower"
  )

  expect_identical(res_by$meta$id_col, "id")
  expect_identical(res_by$meta$group_cols, "condition")
  expect_true("condition" %in% names(res_by$individual_summary))
  expect_identical(nrow(res_by$matrix_data), 12L) # 6 tiers * 2 conditions

  printed_by <- paste(capture.output(print(res_by)), collapse = "\n")
  expect_match(printed_by, "condition = Waitlist")
  expect_match(printed_by, "condition = Intervention")
  expect_match(printed_by, "2 \\(100\\.0%\\)")

  p_faceted <- plot(res_by, mid_select = 2.0)
  expect_s3_class(p_faceted$facet, "FacetWrap")

  p_unfaceted <- plot(res_by, mid_select = 2.0, facet = FALSE)
  expect_s3_class(p_unfaceted$facet, "FacetNull")

  # Also works via dplyr::group_by() and strips .row
  res_grp <- df_grouped |>
    dplyr::mutate(.row = 1L) |>
    dplyr::group_by(id, condition, .row) |>
    individual_certainty(mid = c(2.0, 3.0), better_is = "lower")

  expect_identical(res_grp$individual_summary, res_by$individual_summary)
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
  expect_identical(p_lower$labels$colour, "Signed Certainty")

  p_tier <- plot(res_lower, mid_select = 2.5, color_by = makowski_tier)
  expect_s3_class(p_tier, "ggplot")
  expect_identical(p_tier$labels$colour, "Evidence Tier")

  p_imp <- plot(res_lower, mid_select = 2.5, color_by = "improvement_pct")
  expect_s3_class(p_imp, "ggplot")
  expect_identical(p_imp$labels$colour, "% Improvement")

  # show_ci = FALSE omits the geom_errorbar layer and updates x-axis label
  p_no_ci <- plot(res_lower, mid_select = 2.5, show_ci = FALSE)
  expect_s3_class(p_no_ci, "ggplot")
  expect_identical(length(p_no_ci$layers), length(p_lower$layers) - 1L)
  expect_identical(p_no_ci$labels$x, "Posterior Median")

  expect_error(
    plot(res_lower, show_ci = NA),
    class = "peoplenotp_error_bad_flag"
  )

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
    plot(res_lower, color_by = "invalid_col"),
    class = "peoplenotp_error_bad_color_by"
  )

  expect_error(
    plot(res_lower, point_alpha = 2),
    class = "peoplenotp_error_bad_alpha"
  )

  expect_error(
    plot(res_lower, type = "invalid_type"),
    class = "peoplenotp_error_bad_type"
  )
})

test_that("plot.pnp_individual_certainty(type = 'sensitivity') builds mirrored and single-sided sensitivity plots with grouping", {
  # P1: 100% improvement at MID <= 4 (change = -5, better_is = "lower")
  # P2: 80% improvement at MID = 2, 0% at MID = 4
  # P3: 100% deterioration at MID <= 4 (change = +5, better_is = "lower")
  # P4: 0% improvement/deterioration (change = 0)
  df_draws <- tibble::tibble(
    id = rep(paste0("P", 1:4), each = 100),
    condition = rep(c("GroupA", "GroupA", "GroupB", "GroupB"), each = 100),
    .draw = rep(1:100, times = 4),
    change = c(
      rep(-5, 100),
      c(rep(-3, 80), rep(0, 20)),
      rep(5, 100),
      rep(0, 100)
    )
  )

  res <- individual_certainty(
    df_draws,
    id = id,
    value = change,
    mid = c(2.0, 4.0),
    .by = condition,
    better_is = "lower"
  )

  p_both <- plot(
    res,
    type = "sensitivity",
    thresholds = c(0.5, 0.9)
  )
  expect_s3_class(p_both, "ggplot")
  expect_s3_class(p_both$facet, "FacetWrap")
  expect_s3_class(p_both$layers[[1]]$geom, "GeomTile")
  expect_identical(p_both$labels$x, "Minimal Important Difference (MID)")
  expect_identical(p_both$labels$y, "Certainty Threshold")
  expect_identical(p_both$labels$fill, "Proportion Certain")

  # Check directional annotations in butterfly plot
  sens_labels <- unlist(lapply(p_both$layers, function(l) {
    l$aes_params$label
  }))
  expect_in("\u2191 Improvement", sens_labels)
  expect_in("\u2193 Deterioration", sens_labels)

  # Verify mirrored values in plot data: GroupA has positive prop_signed, GroupB has negative prop_signed
  df_plot_both <- p_both$data
  expect_true(all(
    c("Improvement", "Deterioration") %in% df_plot_both$direction
  ))
  grp_a_imp <- df_plot_both[
    df_plot_both$condition == "GroupA" &
      df_plot_both$direction == "Improvement" &
      df_plot_both$mid == 2.0 &
      df_plot_both$threshold == 0.5,
  ]
  expect_equal(grp_a_imp$prop_signed, 1.0) # Both P1 (100%) and P2 (80%) >= 0.5
  expect_true(grp_a_imp$y_tile > 0)

  grp_b_det <- df_plot_both[
    df_plot_both$condition == "GroupB" &
      df_plot_both$direction == "Deterioration" &
      df_plot_both$mid == 2.0 &
      df_plot_both$threshold == 0.9,
  ]
  expect_equal(grp_b_det$prop_signed, -0.5) # P3 (100%) >= 0.9 out of 2 in GroupB -> -0.5
  expect_true(grp_b_det$y_tile < 0)

  # Single-sided outcome = "improvement" and outcome = "deterioration"
  p_imp <- plot(
    res,
    type = "sensitivity",
    outcome = "improvement",
    color_palette = "viridis",
    facet = FALSE
  )
  expect_s3_class(p_imp, "ggplot")
  expect_s3_class(p_imp$facet, "FacetNull")
  expect_true(all(p_imp$data$prop_certain >= 0))
  expect_identical(unique(p_imp$data$direction), "Improvement")

  p_det <- plot(
    res,
    type = "sensitivity",
    outcome = "deterioration",
    color_palette = "oasis",
    mid_select = 2.0
  )
  expect_s3_class(p_det, "ggplot")
  expect_identical(unique(p_det$data$mid), 2.0)
  expect_identical(unique(p_det$data$direction), "Deterioration")

  # Validation errors for sensitivity plot arguments
  expect_error(
    plot(res, type = "sensitivity", outcome = "invalid"),
    class = "peoplenotp_error_bad_outcome"
  )
  expect_error(
    plot(res, type = "sensitivity", thresholds = c(0, 0.5)),
    class = "peoplenotp_error_bad_thresholds"
  )
  expect_error(
    plot(res, type = "sensitivity", thresholds = c(0.5, 1.2)),
    class = "peoplenotp_error_bad_thresholds"
  )
  expect_error(
    plot(res, type = "sensitivity", color_palette = "nonexistent"),
    class = "peoplenotp_error_bad_color_palette"
  )
  expect_error(
    plot(res, type = "sensitivity", mid_select = 99),
    class = "peoplenotp_error_bad_mid_select"
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

test_that("individual_certainty() supports single-case (N = 1) datasets without an id column", {
  df_single <- tibble::tibble(
    .draw = 1:100,
    .epred = rep(-5, 100)
  )

  res_single <- individual_certainty(
    data = df_single,
    mid = c(2.0, 3.0),
    better_is = "lower"
  )

  expect_s3_class(res_single, "pnp_individual_certainty")
  expect_null(res_single$meta$id_col)
  expect_identical(res_single$meta$n_individuals, 1L)
  expect_identical(nrow(res_single$individual_summary), 2L)
  expect_equal(res_single$individual_summary$improvement_pct, c(1, 1))
})
