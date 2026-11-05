test_that("get_joint_clinical_certainty() works for a single participant with reduced default columns", {
  set.seed(123)
  df_single <- tibble::tibble(
    wellbeing = stats::rnorm(1000, mean = 5, sd = 0.5),
    depression = stats::rnorm(1000, mean = -4, sd = 0.5)
  )

  res <- get_joint_clinical_certainty(
    df_single,
    col_1 = wellbeing,
    col_2 = depression,
    better_is = c("higher", "lower"),
    ranges = list(c(-2, 2), c(-1, 1))
  )

  expect_s3_class(res, "tbl_df")
  expect_identical(nrow(res), 1L)
  expect_named(
    res,
    c(
      "median_1",
      "lower_1",
      "upper_1",
      "median_2",
      "lower_2",
      "upper_2",
      "joint_pd",
      "joint_improvement_pct",
      "joint_deterioration_pct"
    )
  )
  expect_equal(res$joint_pd, 1)
  expect_equal(res$joint_improvement_pct, 1)
  expect_equal(res$joint_deterioration_pct, 0)
  expect_false("joint_clinical_status" %in% names(res))
  expect_false("pd_1" %in% names(res))
  expect_false("pd_2" %in% names(res))
})

test_that("get_joint_clinical_certainty() appends joint_clinical_status when include_status = TRUE", {
  set.seed(123)
  df_single <- tibble::tibble(
    wellbeing = stats::rnorm(1000, mean = 5, sd = 0.5),
    depression = stats::rnorm(1000, mean = -4, sd = 0.5)
  )

  res <- get_joint_clinical_certainty(
    df_single,
    col_1 = wellbeing,
    col_2 = depression,
    better_is = c("higher", "lower"),
    ranges = list(c(-2, 2), c(-1, 1)),
    include_status = TRUE
  )

  expect_true("joint_clinical_status" %in% names(res))
  expect_identical(res$joint_clinical_status, "Significant Improvement")
})

test_that("get_joint_clinical_certainty() appends marginal metrics when include_marginals = TRUE", {
  set.seed(123)
  df_single <- tibble::tibble(
    wellbeing = stats::rnorm(1000, mean = 5, sd = 0.5),
    depression = stats::rnorm(1000, mean = 0, sd = 0.2)
  )

  res <- get_joint_clinical_certainty(
    df_single,
    col_1 = wellbeing,
    col_2 = depression,
    better_is = c("higher", "lower"),
    ranges = list(c(-2, 2), c(-1, 1)),
    include_marginals = TRUE,
    include_status = TRUE
  )

  expect_named(
    res,
    c(
      "median_1",
      "lower_1",
      "upper_1",
      "pd_1",
      "improvement_pct_1",
      "deterioration_pct_1",
      "clinical_status_1",
      "median_2",
      "lower_2",
      "upper_2",
      "pd_2",
      "improvement_pct_2",
      "deterioration_pct_2",
      "clinical_status_2",
      "joint_pd",
      "joint_improvement_pct",
      "joint_deterioration_pct",
      "joint_clinical_status"
    )
  )
  expect_identical(res$clinical_status_1, "Significant Improvement")
  expect_identical(res$clinical_status_2, "Undecided")
  expect_identical(res$joint_clinical_status, "Undecided")
})

test_that("get_joint_clinical_certainty() supports .by for multi-participant summarization", {
  set.seed(42)
  df_multi <- tibble::tibble(
    id = rep(c("P1", "P2", "P3"), each = 1000),
    wellbeing = c(
      stats::rnorm(1000, mean = 5, sd = 0.5),
      stats::rnorm(1000, mean = 0, sd = 0.2),
      stats::rnorm(1000, mean = -5, sd = 0.5)
    ),
    depression = c(
      stats::rnorm(1000, mean = -4, sd = 0.5),
      stats::rnorm(1000, mean = 0, sd = 0.2),
      stats::rnorm(1000, mean = 4, sd = 0.5)
    )
  )

  res <- df_multi |>
    get_joint_clinical_certainty(
      col_1 = wellbeing,
      col_2 = depression,
      .by = id,
      better_is = c("higher", "lower"),
      ranges = list(c(-1, 1), c(-1, 1)),
      include_status = TRUE
    )

  expect_s3_class(res, "tbl_df")
  expect_false(dplyr::is_grouped_df(res))
  expect_identical(nrow(res), 3L)
  expect_identical(res$id, c("P1", "P2", "P3"))
  expect_identical(
    res$joint_clinical_status,
    c("Significant Improvement", "Undecided", "Significant Deterioration")
  )
})

test_that("get_joint_clinical_certainty() respects dplyr::group_by() and strips .row", {
  set.seed(99)
  df_epred <- tibble::tibble(
    id = rep(c("A", "B"), each = 500),
    .row = rep(1:2, each = 500),
    outcome_1 = c(
      stats::rnorm(500, mean = 4, sd = 0.5),
      stats::rnorm(500, mean = -4, sd = 0.5)
    ),
    outcome_2 = c(
      stats::rnorm(500, mean = 4, sd = 0.5),
      stats::rnorm(500, mean = -4, sd = 0.5)
    )
  ) |>
    dplyr::group_by(id, .row)

  res_grouped <- df_epred |>
    get_joint_clinical_certainty(
      col_1 = outcome_1,
      col_2 = outcome_2,
      better_is = c("higher", "higher"),
      ranges = list(c(-1, 1), c(-1, 1)),
      include_status = TRUE
    )

  expect_s3_class(res_grouped, "tbl_df")
  expect_false(dplyr::is_grouped_df(res_grouped))
  expect_identical(nrow(res_grouped), 2L)
  expect_identical(res_grouped$id, c("A", "B"))
  expect_false(".row" %in% names(res_grouped))
  expect_identical(
    res_grouped$joint_clinical_status,
    c("Significant Improvement", "Significant Deterioration")
  )

  # Providing .by overrides existing group_by()
  res_by_override <- df_epred |>
    get_joint_clinical_certainty(
      col_1 = outcome_1,
      col_2 = outcome_2,
      .by = id,
      better_is = c("higher", "higher"),
      ranges = list(c(-1, 1), c(-1, 1)),
      include_status = TRUE
    )

  expect_identical(res_grouped, res_by_override)
})

test_that("get_joint_clinical_certainty() validates arguments with informative cli errors", {
  df_valid <- tibble::tibble(y1 = c(-1, 0, 1), y2 = c(-1, 0, 1))

  expect_error(
    get_joint_clinical_certainty("not_a_df", col_1 = y1, col_2 = y2),
    class = "peoplenotp_error_bad_data"
  )

  expect_error(
    get_joint_clinical_certainty(df_valid, col_1 = missing_col, col_2 = y2),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    get_joint_clinical_certainty(df_valid, col_1 = y1, col_2 = missing_col),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    get_joint_clinical_certainty(
      tibble::tibble(y1 = c("a", "b"), y2 = c(1, 2)),
      col_1 = y1,
      col_2 = y2
    ),
    class = "peoplenotp_error_bad_column_type"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      better_is = "higher"
    ),
    class = "peoplenotp_error_bad_better_is"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      ranges = list(c(2, -2), c(-1, 1))
    ),
    class = "peoplenotp_error_bad_range"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      ci_width = 1.5
    ),
    class = "peoplenotp_error_bad_ci_width"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      certainty_boundaries = c(0.99, 0.95)
    ),
    class = "peoplenotp_error_bad_certainty_boundaries"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      include_marginals = NA
    ),
    class = "peoplenotp_error_bad_include_marginals"
  )

  expect_error(
    get_joint_clinical_certainty(
      df_valid,
      col_1 = y1,
      col_2 = y2,
      include_status = NA
    ),
    class = "peoplenotp_error_bad_include_status"
  )
})
