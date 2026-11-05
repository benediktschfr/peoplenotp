test_that("get_clinical_certainty() works for a single participant with reduced default columns", {
  set.seed(123)
  df_single <- tibble::tibble(
    iteration = 1:1000,
    score_difference = stats::rnorm(1000, mean = 5, sd = 0.5)
  )

  res <- get_clinical_certainty(
    df_single,
    col = score_difference,
    better_is = "higher",
    range = c(-2, 2)
  )

  expect_s3_class(res, "tbl_df")
  expect_identical(nrow(res), 1L)
  expect_named(
    res,
    c(
      "median",
      "lower",
      "upper",
      "pd",
      "improvement_pct",
      "deterioration_pct"
    )
  )
  expect_equal(res$pd, 1)
  expect_equal(res$improvement_pct, 1)
  expect_equal(res$deterioration_pct, 0)
  expect_false("clinical_status" %in% names(res))
  expect_false(any(c("rope_low", "rope_high", "rope_pct") %in% names(res)))
})

test_that("get_clinical_certainty() includes clinical_status when include_status = TRUE", {
  set.seed(123)
  df_single <- tibble::tibble(
    score_difference = stats::rnorm(1000, mean = 5, sd = 0.5)
  )

  res <- get_clinical_certainty(
    df_single,
    col = score_difference,
    range = c(-2, 2),
    include_status = TRUE
  )

  expect_true("clinical_status" %in% names(res))
  expect_identical(res$clinical_status, "Significant Improvement")
})

test_that("get_clinical_certainty() includes ROPE bounds when include_rope = TRUE", {
  set.seed(123)
  df_single <- tibble::tibble(
    score_difference = stats::rnorm(1000, mean = 5, sd = 0.5)
  )

  res <- get_clinical_certainty(
    df_single,
    col = score_difference,
    range = c(-2, 2),
    include_rope = TRUE
  )

  expect_true(all(c("rope_low", "rope_high", "rope_pct") %in% names(res)))
  expect_identical(res$rope_low, -2)
  expect_identical(res$rope_high, 2)
  expect_equal(res$rope_pct, 0)
})

test_that("get_clinical_certainty() defaults col to .epred and supports .by", {
  set.seed(42)
  df_multi <- tibble::tibble(
    id = rep(c("P1", "P2", "P3"), each = 1000),
    .epred = c(
      stats::rnorm(1000, mean = 5, sd = 0.5),
      stats::rnorm(1000, mean = 0, sd = 0.2),
      stats::rnorm(1000, mean = -5, sd = 0.5)
    )
  )

  res <- df_multi |>
    get_clinical_certainty(
      .by = id,
      better_is = "higher",
      range = c(-1, 1),
      include_status = TRUE
    )

  expect_s3_class(res, "tbl_df")
  expect_false(dplyr::is_grouped_df(res))
  expect_identical(nrow(res), 3L)
  expect_identical(res$id, c("P1", "P2", "P3"))
  expect_identical(
    res$clinical_status,
    c("Significant Improvement", "Undecided", "Significant Deterioration")
  )
})

test_that("get_clinical_certainty() respects dplyr::group_by() and strips .row from add_epred_draws()", {
  set.seed(99)
  df_epred <- tibble::tibble(
    id = rep(c("A", "B"), each = 500),
    .row = rep(1:2, each = 500),
    .epred = c(
      stats::rnorm(500, mean = 4, sd = 0.5),
      stats::rnorm(500, mean = -4, sd = 0.5)
    )
  ) |>
    dplyr::group_by(id, .row)

  res_grouped <- df_epred |>
    get_clinical_certainty(range = c(-1, 1), include_status = TRUE)

  expect_s3_class(res_grouped, "tbl_df")
  expect_false(dplyr::is_grouped_df(res_grouped))
  expect_identical(nrow(res_grouped), 2L)
  expect_identical(res_grouped$id, c("A", "B"))
  expect_false(".row" %in% names(res_grouped))
  expect_identical(
    res_grouped$clinical_status,
    c("Significant Improvement", "Significant Deterioration")
  )

  # Providing .by overrides existing group_by()
  res_by_override <- df_epred |>
    get_clinical_certainty(.by = id, range = c(-1, 1), include_status = TRUE)

  expect_identical(res_grouped, res_by_override)
})

test_that("get_clinical_certainty() handles better_is = 'lower'", {
  set.seed(101)
  df_pain <- tibble::tibble(
    .epred = stats::rnorm(1000, mean = -4, sd = 0.5)
  )

  res <- get_clinical_certainty(
    df_pain,
    better_is = "lower",
    range = c(-1, 1),
    include_status = TRUE
  )

  expect_identical(res$clinical_status, "Significant Improvement")
  expect_false(any(c("rope_low", "rope_high", "rope_pct") %in% names(res)))
  expect_equal(res$improvement_pct, 1)
  expect_equal(res$deterioration_pct, 0)
})

test_that("get_clinical_certainty() validates arguments with informative cli errors", {
  df_valid <- tibble::tibble(.epred = c(-1, 0, 1))

  expect_error(
    get_clinical_certainty("not_a_df"),
    class = "peoplenotp_error_bad_data"
  )

  expect_error(
    get_clinical_certainty(df_valid, col = missing_col),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    get_clinical_certainty(tibble::tibble(.epred = c("a", "b"))),
    class = "peoplenotp_error_bad_column_type"
  )

  expect_error(
    get_clinical_certainty(df_valid, better_is = "invalid"),
    class = "peoplenotp_error_bad_better_is"
  )

  expect_error(
    get_clinical_certainty(df_valid, range = c(2, -2)),
    class = "peoplenotp_error_bad_range"
  )

  expect_error(
    get_clinical_certainty(df_valid, ci_width = 1.5),
    class = "peoplenotp_error_bad_ci_width"
  )

  expect_error(
    get_clinical_certainty(df_valid, certainty_boundaries = c(0.99, 0.95)),
    class = "peoplenotp_error_bad_certainty_boundaries"
  )

  expect_error(
    get_clinical_certainty(df_valid, include_rope = NA),
    class = "peoplenotp_error_bad_include_rope"
  )

  expect_error(
    get_clinical_certainty(df_valid, include_status = NA),
    class = "peoplenotp_error_bad_include_status"
  )
})
