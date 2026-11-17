test_that("plot_clinical_sensitivity() visualizes sensitivity grid output from get_clinical_sensitivity_grid()", {
  set.seed(123)
  draws_df <- tibble::tibble(
    id = rep(1:10, each = 50),
    .epred = stats::rnorm(500, mean = 12, sd = 3)
  )

  grid_data <- get_clinical_sensitivity_grid(
    draws_df,
    rope_range = c(8, 12),
    rope_step = 2,
    threshold_range = c(0.8, 0.9),
    threshold_step = 0.1
  )

  p <- plot_clinical_sensitivity(grid_data, color_palette = "viridis")
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2L)

  p_custom <- plot_clinical_sensitivity(grid_data, color_palette = "oasis")
  expect_s3_class(p_custom, "ggplot")
})

test_that("plot_clinical_sensitivity() validates arguments with informative cli errors", {
  valid_grid <- tibble::tibble(
    rope = c(1, 2),
    threshold = c(0.8, 0.9),
    prop_certain = c(0.7, 0.5)
  )

  expect_error(
    plot_clinical_sensitivity("not_a_df"),
    class = "peoplenotp_error_bad_data"
  )

  expect_error(
    plot_clinical_sensitivity(tibble::tibble(rope = 1, threshold = 0.8)),
    class = "peoplenotp_error_missing_column"
  )

  expect_error(
    plot_clinical_sensitivity(valid_grid, color_palette = "nonexistent_pal"),
    class = "peoplenotp_error_bad_color_palette"
  )
})
