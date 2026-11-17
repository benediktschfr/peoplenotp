test_that("group_certainty() summarizes fixed effects (rope_pct, pd, CrI) and extracts random effects", {
  set.seed(123)
  mock_data <- tibble::tibble(
    day = 1:10,
    score = stats::rnorm(10)
  )
  mock_model <- stats::lm(score ~ day, data = mock_data)

  all_model_vars <- c(
    "b_Intercept",
    "b_time",
    "sd_id__Intercept",
    "sd_id__time",
    "cor_id__Intercept__time",
    "sigma",
    "r_id[1,Intercept]"
  )

  draws_lookup <- list(
    b_Intercept = stats::rnorm(1000, mean = 20, sd = 1),
    b_time = stats::rnorm(1000, mean = -0.5, sd = 0.2),
    sd_id__Intercept = stats::rnorm(1000, mean = 3, sd = 0.2),
    sd_id__time = stats::rnorm(1000, mean = 1.2, sd = 0.1),
    cor_id__Intercept__time = stats::rnorm(1000, mean = -0.3, sd = 0.05)
  )

  local_mocked_bindings(
    get_variables = function(model, ...) {
      all_model_vars
    },
    gather_draws = function(model, ..., regex = FALSE) {
      dots <- rlang::enquos(...)
      pat <- rlang::as_name(dots[[1]])
      matched_vars <- grep(pat, names(draws_lookup), value = TRUE)
      purrr::map_dfr(matched_vars, function(v) {
        tibble::tibble(
          .variable = v,
          .draw = seq_along(draws_lookup[[v]]),
          .value = draws_lookup[[v]]
        )
      })
    },
    .package = "tidybayes"
  )

  res <- group_certainty(
    mock_model,
    rope = c(-2.5, 2.5),
    ci = 0.95
  )

  expect_s3_class(res, c("pnp_group_certainty", "data.frame"))
  expect_named(
    res,
    c(
      "parameter",
      "effect",
      "median",
      "CI_low",
      "CI_high",
      "pd",
      "rope_pct"
    )
  )
  expect_false(any(c("improvement_pct", "deterioration_pct") %in% names(res)))

  fixed_rows <- res[res$effect == "fixed", ]
  random_rows <- res[res$effect == "random", ]

  expect_identical(fixed_rows$parameter, c("b_Intercept", "b_time"))
  expect_identical(
    random_rows$parameter,
    c("sd_id__Intercept", "sd_id__time", "cor_id__Intercept__time")
  )

  # b_Intercept ~ N(20, 1) lies outside [-2.5, 2.5]; b_time ~ N(-0.5, 0.2) lies inside [-2.5, 2.5]
  expect_equal(fixed_rows$rope_pct, c(0, 1))

  # Random effects have median & CrI, and NA for pd/rope_pct
  expect_true(all(!is.na(random_rows$median)))
  expect_true(all(!is.na(random_rows$CI_low)))
  expect_true(all(!is.na(random_rows$CI_high)))
  expect_true(all(is.na(random_rows$rope_pct)))
})

test_that("group_certainty() supports parameter filtering and print()", {
  set.seed(42)
  draws_df <- tibble::tibble(
    .draw = 1:1000,
    b_Intercept = stats::rnorm(1000, mean = 10, sd = 1),
    b_treatment = stats::rnorm(1000, mean = 4, sd = 0.4),
    b_age = stats::rnorm(1000, mean = 0, sd = 0.2),
    sd_id__Intercept = stats::rnorm(1000, mean = 2, sd = 0.1)
  )

  res <- group_certainty(
    draws_df,
    rope = c(-2, 2),
    parameters = "treatment"
  )

  expect_s3_class(res, c("pnp_group_certainty", "data.frame"))
  fixed_rows <- res[res$effect == "fixed", ]
  expect_identical(fixed_rows$parameter, "b_treatment")
  expect_equal(fixed_rows$rope_pct, 0)

  printed <- capture.output(ret <- print(res))
  expect_identical(ret, res)
  expect_match(
    paste(printed, collapse = "\n"),
    "Group-Level Clinical Certainty \\(Fixed Effects\\)"
  )
  expect_match(
    paste(printed, collapse = "\n"),
    "Individual Heterogeneity \\(Random Effects\\)"
  )
})

test_that("group_certainty() validates arguments with informative cli errors", {
  mock_data <- tibble::tibble(x = 1:5, y = 1:5)
  valid_model <- stats::lm(y ~ x, data = mock_data)

  local_mocked_bindings(
    get_variables = function(model, ...) c("b_Intercept", "b_x"),
    .package = "tidybayes"
  )

  expect_error(
    group_certainty("not_a_model", rope = c(-1, 1)),
    class = "peoplenotp_error_bad_model"
  )

  expect_error(
    group_certainty(valid_model, rope = c(2, -2)),
    class = "peoplenotp_error_bad_rope"
  )

  expect_error(
    group_certainty(valid_model, rope = 2),
    class = "peoplenotp_error_bad_rope"
  )

  expect_error(
    group_certainty(valid_model, rope = c(-1, 1), parameters = ""),
    class = "peoplenotp_error_bad_parameters"
  )

  expect_error(
    group_certainty(valid_model, rope = c(-1, 1), ci = 1.2),
    class = "peoplenotp_error_bad_ci"
  )

  expect_error(
    group_certainty(
      valid_model,
      rope = c(-1, 1),
      parameters = "^nonexistent$"
    ),
    class = "peoplenotp_error_no_matching_variables"
  )
})
