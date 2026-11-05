test_that("get_traceplot() creates sectioned patchwork plots for enabled parameter groups and omits random effects", {
  set.seed(123)
  mock_data <- tibble::tibble(
    day = 1:10,
    severity_score = stats::rnorm(10, mean = 30, sd = 5)
  )
  mock_model <- stats::lm(severity_score ~ day, data = mock_data)

  all_model_vars <- c(
    "b_Intercept",
    "b_day",
    "sd_id__Intercept",
    "sigma",
    "ar[1]",
    "r_id[1,Intercept]",
    "r_id[2,Intercept]",
    "lp__"
  )

  queried_patterns <- character(0)

  local_mocked_bindings(
    get_variables = function(model, ...) {
      all_model_vars
    },
    gather_draws = function(model, ..., regex = FALSE) {
      dots <- rlang::enquos(...)
      pat <- rlang::as_name(dots[[1]])
      queried_patterns <<- c(queried_patterns, pat)
      matched_vars <- grep(pat, all_model_vars, value = TRUE)
      tidyr::expand_grid(
        .variable = matched_vars,
        .chain = 1:2,
        .iteration = 1:20
      ) |>
        dplyr::mutate(
          .draw = dplyr::row_number(),
          .value = stats::rnorm(dplyr::n())
        )
    },
    .package = "tidybayes"
  )

  p <- get_traceplot(
    mock_model,
    alpha = 0.7,
    linewidth = 0.8,
    scales = "free",
    facet_ncol = 2,
    patchwork_ncol = 1
  )

  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "ggplot")
  expect_identical(
    queried_patterns,
    c("^b_", "^sd_", "^sigma", "^(ar(\\[|$)|cor_ar)")
  )
  # 4 sections: 3 in p$patches$plots and the last active plot in p itself
  expect_length(p$patches$plots, 3L)
  expect_equal(p$layers[[1]]$aes_params$alpha, 0.7)
  expect_equal(p$layers[[1]]$aes_params$linewidth, 0.8)
  expect_equal(p$facet$params$ncol, 2)
})

test_that("get_traceplot() silently skips enabled parameter groups that do not exist in the model", {
  mock_data <- tibble::tibble(
    day = 1:10,
    severity_score = stats::rnorm(10)
  )
  mock_model <- stats::lm(severity_score ~ day, data = mock_data)

  # Only fixed effects and sigma present (no sd_ or ar)
  simple_vars <- c("b_Intercept", "b_day", "sigma", "lp__")
  queried_patterns <- character(0)

  local_mocked_bindings(
    get_variables = function(model, ...) {
      simple_vars
    },
    gather_draws = function(model, ..., regex = FALSE) {
      dots <- rlang::enquos(...)
      pat <- rlang::as_name(dots[[1]])
      queried_patterns <<- c(queried_patterns, pat)
      matched_vars <- grep(pat, simple_vars, value = TRUE)
      tidyr::expand_grid(
        .variable = matched_vars,
        .chain = 1:2,
        .iteration = 1:10
      ) |>
        dplyr::mutate(
          .draw = dplyr::row_number(),
          .value = stats::rnorm(dplyr::n())
        )
    },
    .package = "tidybayes"
  )

  p <- get_traceplot(mock_model)

  expect_s3_class(p, "patchwork")
  expect_identical(queried_patterns, c("^b_", "^sigma"))
  expect_length(p$patches$plots, 1L)
})

test_that("get_traceplot() respects individual boolean flags and custom regex override", {
  mock_data <- tibble::tibble(
    day = 1:10,
    severity_score = stats::rnorm(10)
  )
  mock_model <- stats::lm(severity_score ~ day, data = mock_data)

  all_model_vars <- c(
    "b_Intercept",
    "b_day",
    "sd_id__Intercept",
    "sigma",
    "ar[1]"
  )

  queried_patterns <- character(0)

  local_mocked_bindings(
    get_variables = function(model, ...) {
      all_model_vars
    },
    gather_draws = function(model, ..., regex = FALSE) {
      dots <- rlang::enquos(...)
      pat <- rlang::as_name(dots[[1]])
      queried_patterns <<- c(queried_patterns, pat)
      matched_vars <- grep(pat, all_model_vars, value = TRUE)
      tidyr::expand_grid(
        .variable = matched_vars,
        .chain = 1:2,
        .iteration = 1:10
      ) |>
        dplyr::mutate(
          .draw = dplyr::row_number(),
          .value = stats::rnorm(dplyr::n())
        )
    },
    .package = "tidybayes"
  )

  p_flags <- get_traceplot(
    mock_model,
    fixed = TRUE,
    sd = FALSE,
    sigma = FALSE,
    ar = FALSE
  )
  expect_s3_class(p_flags, "patchwork")
  expect_identical(queried_patterns, "^b_")

  queried_patterns <- character(0)
  p_regex <- get_traceplot(
    mock_model,
    fixed = FALSE,
    sd = FALSE,
    sigma = FALSE,
    ar = FALSE,
    regex = "^(b_day|sigma)$"
  )
  expect_s3_class(p_regex, "patchwork")
  expect_identical(queried_patterns, "^(b_day|sigma)$")
})

test_that("get_traceplot() validates arguments and throws informative cli errors", {
  mock_data <- tibble::tibble(
    day = 1:10,
    severity_score = stats::rnorm(10)
  )
  valid_model <- stats::lm(severity_score ~ day, data = mock_data)

  local_mocked_bindings(
    get_variables = function(model, ...) {
      c("b_Intercept", "sigma")
    },
    .package = "tidybayes"
  )

  expect_error(
    get_traceplot("not_a_model"),
    class = "peoplenotp_error_bad_model"
  )

  expect_error(
    get_traceplot(valid_model, fixed = "yes"),
    class = "peoplenotp_error_bad_flag"
  )

  expect_error(
    get_traceplot(valid_model, sd = NA),
    class = "peoplenotp_error_bad_flag"
  )

  expect_error(
    get_traceplot(
      valid_model,
      fixed = FALSE,
      sd = FALSE,
      sigma = FALSE,
      ar = FALSE
    ),
    class = "peoplenotp_error_no_parameter_groups"
  )

  expect_error(
    get_traceplot(valid_model, regex = ""),
    class = "peoplenotp_error_bad_regex"
  )

  expect_error(
    get_traceplot(valid_model, alpha = -0.1),
    class = "peoplenotp_error_bad_alpha"
  )

  expect_error(
    get_traceplot(valid_model, linewidth = 0),
    class = "peoplenotp_error_bad_linewidth"
  )

  expect_error(
    get_traceplot(valid_model, scales = "invalid_scale"),
    class = "peoplenotp_error_bad_scales"
  )

  expect_error(
    get_traceplot(valid_model, facet_ncol = 0),
    class = "peoplenotp_error_bad_facet_ncol"
  )

  expect_error(
    get_traceplot(valid_model, patchwork_ncol = 1.5),
    class = "peoplenotp_error_bad_patchwork_ncol"
  )

  expect_error(
    get_traceplot(valid_model, regex = "^nonexistent_param_"),
    class = "peoplenotp_error_no_matching_variables"
  )
})
