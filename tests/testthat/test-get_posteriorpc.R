test_that("get_posteriorpc() builds a ggplot with density layers for standard models", {
  set.seed(123)
  mock_data <- tibble::tibble(
    day = 1:20,
    severity_score = stats::rnorm(20, mean = 30, sd = 5)
  )

  mock_model <- stats::lm(severity_score ~ day, data = mock_data)

  captured_newdata <- NULL

  local_mocked_bindings(
    add_predicted_draws = function(newdata, object, ..., ndraws = NULL) {
      captured_newdata <<- newdata
      n_draws <- ndraws %||% 10L
      tidyr::expand_grid(
        newdata,
        .draw = seq_len(n_draws)
      ) |>
        dplyr::mutate(
          .prediction = stats::rnorm(dplyr::n(), mean = 30, sd = 5)
        )
    },
    .package = "tidybayes"
  )

  p <- get_posteriorpc(
    mock_model,
    ndraws = 15,
    alpha = 0.3,
    accent_color = "#E76F51",
    linewidth = 2.5
  )

  expect_s3_class(p, "ggplot")
  expect_true("severity_score" %in% names(captured_newdata))
  expect_length(p$layers, 2L)
  expect_equal(p$layers[[1]]$aes_params$alpha, 0.3)
  expect_equal(p$layers[[2]]$aes_params$colour, "#E76F51")
  expect_equal(p$layers[[2]]$aes_params$linewidth, 2.5)
})

test_that("get_posteriorpc() removes outcome and informs when ar() is present unless keep_outcome = TRUE", {
  set.seed(42)
  mock_data <- tibble::tibble(
    day = 1:20,
    severity_score = stats::rnorm(20, mean = 30, sd = 5)
  )

  mock_ar_model <- stats::lm(severity_score ~ day, data = mock_data)

  captured_newdata <- NULL

  local_mocked_bindings(
    find_formula = function(x, ...) {
      list(conditional = severity_score ~ day + ar(time = day))
    },
    .package = "insight"
  )

  local_mocked_bindings(
    add_predicted_draws = function(newdata, object, ..., ndraws = NULL) {
      captured_newdata <<- newdata
      n_draws <- ndraws %||% 10L
      tidyr::expand_grid(
        newdata,
        .draw = seq_len(n_draws)
      ) |>
        dplyr::mutate(
          .prediction = stats::rnorm(dplyr::n(), mean = 30, sd = 5)
        )
    },
    .package = "tidybayes"
  )

  expect_message(
    p_dropped <- get_posteriorpc(mock_ar_model, ndraws = 10),
    class = "peoplenotp_message_ar_outcome_removed"
  )

  expect_s3_class(p_dropped, "ggplot")
  expect_false("severity_score" %in% names(captured_newdata))
  expect_true("severity_score" %in% names(p_dropped$layers[[2]]$data))

  expect_no_message(
    p_kept <- get_posteriorpc(mock_ar_model, ndraws = 10, keep_outcome = TRUE)
  )
  expect_s3_class(p_kept, "ggplot")
  expect_true("severity_score" %in% names(captured_newdata))
})

test_that("get_posteriorpc() validates arguments with informative cli errors", {
  mock_data <- tibble::tibble(
    day = 1:10,
    severity_score = stats::rnorm(10)
  )
  valid_model <- stats::lm(severity_score ~ day, data = mock_data)

  expect_error(
    get_posteriorpc("not_a_model"),
    class = "peoplenotp_error_bad_model"
  )

  expect_error(
    get_posteriorpc(list(data = "not_a_df")),
    class = "peoplenotp_error_bad_model"
  )

  expect_error(
    get_posteriorpc(valid_model, ndraws = -5),
    class = "peoplenotp_error_bad_ndraws"
  )

  expect_error(
    get_posteriorpc(valid_model, ndraws = 2.5),
    class = "peoplenotp_error_bad_ndraws"
  )

  expect_error(
    get_posteriorpc(valid_model, alpha = 1.5),
    class = "peoplenotp_error_bad_alpha"
  )

  expect_error(
    get_posteriorpc(valid_model, accent_color = ""),
    class = "peoplenotp_error_bad_accent_color"
  )

  expect_error(
    get_posteriorpc(valid_model, linewidth = 0),
    class = "peoplenotp_error_bad_linewidth"
  )

  expect_error(
    get_posteriorpc(valid_model, keep_outcome = NA),
    class = "peoplenotp_error_bad_keep_outcome"
  )
})
