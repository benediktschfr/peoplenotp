#' Calculate clinical certainty based on ROPE
#'
#' @description
#' Evaluates posterior draws (such as differences between measurement points or
#' predictions from `tidybayes::add_epred_draws()`) and categorizes the effect
#' into clinical statuses. Calculates the median, credible interval, and
#' Probability of Direction (pd), and compares the posterior distribution
#' against a defined *Region of Practical Equivalence (ROPE)*.
#'
#' Supports both ungrouped data frames (single participant) and grouped data
#' frames or `.by` specifications (multiple participants), allowing it to be
#' piped directly after `tidybayes::add_epred_draws()`.
#'
#' @param data A data frame or tibble containing the posterior draws. If
#'   grouped (e.g., via [dplyr::group_by()]) and `.by` is not supplied,
#'   existing grouping variables except `.row` (automatically added by
#'   `tidybayes::add_epred_draws()`) are preserved and used for group-wise
#'   summarization.
#' @param col <[`data-masking`][rlang::args_data_masking]> The column in `data`
#'   containing the numeric posterior draws. Defaults to `.epred`.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional grouping
#'   columns for summarizing across multiple participants or conditions (e.g.,
#'   `.by = id` or `.by = c(id, group)`). If provided, any existing grouping on
#'   `data` is ignored. Defaults to `NULL`.
#' @param better_is A single character string indicating the desired direction
#'   of the effect. Allowed values are `"higher"` (higher values are better,
#'   e.g., well-being) or `"lower"` (lower values are better, e.g., pain).
#'   Defaults to `"higher"`.
#' @param range A numeric vector of length 2 defining the lower and upper
#'   boundaries of the ROPE (Minimum Clinically Important Difference), in the
#'   form `c(lower, upper)`. Defaults to `c(-1, 1)`.
#' @param ci_width A single numeric value strictly between 0 and 1 defining the
#'   width of the credible interval returned for the median effect. Defaults to
#'   `0.95`.
#' @param certainty_boundaries A numeric vector of length 2 defining the
#'   probability thresholds (between 0 and 1) for assigning clinical categories
#'   in the order `c(likely_boundary, significant_boundary)`. Defaults to
#'   `c(0.975, 0.99)`.
#' @param include_rope A single logical flag indicating whether the ROPE
#'   boundaries (`rope_low`, `rope_high`) and the proportion of draws within
#'   the ROPE (`rope_pct`) should be included in the output. Defaults to `TRUE`.
#'
#' @return A [tibble::tibble()] with one row per group (or a single row if
#'   ungrouped) containing any grouping columns along with:
#'   * `median`: Posterior median of the draws.
#'   * `lower`: Lower bound of the credible interval.
#'   * `upper`: Upper bound of the credible interval.
#'   * `pd`: Probability of Direction (maximum of the proportion of draws
#'     strictly greater than 0 or strictly less than 0).
#'   * `rope_low`, `rope_high`, `rope_pct`: Lower ROPE bound, upper ROPE bound,
#'     and proportion of draws within the ROPE (included when
#'     `include_rope = TRUE`).
#'   * `improvement_pct`: Proportion of draws exceeding the ROPE in the
#'     beneficial direction.
#'   * `deterioration_pct`: Proportion of draws exceeding the ROPE in the
#'     detrimental direction.
#'   * `clinical_status`: Character classification (`"Significant Improvement"`,
#'     `"Likely Improvement"`, `"Significant Deterioration"`,
#'     `"Likely Deterioration"`, or `"Undecided"`).
#'
#' @examples
#' # 1. Single participant workflow
#' set.seed(123)
#' df_patient <- tibble::tibble(
#'   iteration = 1:4000,
#'   score_difference = stats::rnorm(4000, mean = 5, sd = 1.5)
#' )
#'
#' df_patient |>
#'   get_clinical_certainty(
#'     col = score_difference,
#'     better_is = "higher",
#'     range = c(-2, 2)
#'   )
#'
#' # 2. Multiple participants using `.by` (or directly after add_epred_draws)
#' df_multi <- tibble::tibble(
#'   id = rep(c("P1", "P2", "P3"), each = 1000),
#'   .epred = c(
#'     stats::rnorm(1000, mean = 4.5, sd = 1.0),
#'     stats::rnorm(1000, mean = 0.2, sd = 0.8),
#'     stats::rnorm(1000, mean = -3.8, sd = 1.1)
#'   )
#' )
#'
#' df_multi |>
#'   get_clinical_certainty(
#'     .by = id,
#'     better_is = "higher",
#'     range = c(-2, 2)
#'   )
#'
#' @export
get_clinical_certainty <- function(
  data,
  col = .epred,
  .by = NULL,
  better_is = "higher",
  range = c(-1, 1),
  ci_width = 0.95,
  certainty_boundaries = c(0.975, 0.99),
  include_rope = TRUE
) {
  if (!is.data.frame(data)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg data} must be a data frame or tibble.",
        "i" = "Received an object of class {.cls {class(data)}}."
      ),
      class = "peoplenotp_error_bad_data"
    )
  }

  if (
    !is.character(better_is) ||
      length(better_is) != 1L ||
      is.na(better_is) ||
      !better_is %in% c("higher", "lower")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg better_is} must be either {.val higher} or {.val lower}.",
        "i" = "Received {.val {better_is}}."
      ),
      class = "peoplenotp_error_bad_better_is"
    )
  }

  if (
    !is.numeric(range) ||
      length(range) != 2L ||
      any(is.na(range)) ||
      range[1] > range[2]
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg range} must be a numeric vector of length 2 with `range[1] <= range[2]`.",
        "i" = "Received {.val {range}}."
      ),
      class = "peoplenotp_error_bad_range"
    )
  }

  if (
    !is.numeric(ci_width) ||
      length(ci_width) != 1L ||
      is.na(ci_width) ||
      ci_width <= 0 ||
      ci_width >= 1
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg ci_width} must be a single numeric value strictly between 0 and 1.",
        "i" = "Received {.val {ci_width}}."
      ),
      class = "peoplenotp_error_bad_ci_width"
    )
  }

  if (
    !is.numeric(certainty_boundaries) ||
      length(certainty_boundaries) != 2L ||
      any(is.na(certainty_boundaries)) ||
      any(certainty_boundaries <= 0) ||
      any(certainty_boundaries >= 1) ||
      certainty_boundaries[1] > certainty_boundaries[2]
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg certainty_boundaries} must be a numeric vector of length 2 between 0 and 1 with `certainty_boundaries[1] <= certainty_boundaries[2]`.",
        "i" = "Received {.val {certainty_boundaries}}."
      ),
      class = "peoplenotp_error_bad_certainty_boundaries"
    )
  }

  if (
    !is.logical(include_rope) ||
      length(include_rope) != 1L ||
      is.na(include_rope)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg include_rope} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {include_rope}}."
      ),
      class = "peoplenotp_error_bad_include_rope"
    )
  }

  col_quo <- rlang::enquo(col)
  by_quo <- rlang::enquo(.by)

  if (!rlang::quo_is_null(by_quo)) {
    data <- dplyr::ungroup(data)
  } else if (dplyr::is_grouped_df(data)) {
    grp_vars <- dplyr::group_vars(data)
    if (".row" %in% grp_vars && length(grp_vars) > 1L) {
      keep_grps <- setdiff(grp_vars, ".row")
      data <- data |>
        dplyr::ungroup() |>
        dplyr::group_by(dplyr::across(dplyr::all_of(keep_grps)))
    }
  }

  col_vals <- tryCatch(
    dplyr::pull(dplyr::ungroup(data), !!col_quo),
    error = function(cnd) {
      cli::cli_abort(
        c(
          "x" = "Column {.code {rlang::as_label(col_quo)}} was not found in {.arg data}.",
          "i" = "Available columns: {.val {names(data)}}."
        ),
        class = "peoplenotp_error_missing_column",
        parent = cnd
      )
    }
  )

  if (!is.numeric(col_vals)) {
    cli::cli_abort(
      c(
        "x" = "Column {.code {rlang::as_label(col_quo)}} must be numeric.",
        "i" = "Column has class {.cls {class(col_vals)}}."
      ),
      class = "peoplenotp_error_bad_column_type"
    )
  }

  res <- dplyr::summarise(
    data,
    compute_certainty_draws(
      posterior_draws = !!col_quo,
      better_is = better_is,
      range = range,
      ci_width = ci_width,
      certainty_boundaries = certainty_boundaries,
      include_rope = include_rope
    ),
    .by = !!by_quo
  )

  tibble::as_tibble(dplyr::ungroup(res))
}

#' Compute clinical certainty metrics for a single vector of posterior draws
#'
#' @param posterior_draws A numeric vector of posterior draws.
#' @param better_is Character string: `"higher"` or `"lower"`.
#' @param range Numeric vector of length 2 for ROPE bounds.
#' @param ci_width Numeric scalar for credible interval width.
#' @param certainty_boundaries Numeric vector of length 2 for thresholds.
#' @param include_rope Logical scalar indicating whether to include ROPE columns.
#'
#' @return A single-row [tibble::tibble()] of summary metrics.
#' @noRd
compute_certainty_draws <- function(
  posterior_draws,
  better_is,
  range,
  ci_width,
  certainty_boundaries,
  include_rope
) {
  pd_val <- max(
    mean(posterior_draws > 0, na.rm = TRUE),
    mean(posterior_draws < 0, na.rm = TRUE)
  )

  if (better_is == "higher") {
    imp_pct <- mean(posterior_draws > range[2], na.rm = TRUE)
    det_pct <- mean(posterior_draws < range[1], na.rm = TRUE)
  } else {
    imp_pct <- mean(posterior_draws < range[1], na.rm = TRUE)
    det_pct <- mean(posterior_draws > range[2], na.rm = TRUE)
  }

  status <- assign_clinical_status(
    imp_pct = imp_pct,
    det_pct = det_pct,
    bound_likely = certainty_boundaries[[1]],
    bound_signif = certainty_boundaries[[2]]
  )

  effect <- tidybayes::median_qi(posterior_draws, .width = ci_width)

  if (include_rope) {
    rope_pct <- mean(
      posterior_draws >= range[1] & posterior_draws <= range[2],
      na.rm = TRUE
    )

    tibble::tibble(
      median = effect$y,
      lower = effect$ymin,
      upper = effect$ymax,
      pd = pd_val,
      rope_low = range[1],
      rope_high = range[2],
      rope_pct = rope_pct,
      improvement_pct = imp_pct,
      deterioration_pct = det_pct,
      clinical_status = status
    )
  } else {
    tibble::tibble(
      median = effect$y,
      lower = effect$ymin,
      upper = effect$ymax,
      pd = pd_val,
      improvement_pct = imp_pct,
      deterioration_pct = det_pct,
      clinical_status = status
    )
  }
}
