#' Calculate joint clinical certainty for two outcomes
#'
#' @description
#' Calculates the joint probability of clinical success (or deterioration)
#' across two co-primary endpoints. By evaluating the condition row-wise (per
#' MCMC draw), it naturally accounts for the correlation between the two
#' posterior distributions. Calculates the joint Probability of Direction
#' (`joint_pd`), joint improvement percentage, and joint deterioration
#' percentage. Optionally, marginal statistics and clinical statuses can be
#' returned.
#'
#' Supports both ungrouped data frames (single participant) and grouped data
#' frames or `.by` specifications (multiple participants), allowing it to be
#' piped directly after `tidybayes` draw-generating functions.
#'
#' @param data A data frame or tibble containing the posterior draws. If
#'   grouped (e.g., via [dplyr::group_by()]) and `.by` is not supplied,
#'   existing grouping variables except `.row` (automatically added by
#'   `tidybayes::add_epred_draws()`) are preserved and used for group-wise
#'   summarization.
#' @param col_1 <[`data-masking`][rlang::args_data_masking]> The first numeric
#'   column in `data` containing posterior draws.
#' @param col_2 <[`data-masking`][rlang::args_data_masking]> The second numeric
#'   column in `data` containing posterior draws.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional grouping
#'   columns for summarizing across multiple participants or conditions (e.g.,
#'   `.by = id` or `.by = c(id, group)`). If provided, any existing grouping on
#'   `data` is ignored. Defaults to `NULL`.
#' @param better_is A character vector of length 2 indicating the desired
#'   direction for `col_1` and `col_2`. Allowed values are `"higher"` or
#'   `"lower"`. Defaults to `c("higher", "higher")`.
#' @param ranges A list of two numeric vectors, each of length 2, defining the
#'   lower and upper ROPE boundaries for `col_1` and `col_2` in the form
#'   `list(c(lower_1, upper_1), c(lower_2, upper_2))`. Defaults to
#'   `list(c(-1, 1), c(-1, 1))`.
#' @param ci_width A single numeric value strictly between 0 and 1 defining the
#'   width of the credible interval returned for the medians. Defaults to
#'   `0.95`.
#' @param certainty_boundaries A numeric vector of length 2 defining the
#'   probability thresholds (between 0 and 1) for assigning `"Likely"` and
#'   `"Significant"` categories, in the order
#'   `c(likely_boundary, significant_boundary)`. Defaults to `c(0.975, 0.99)`.
#' @param include_marginals A single logical flag. If `TRUE`, the output also
#'   includes marginal statistics (`pd`, `improvement_pct`, and
#'   `deterioration_pct`) for each individual outcome. Defaults to `FALSE`.
#' @param include_status A single logical flag indicating whether categorical
#'   clinical status classifications (`joint_clinical_status`, and marginal
#'   statuses if `include_marginals = TRUE`) should be included in the output.
#'   Defaults to `FALSE`.
#'
#' @return A [tibble::tibble()] with one row per group (or a single row if
#'   ungrouped) containing any grouping columns along with:
#'   * `median_1`, `lower_1`, `upper_1`: Posterior median and credible
#'     interval bounds for `col_1`.
#'   * `pd_1`, `improvement_pct_1`, `deterioration_pct_1`:
#'     Marginal Probability of Direction and ROPE probabilities for `col_1`
#'     (included when `include_marginals = TRUE`).
#'   * `clinical_status_1`: Marginal clinical status for `col_1` (included when
#'     `include_marginals = TRUE` and `include_status = TRUE`).
#'   * `median_2`, `lower_2`, `upper_2`: Posterior median and credible
#'     interval bounds for `col_2`.
#'   * `pd_2`, `improvement_pct_2`, `deterioration_pct_2`:
#'     Marginal Probability of Direction and ROPE probabilities for `col_2`
#'     (included when `include_marginals = TRUE`).
#'   * `clinical_status_2`: Marginal clinical status for `col_2` (included when
#'     `include_marginals = TRUE` and `include_status = TRUE`).
#'   * `joint_pd`: Joint Probability of Direction (proportion of draws where
#'     both outcomes simultaneously lie on their respective dominant sides).
#'   * `joint_improvement_pct`: Proportion of draws where both outcomes
#'     simultaneously exceed their ROPE in the beneficial direction.
#'   * `joint_deterioration_pct`: Proportion of draws where both outcomes
#'     simultaneously exceed their ROPE in the detrimental direction.
#'   * `joint_clinical_status`: Character classification of joint status
#'     (`"Significant Improvement"`, `"Likely Improvement"`,
#'     `"Significant Deterioration"`, `"Likely Deterioration"`, or
#'     `"Undecided"`, included when `include_status = TRUE`).
#'
#' @examples
#' # 1. Single participant workflow
#' set.seed(123)
#' df_patient <- tibble::tibble(
#'   wellbeing = stats::rnorm(4000, mean = 5, sd = 2),
#'   depression = wellbeing * -0.8 + stats::rnorm(4000, 0, 1)
#' )
#'
#' # Default: Medians, CrIs, and joint statistics
#' df_patient |>
#'   get_joint_clinical_certainty(
#'     col_1 = wellbeing,
#'     col_2 = depression,
#'     better_is = c("higher", "lower"),
#'     ranges = list(c(-2, 2), c(-1, 1))
#'   )
#'
#' # Extended: Include marginal statistics and clinical status
#' df_patient |>
#'   get_joint_clinical_certainty(
#'     col_1 = wellbeing,
#'     col_2 = depression,
#'     better_is = c("higher", "lower"),
#'     ranges = list(c(-2, 2), c(-1, 1)),
#'     include_marginals = TRUE,
#'     include_status = TRUE
#'   )
#'
#' # 2. Multiple participants using `.by` or `dplyr::group_by()`
#' df_multi <- tibble::tibble(
#'   id = rep(c("P1", "P2"), each = 1000),
#'   wellbeing = c(
#'     stats::rnorm(1000, mean = 5, sd = 1),
#'     stats::rnorm(1000, mean = 0, sd = 1)
#'   ),
#'   depression = c(
#'     stats::rnorm(1000, mean = -4, sd = 1),
#'     stats::rnorm(1000, mean = 0, sd = 1)
#'   )
#' )
#'
#' df_multi |>
#'   get_joint_clinical_certainty(
#'     col_1 = wellbeing,
#'     col_2 = depression,
#'     .by = id,
#'     better_is = c("higher", "lower"),
#'     ranges = list(c(-2, 2), c(-1, 1))
#'   )
#'
#' @export
get_joint_clinical_certainty <- function(
  data,
  col_1,
  col_2,
  .by = NULL,
  better_is = c("higher", "higher"),
  ranges = list(c(-1, 1), c(-1, 1)),
  ci_width = 0.95,
  certainty_boundaries = c(0.975, 0.99),
  include_marginals = FALSE,
  include_status = FALSE
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
      length(better_is) != 2L ||
      any(is.na(better_is)) ||
      !all(better_is %in% c("higher", "lower"))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg better_is} must be a character vector of length 2 containing only {.val higher} or {.val lower}.",
        "i" = "Received {.val {better_is}}."
      ),
      class = "peoplenotp_error_bad_better_is"
    )
  }

  valid_ranges <- is.list(ranges) &&
    length(ranges) == 2L &&
    all(vapply(
      ranges,
      \(r) is.numeric(r) && length(r) == 2L && !any(is.na(r)) && r[1] <= r[2],
      logical(1)
    ))

  if (!valid_ranges) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg ranges} must be a list of two numeric vectors of length 2, each with `range[1] <= range[2]`.",
        "i" = "Example: {.code list(c(-1, 1), c(-1, 1))}."
      ),
      class = c("peoplenotp_error_bad_ranges", "peoplenotp_error_bad_range")
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
    !is.logical(include_marginals) ||
      length(include_marginals) != 1L ||
      is.na(include_marginals)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg include_marginals} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {include_marginals}}."
      ),
      class = "peoplenotp_error_bad_include_marginals"
    )
  }

  if (
    !is.logical(include_status) ||
      length(include_status) != 1L ||
      is.na(include_status)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg include_status} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {include_status}}."
      ),
      class = "peoplenotp_error_bad_include_status"
    )
  }

  col_1_quo <- rlang::enquo(col_1)
  col_2_quo <- rlang::enquo(col_2)
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

  ungrouped_data <- dplyr::ungroup(data)

  for (col_quo in list(col_1_quo, col_2_quo)) {
    col_vals <- tryCatch(
      dplyr::pull(ungrouped_data, !!col_quo),
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
  }

  res <- dplyr::summarise(
    data,
    compute_joint_certainty_draws(
      draws1 = !!col_1_quo,
      draws2 = !!col_2_quo,
      better_is = better_is,
      ranges = ranges,
      ci_width = ci_width,
      certainty_boundaries = certainty_boundaries,
      include_marginals = include_marginals,
      include_status = include_status
    ),
    .by = !!by_quo
  )

  tibble::as_tibble(dplyr::ungroup(res))
}

#' Compute joint clinical certainty metrics for two vectors of posterior draws
#'
#' @param draws1 A numeric vector of posterior draws for outcome 1.
#' @param draws2 A numeric vector of posterior draws for outcome 2.
#' @param better_is Character vector of length 2: `"higher"` or `"lower"`.
#' @param ranges List of two numeric vectors of length 2 for ROPE bounds.
#' @param ci_width Numeric scalar for credible interval width.
#' @param certainty_boundaries Numeric vector of length 2 for thresholds.
#' @param include_marginals Logical scalar indicating whether to include
#'   marginal statistics.
#' @param include_status Logical scalar indicating whether to include
#'   clinical statuses.
#'
#' @return A single-row [tibble::tibble()] of summary metrics.
#' @noRd
compute_joint_certainty_draws <- function(
  draws1,
  draws2,
  better_is,
  ranges,
  ci_width,
  certainty_boundaries,
  include_marginals,
  include_status
) {
  pd_1_pos <- mean(draws1 > 0, na.rm = TRUE)
  pd_1_neg <- mean(draws1 < 0, na.rm = TRUE)
  pd_1 <- max(pd_1_pos, pd_1_neg)
  mask_pd_1 <- if (pd_1_pos >= pd_1_neg) draws1 > 0 else draws1 < 0

  pd_2_pos <- mean(draws2 > 0, na.rm = TRUE)
  pd_2_neg <- mean(draws2 < 0, na.rm = TRUE)
  pd_2 <- max(pd_2_pos, pd_2_neg)
  mask_pd_2 <- if (pd_2_pos >= pd_2_neg) draws2 > 0 else draws2 < 0

  joint_pd <- mean(mask_pd_1 & mask_pd_2, na.rm = TRUE)

  get_flags <- function(draws, direction, rope) {
    if (direction == "higher") {
      imp <- draws > rope[2]
      det <- draws < rope[1]
    } else {
      imp <- draws < rope[1]
      det <- draws > rope[2]
    }
    list(imp = imp, det = det)
  }

  flags1 <- get_flags(draws1, better_is[1], ranges[[1]])
  flags2 <- get_flags(draws2, better_is[2], ranges[[2]])

  imp_pct_1 <- mean(flags1$imp, na.rm = TRUE)
  det_pct_1 <- mean(flags1$det, na.rm = TRUE)

  imp_pct_2 <- mean(flags2$imp, na.rm = TRUE)
  det_pct_2 <- mean(flags2$det, na.rm = TRUE)

  joint_imp_pct <- mean(flags1$imp & flags2$imp, na.rm = TRUE)
  joint_det_pct <- mean(flags1$det & flags2$det, na.rm = TRUE)

  bound_likely <- certainty_boundaries[[1]]
  bound_signif <- certainty_boundaries[[2]]

  e1 <- tidybayes::median_qi(draws1, .width = ci_width)
  e2 <- tidybayes::median_qi(draws2, .width = ci_width)

  res_list <- list(
    median_1 = e1$y,
    lower_1 = e1$ymin,
    upper_1 = e1$ymax
  )

  if (include_marginals) {
    res_list$pd_1 <- pd_1
    res_list$improvement_pct_1 <- imp_pct_1
    res_list$deterioration_pct_1 <- det_pct_1
    if (include_status) {
      res_list$clinical_status_1 <- assign_clinical_status(
        imp_pct = imp_pct_1,
        det_pct = det_pct_1,
        bound_likely = bound_likely,
        bound_signif = bound_signif
      )
    }
  }

  res_list$median_2 <- e2$y
  res_list$lower_2 <- e2$ymin
  res_list$upper_2 <- e2$ymax

  if (include_marginals) {
    res_list$pd_2 <- pd_2
    res_list$improvement_pct_2 <- imp_pct_2
    res_list$deterioration_pct_2 <- det_pct_2
    if (include_status) {
      res_list$clinical_status_2 <- assign_clinical_status(
        imp_pct = imp_pct_2,
        det_pct = det_pct_2,
        bound_likely = bound_likely,
        bound_signif = bound_signif
      )
    }
  }

  res_list$joint_pd <- joint_pd
  res_list$joint_improvement_pct <- joint_imp_pct
  res_list$joint_deterioration_pct <- joint_det_pct

  if (include_status) {
    res_list$joint_clinical_status <- assign_clinical_status(
      imp_pct = joint_imp_pct,
      det_pct = joint_det_pct,
      bound_likely = bound_likely,
      bound_signif = bound_signif
    )
  }

  tibble::as_tibble(res_list)
}
