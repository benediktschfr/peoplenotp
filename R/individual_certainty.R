#' Evaluate individual-level Bayesian clinical certainty across MID thresholds
#'
#' @description
#' Summarizes individual-level posterior draws across one or more Minimal
#' Important Difference (`mid`) thresholds. Computes posterior medians, credible
#' intervals, Probability of Direction (`pd`), full directional ROPE
#' partitioning (`rope_pct`, `improvement_pct`, `deterioration_pct`),
#' `directed_certainty`, and ordered Makowski evidence tiers (`makowski_tier`),
#' and tabulates an inverted Clinical Significance Evidence Matrix.
#'
#' @param data A tidy data frame or tibble containing posterior draws per
#'   individual.
#' @param id <[`data-masking`][rlang::args_data_masking]> Unquoted column name
#'   uniquely identifying individuals.
#' @param value <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name containing the individual target estimate per draw (e.g., `.epred`,
#'   `net_change`).
#' @param mid Numeric scalar (e.g., `2.5`) or numeric vector of positive values
#'   (e.g., `seq(1.5, 3.5, by = 0.5)`) defining the Minimal Important
#'   Difference threshold(s).
#' @param better_is Character string indicating clinical polarity (`"lower"` or
#'   `"higher"`). Defaults to `"lower"`.
#' @param draw_col <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name identifying the MCMC draw index (default `.draw`).
#' @param ci Credible interval width for individual summaries (default `0.95`).
#'
#' @return An S3 list object of class `"pnp_individual_certainty"` containing:
#'   * `$individual_summary`: A [tibble::tibble()] with one row per individual
#'     and `mid` value containing `<id>`, `mid`, `median`, `ci_low`, `ci_high`,
#'     `pd`, `rope_pct`, `improvement_pct`, `deterioration_pct`,
#'     `directed_certainty`, and `makowski_tier` (an ordered factor).
#'   * `$matrix_data`: A [tibble::tibble()] representing the inverted Makowski
#'     Evidence Matrix crossing `makowski_tier` (rows) by `mid` thresholds
#'     (columns).
#'   * `$meta`: A list of metadata (`mid`, `better_is`, `ci`, `id_col`,
#'     `value_col`, `draw_col`, `n_individuals`, and `matrix_long`).
#'
#' @examples
#' set.seed(123)
#' draws_df <- tibble::tibble(
#'   id = rep(paste0("P", 1:10), each = 200),
#'   .draw = rep(1:200, times = 10),
#'   net_change = stats::rnorm(2000, mean = rep(seq(-5, 2, length.out = 10), each = 200), sd = 1)
#' )
#'
#' res <- individual_certainty(
#'   data = draws_df,
#'   id = id,
#'   value = net_change,
#'   mid = c(2.0, 2.5, 3.0),
#'   better_is = "lower"
#' )
#'
#' print(res)
#' plot(res, mid_select = 2.5)
#'
#' @export
individual_certainty <- function(
  data,
  id,
  value,
  mid,
  better_is = c("lower", "higher"),
  draw_col = .draw,
  ci = 0.95
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

  if (missing(id)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg id} is required to identify individuals."
      ),
      class = "peoplenotp_error_missing_column"
    )
  }

  if (missing(value)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg value} is required to specify the posterior draw column."
      ),
      class = "peoplenotp_error_missing_column"
    )
  }

  if (
    missing(mid) ||
      !is.numeric(mid) ||
      length(mid) == 0L ||
      any(is.na(mid)) ||
      any(is.infinite(mid)) ||
      any(mid <= 0)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid} must be a positive numeric scalar or vector.",
        "i" = "Received {.val {if (missing(mid)) NULL else mid}}."
      ),
      class = "peoplenotp_error_bad_mid"
    )
  }

  if (identical(better_is, c("lower", "higher"))) {
    better_is <- "lower"
  }

  if (
    !is.character(better_is) ||
      length(better_is) != 1L ||
      is.na(better_is) ||
      !better_is %in% c("lower", "higher")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg better_is} must be either {.val lower} or {.val higher}.",
        "i" = "Received {.val {better_is}}."
      ),
      class = "peoplenotp_error_bad_better_is"
    )
  }

  if (
    !is.numeric(ci) ||
      length(ci) != 1L ||
      is.na(ci) ||
      ci <= 0 ||
      ci >= 1
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg ci} must be a single numeric value strictly between 0 and 1.",
        "i" = "Received {.val {ci}}."
      ),
      class = c("peoplenotp_error_bad_ci", "peoplenotp_error_bad_ci_width")
    )
  }

  ungrouped_data <- dplyr::ungroup(data)

  id_quo <- rlang::enquo(id)
  value_quo <- rlang::enquo(value)
  draw_quo <- rlang::enquo(draw_col)

  id_name <- resolve_column_name(ungrouped_data, id_quo)
  value_name <- resolve_column_name(ungrouped_data, value_quo)

  if (!missing(draw_col)) {
    draw_name <- resolve_column_name(ungrouped_data, draw_quo)
  } else {
    draw_label <- rlang::as_label(draw_quo)
    draw_name <- if (draw_label %in% names(ungrouped_data)) draw_label else NULL
  }

  val_vec <- ungrouped_data[[value_name]]
  if (!is.numeric(val_vec)) {
    cli::cli_abort(
      c(
        "x" = "Column {.code {value_name}} must be numeric.",
        "i" = "Column has class {.cls {class(val_vec)}}."
      ),
      class = "peoplenotp_error_bad_column_type"
    )
  }

  mid_vals <- unique(mid)

  base_summary <- dplyr::summarise(
    ungrouped_data,
    compute_individual_base(
      draws = .data[[value_name]],
      ci = ci
    ),
    .by = dplyr::all_of(id_name)
  )

  mid_summaries <- lapply(mid_vals, function(m) {
    rope_df <- dplyr::summarise(
      ungrouped_data,
      compute_individual_mid(
        draws = .data[[value_name]],
        m = m,
        better_is = better_is
      ),
      .by = dplyr::all_of(id_name)
    )

    merged <- dplyr::left_join(base_summary, rope_df, by = id_name)
    dplyr::select(
      merged,
      dplyr::all_of(id_name),
      mid,
      median,
      ci_low,
      ci_high,
      pd,
      rope_pct,
      improvement_pct,
      deterioration_pct,
      directed_certainty,
      makowski_tier
    )
  })

  individual_summary <- tibble::as_tibble(dplyr::bind_rows(mid_summaries))
  n_individuals <- dplyr::n_distinct(ungrouped_data[[id_name]])

  matrix_long <- build_evidence_matrix_long(
    individual_summary = individual_summary,
    mid_vals = mid_vals,
    n_individuals = n_individuals
  )

  matrix_data <- format_evidence_matrix_wide(
    matrix_long = matrix_long,
    mid_vals = mid_vals,
    format = "both",
    digits = 1
  )
  attr(matrix_data, "long") <- matrix_long

  meta <- list(
    mid = mid_vals,
    better_is = better_is,
    ci = ci,
    id_col = id_name,
    value_col = value_name,
    draw_col = draw_name,
    n_individuals = n_individuals,
    matrix_long = matrix_long
  )

  structure(
    list(
      individual_summary = individual_summary,
      matrix_data = matrix_data,
      meta = meta
    ),
    class = "pnp_individual_certainty"
  )
}

#' Print method for `pnp_individual_certainty` objects
#'
#' @param x An object of class `"pnp_individual_certainty"`.
#' @param format Character string specifying cell format: `"both"` (default,
#'   `n (%)`), `"count"` (`n`), or `"prop"` (`%`).
#' @param digits Number of decimal places for percentages (default `1`).
#' @param ... Additional arguments passed to [insight::export_table()].
#'
#' @return Invisibly returns `x`.
#' @export
print.pnp_individual_certainty <- function(
  x,
  format = c("both", "count", "prop"),
  digits = 1,
  ...
) {
  if (identical(format, c("both", "count", "prop"))) {
    format <- "both"
  }

  if (
    !is.character(format) ||
      length(format) != 1L ||
      is.na(format) ||
      !format %in% c("both", "count", "prop")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg format} must be one of {.val both}, {.val count}, or {.val prop}.",
        "i" = "Received {.val {format}}."
      ),
      class = "peoplenotp_error_bad_format"
    )
  }

  wide_df <- format_evidence_matrix_wide(
    matrix_long = x$meta$matrix_long,
    mid_vals = x$meta$mid,
    format = format,
    digits = digits
  )

  display_df <- as.data.frame(
    wide_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  display_df[[1]] <- as.character(display_df[[1]])
  names(display_df)[1] <- "Evidence Tier"

  out <- insight::export_table(
    display_df,
    title = "# Clinical Significance Evidence Matrix",
    ...
  )

  cat(out, "\n", sep = "")
  invisible(x)
}

#' Summary method for `pnp_individual_certainty` objects
#'
#' @param object An object of class `"pnp_individual_certainty"`.
#' @param format Character string passed to [print.pnp_individual_certainty()].
#' @param digits Integer passed to [print.pnp_individual_certainty()].
#' @param ... Additional arguments passed to [print.pnp_individual_certainty()].
#'
#' @return Invisibly returns `object`.
#' @export
summary.pnp_individual_certainty <- function(
  object,
  format = c("both", "count", "prop"),
  digits = 1,
  ...
) {
  print.pnp_individual_certainty(object, format = format, digits = digits, ...)
}

#' Plot a Bayesian Clinical Volcano Plot for `pnp_individual_certainty` objects
#'
#' @description
#' Visualizes individual posterior medians and credible intervals on the x-axis
#' against `directed_certainty` (0% to 100%) on the y-axis, with vertical ROPE
#' shading via [geom_rope()] and automatic directional annotations derived from
#' `x$meta$better_is`.
#'
#' @param x An object of class `"pnp_individual_certainty"`.
#' @param mid_select Optional single numeric value selecting which `mid`
#'   threshold from `x$meta$mid` to plot. Defaults to `NULL` (uses the first
#'   `mid` value in `x$meta$mid`).
#' @param point_alpha Numeric in `[0, 1]` controlling point opacity (default
#'   `0.85`).
#' @param interval_alpha Numeric in `[0, 1]` controlling horizontal error bar
#'   opacity (default `0.5`).
#' @param ... Additional arguments passed to [ggplot2::geom_point()].
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.pnp_individual_certainty <- function(
  x,
  mid_select = NULL,
  point_alpha = 0.85,
  interval_alpha = 0.5,
  ...
) {
  if (
    !is.numeric(point_alpha) ||
      length(point_alpha) != 1L ||
      is.na(point_alpha) ||
      point_alpha < 0 ||
      point_alpha > 1
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg point_alpha} must be a single numeric value between 0 and 1.",
        "i" = "Received {.val {point_alpha}}."
      ),
      class = "peoplenotp_error_bad_alpha"
    )
  }

  if (
    !is.numeric(interval_alpha) ||
      length(interval_alpha) != 1L ||
      is.na(interval_alpha) ||
      interval_alpha < 0 ||
      interval_alpha > 1
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg interval_alpha} must be a single numeric value between 0 and 1.",
        "i" = "Received {.val {interval_alpha}}."
      ),
      class = "peoplenotp_error_bad_alpha"
    )
  }

  available_mids <- x$meta$mid
  if (is.null(mid_select)) {
    selected_mid <- available_mids[[1]]
  } else {
    if (
      !is.numeric(mid_select) ||
        length(mid_select) != 1L ||
        is.na(mid_select) ||
        !mid_select %in% available_mids
    ) {
      cli::cli_abort(
        c(
          "x" = "Argument {.arg mid_select} must be a single value present in `x$meta$mid` ({.val {available_mids}}).",
          "i" = "Received {.val {mid_select}}."
        ),
        class = "peoplenotp_error_bad_mid_select"
      )
    }
    selected_mid <- mid_select
  }

  plot_df <- dplyr::filter(x$individual_summary, mid == selected_mid)

  if (x$meta$better_is == "lower") {
    left_label <- "\u2190 Favors Improvement"
    right_label <- "Favors Deterioration \u2192"
  } else {
    left_label <- "\u2190 Favors Deterioration"
    right_label <- "Favors Improvement \u2192"
  }

  x_min_bound <- min(c(plot_df$ci_low, -selected_mid * 1.5), na.rm = TRUE)
  x_max_bound <- max(c(plot_df$ci_high, selected_mid * 1.5), na.rm = TRUE)
  left_x <- min(-selected_mid * 1.1, (x_min_bound - selected_mid) / 2)
  right_x <- max(selected_mid * 1.1, (x_max_bound + selected_mid) / 2)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = median, y = directed_certainty)
  ) +
    geom_rope(x_range = c(-selected_mid, selected_mid)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0,
      alpha = interval_alpha
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = makowski_tier),
      alpha = point_alpha,
      ...
    ) +
    ggplot2::annotate(
      "text",
      x = left_x,
      y = 1.03,
      label = left_label,
      hjust = 0.5,
      size = 3.5
    ) +
    ggplot2::annotate(
      "text",
      x = right_x,
      y = 1.03,
      label = right_label,
      hjust = 0.5,
      size = 3.5
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      limits = c(0, 1.08),
      breaks = seq(0, 1, by = 0.25)
    ) +
    ggplot2::labs(
      x = "Posterior Median (with Credible Interval)",
      y = "Directed Clinical Certainty",
      color = "Evidence Tier"
    ) +
    ggplot2::theme_minimal()
}

#' Ordered Makowski tier factor levels
#'
#' @return Character vector of length 5.
#' @noRd
makowski_tier_levels <- function() {
  c(
    "Negligible (>99% ROPE)",
    "Prob. Negligible (>97.5%)",
    "Undecided",
    "Prob. Significant (<2.5%)",
    "Significant (<1% ROPE)"
  )
}

#' Classify ROPE probabilities into ordered Makowski evidence tiers
#'
#' @param rope_pct Numeric vector of ROPE probabilities in `[0, 1]`.
#'
#' @return Ordered factor.
#' @noRd
assign_makowski_tier <- function(rope_pct) {
  tiers <- dplyr::case_when(
    rope_pct > 0.99 ~ "Negligible (>99% ROPE)",
    rope_pct > 0.975 ~ "Prob. Negligible (>97.5%)",
    rope_pct < 0.01 ~ "Significant (<1% ROPE)",
    rope_pct < 0.025 ~ "Prob. Significant (<2.5%)",
    TRUE ~ "Undecided"
  )

  factor(tiers, levels = makowski_tier_levels(), ordered = TRUE)
}

#' Resolve an unquoted or string column selector against a data frame
#'
#' @param data A data frame.
#' @param col_quo Quosure capturing the column argument.
#'
#' @return Single character string column name.
#' @noRd
resolve_column_name <- function(data, col_quo) {
  expr <- rlang::quo_get_expr(col_quo)
  col_name <- if (is.character(expr) && length(expr) == 1L) {
    expr
  } else if (is.symbol(expr)) {
    rlang::as_string(expr)
  } else {
    rlang::as_label(col_quo)
  }

  if (!col_name %in% names(data)) {
    cli::cli_abort(
      c(
        "x" = "Column {.code {col_name}} was not found in {.arg data}.",
        "i" = "Available columns: {.val {names(data)}}."
      ),
      class = "peoplenotp_error_missing_column"
    )
  }

  col_name
}

#' Compute MID-invariant individual summary metrics
#'
#' @param draws Numeric vector of posterior draws.
#' @param ci Credible interval width.
#'
#' @return A single-row tibble with `median`, `ci_low`, `ci_high`, and `pd`.
#' @noRd
compute_individual_base <- function(draws, ci) {
  qi <- tidybayes::median_qi(draws, .width = ci)
  pd_val <- max(
    mean(draws > 0, na.rm = TRUE),
    mean(draws < 0, na.rm = TRUE)
  )

  tibble::tibble(
    median = qi$y,
    ci_low = qi$ymin,
    ci_high = qi$ymax,
    pd = pd_val
  )
}

#' Compute MID-specific individual ROPE metrics
#'
#' @param draws Numeric vector of posterior draws.
#' @param m Positive numeric MID scalar.
#' @param better_is `"lower"` or `"higher"`.
#'
#' @return A single-row tibble.
#' @noRd
compute_individual_mid <- function(draws, m, better_is) {
  rope_pct <- mean(draws >= -m & draws <= m, na.rm = TRUE)

  if (better_is == "lower") {
    imp_pct <- mean(draws < -m, na.rm = TRUE)
    det_pct <- mean(draws > m, na.rm = TRUE)
  } else {
    imp_pct <- mean(draws > m, na.rm = TRUE)
    det_pct <- mean(draws < -m, na.rm = TRUE)
  }

  dir_cert <- pmax(imp_pct, det_pct)

  tibble::tibble(
    mid = m,
    rope_pct = rope_pct,
    improvement_pct = imp_pct,
    deterioration_pct = det_pct,
    directed_certainty = dir_cert,
    makowski_tier = assign_makowski_tier(rope_pct)
  )
}

#' Build long evidence matrix counts and proportions
#'
#' @param individual_summary Individual summary tibble.
#' @param mid_vals Numeric vector of MID values.
#' @param n_individuals Total number of individuals.
#'
#' @return A tibble with `makowski_tier`, `mid`, `n`, `prop`, `pct`.
#' @noRd
build_evidence_matrix_long <- function(
  individual_summary,
  mid_vals,
  n_individuals
) {
  all_tiers <- factor(
    makowski_tier_levels(),
    levels = makowski_tier_levels(),
    ordered = TRUE
  )

  grid <- tidyr::expand_grid(
    makowski_tier = all_tiers,
    mid = mid_vals
  )

  counts <- dplyr::summarise(
    individual_summary,
    n = dplyr::n(),
    .by = c(makowski_tier, mid)
  )

  dplyr::left_join(grid, counts, by = c("makowski_tier", "mid")) |>
    dplyr::mutate(
      n = dplyr::coalesce(n, 0L),
      prop = if (n_individuals > 0L) n / n_individuals else 0,
      pct = prop * 100
    )
}

#' Pivot long evidence matrix into wide inverted matrix format
#'
#' @param matrix_long Long evidence matrix tibble.
#' @param mid_vals Numeric vector of MID values.
#' @param format `"both"`, `"count"`, or `"prop"`.
#' @param digits Number of decimal digits for percentages.
#'
#' @return A wide tibble with `makowski_tier` rows and `MID = ...` columns.
#' @noRd
format_evidence_matrix_wide <- function(
  matrix_long,
  mid_vals,
  format = "both",
  digits = 1
) {
  mid_labels <- stats::setNames(
    paste0("MID = ", format(mid_vals, trim = TRUE, nsmall = 1)),
    as.character(mid_vals)
  )

  cells <- if (format == "both") {
    sprintf(paste0("%d (%.", digits, "f%%)"), matrix_long$n, matrix_long$pct)
  } else if (format == "count") {
    as.character(matrix_long$n)
  } else {
    sprintf(paste0("%.", digits, "f%%"), matrix_long$pct)
  }

  formatted_long <- dplyr::mutate(
    matrix_long,
    mid_col = unname(mid_labels[as.character(mid)]),
    cell = cells
  )

  tidyr::pivot_wider(
    dplyr::select(formatted_long, makowski_tier, mid_col, cell),
    names_from = mid_col,
    values_from = cell
  )
}
