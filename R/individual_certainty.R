#' Evaluate individual-level Bayesian clinical certainty across MID thresholds
#'
#' @description
#' Summarizes individual-level posterior draws across one or more Minimal
#' Important Difference (`mid`) thresholds. Computes posterior medians, credible
#' intervals, Probability of Direction (`pd`), full directional ROPE
#' partitioning (`rope_pct`, `improvement_pct`, `deterioration_pct`),
#' `directed_certainty`, `signed_certainty`, and ordered Makowski evidence tiers
#' (`makowski_tier`), and tabulates an inverted Clinical Significance Evidence
#' Matrix (optionally stratified by `.by` grouping variables such as therapy
#' conditions).
#'
#' @param data A tidy data frame or tibble containing posterior draws per
#'   individual.
#' @param id <[`data-masking`][rlang::args_data_masking]> Unquoted column name
#'   uniquely identifying individuals (default `id`), or `NULL` for single-case
#'   ($N = 1$) datasets without an ID column.
#' @param value <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name containing the individual target estimate per draw (e.g., `.epred`,
#'   `net_change`). Defaults to `.epred`.
#' @param mid Minimal Important Difference threshold(s). Can be supplied as:
#'   * A positive numeric scalar (e.g., `2.5`), defining a symmetric ROPE
#'     `c(-2.5, 2.5)`.
#'   * A numeric vector of length 2 with `mid[1] < 0` and `mid[2] > 0` (e.g.,
#'     `c(-10, 5)`), defining a single asymmetric ROPE range.
#'   * A numeric vector of positive values (e.g., `seq(1.5, 3.5, by = 0.5)`),
#'     defining multiple symmetric MID thresholds for sensitivity analyses.
#'   * A list of MID thresholds (e.g., `list(c(-10, 5), c(-12, 6))`), where each
#'     element is either a positive scalar or a numeric vector of length 2 with
#'     `mid[1] < 0 < mid[2]`, defining multiple MID steps (symmetric or
#'     asymmetric).
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional additional
#'   grouping column(s) such as therapy condition (e.g., `.by = group` or
#'   `.by = c(group, site)`). When provided (or when `data` is grouped via
#'   [dplyr::group_by()]), `$individual_summary`, `$matrix_data`, `print()`, and
#'   `plot()` are stratified across these groups. Defaults to `NULL`.
#' @param better_is Character string indicating clinical polarity (`"lower"` or
#'   `"higher"`). Defaults to `"lower"`.
#' @param draw_col <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name identifying the MCMC draw index (default `.draw`).
#' @param ci Credible interval width for individual summaries (default `0.95`).
#'
#' @return An S3 list object of class `"pnp_individual_certainty"` containing:
#'   * `$individual_summary`: A [tibble::tibble()] with one row per individual
#'     and `mid` step containing `<id>`, any `.by` grouping columns, `mid`,
#'     `rope_low`, `rope_high`, `median`, `ci_low`, `ci_high`, `pd`, `rope_pct`,
#'     `improvement_pct`, `deterioration_pct`, `directed_certainty`,
#'     `signed_certainty`, and `makowski_tier` (an ordered factor).
#'   * `$matrix_data`: A [tibble::tibble()] representing the inverted Makowski
#'     Evidence Matrix crossing `makowski_tier` (rows) by `mid` thresholds
#'     (columns), stratified by `.by` if supplied.
#'   * `$meta`: A list of metadata (`mid`, `mid_bounds`, `better_is`, `ci`,
#'     `id_col`, `group_cols`, `value_col`, `draw_col`, `n_individuals`, and
#'     `matrix_long`).
#'
#' @examples
#' set.seed(123)
#' draws_df <- tibble::tibble(
#'   id = rep(paste0("P", 1:10), each = 200),
#'   group = rep(rep(c("Waitlist", "Intervention"), each = 5), each = 200),
#'   .draw = rep(1:200, times = 10),
#'   net_change = stats::rnorm(2000, mean = rep(seq(-5, 2, length.out = 10), each = 200), sd = 1)
#' )
#'
#' res <- individual_certainty(
#'   data = draws_df,
#'   id = id,
#'   value = net_change,
#'   mid = c(2.0, 2.5, 3.0),
#'   .by = group,
#'   better_is = "lower"
#' )
#'
#' print(res)
#' plot(res, mid_select = 2.5)
#'
#' @export
individual_certainty <- function(
  data,
  id = id,
  value = .epred,
  mid,
  .by = NULL,
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

  mid_info <- resolve_individual_mid(mid)
  mid_vals <- mid_info$mid_vals
  mid_bounds <- mid_info$mid_bounds

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
  by_quo <- rlang::enquo(.by)
  draw_quo <- rlang::enquo(draw_col)

  id_name <- resolve_id_column(
    data = ungrouped_data,
    id_quo = id_quo,
    id_missing = missing(id)
  )
  value_name <- resolve_column_name(ungrouped_data, value_quo)

  group_cols <- resolve_group_columns(
    data = data,
    ungrouped_data = ungrouped_data,
    by_quo = by_quo,
    id_name = id_name
  )

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

  summary_by <- c(id_name, group_cols)

  base_summary <- dplyr::summarise(
    ungrouped_data,
    compute_individual_base(
      draws = .data[[value_name]],
      ci = ci
    ),
    .by = dplyr::all_of(summary_by)
  )

  mid_summaries <- lapply(seq_len(nrow(mid_bounds)), function(step_idx) {
    m_val <- mid_bounds$mid[[step_idx]]
    r_low <- mid_bounds$rope_low[[step_idx]]
    r_high <- mid_bounds$rope_high[[step_idx]]

    rope_df <- dplyr::summarise(
      ungrouped_data,
      compute_individual_mid(
        draws = .data[[value_name]],
        m = m_val,
        rope_low = r_low,
        rope_high = r_high,
        better_is = better_is
      ),
      .by = dplyr::all_of(summary_by)
    )

    merged <- if (length(summary_by) == 0L) {
      dplyr::bind_cols(base_summary, rope_df)
    } else {
      dplyr::left_join(base_summary, rope_df, by = summary_by)
    }
    dplyr::select(
      merged,
      dplyr::all_of(summary_by),
      mid,
      rope_low,
      rope_high,
      median,
      ci_low,
      ci_high,
      pd,
      rope_pct,
      improvement_pct,
      deterioration_pct,
      directed_certainty,
      signed_certainty,
      makowski_tier
    )
  })

  individual_summary <- tibble::as_tibble(dplyr::bind_rows(mid_summaries))
  n_individuals <- nrow(base_summary)

  matrix_long <- build_evidence_matrix_long(
    individual_summary = individual_summary,
    mid_vals = mid_vals,
    group_cols = group_cols
  )

  matrix_data <- format_evidence_matrix_wide(
    matrix_long = matrix_long,
    mid_vals = mid_vals,
    group_cols = group_cols,
    format = "both",
    digits = 1
  )
  attr(matrix_data, "long") <- matrix_long

  meta <- list(
    mid = mid_vals,
    mid_bounds = mid_bounds,
    better_is = better_is,
    ci = ci,
    id_col = id_name,
    group_cols = group_cols,
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
#' @param mid_select Optional numeric scalar or vector selecting which `mid`
#'   threshold(s) from `x$meta$mid` to display as columns in the printed
#'   Evidence Matrix. When `NULL` (default), displays the single `mid` if only
#'   one was computed, both `min` and `max` if two were computed, or the `min`,
#'   `middle`, and `max` values of the `mid` sequence if three or more were
#'   computed.
#' @param ... Additional arguments passed to [insight::export_table()].
#'
#' @return Invisibly returns `x`.
#' @export
print.pnp_individual_certainty <- function(
  x,
  format = c("both", "count", "prop"),
  digits = 1,
  mid_select = NULL,
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

  selected_mids <- resolve_print_mids(
    available_mids = x$meta$mid,
    mid_select = mid_select,
    mid_bounds = x$meta$mid_bounds
  )

  group_cols <- x$meta$group_cols %||% character(0)
  filtered_long <- dplyr::filter(
    x$meta$matrix_long,
    mid %in% selected_mids
  )

  wide_df <- format_evidence_matrix_wide(
    matrix_long = filtered_long,
    mid_vals = selected_mids,
    group_cols = group_cols,
    format = format,
    digits = digits
  )

  if (length(group_cols) == 0L) {
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
  } else {
    distinct_grps <- dplyr::distinct(
      wide_df,
      dplyr::across(dplyr::all_of(group_cols))
    )
    table_list <- vector("list", nrow(distinct_grps))
    title_list <- vector("list", nrow(distinct_grps))

    for (i in seq_len(nrow(distinct_grps))) {
      grp_row <- distinct_grps[i, , drop = FALSE]
      sub_wide <- dplyr::inner_join(wide_df, grp_row, by = group_cols) |>
        dplyr::select(-dplyr::all_of(group_cols))

      sub_df <- as.data.frame(
        sub_wide,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      sub_df[[1]] <- as.character(sub_df[[1]])
      names(sub_df)[1] <- "Evidence Tier"
      rownames(sub_df) <- NULL

      grp_label <- paste(
        vapply(
          group_cols,
          function(col) paste0(col, " = ", as.character(grp_row[[col]])),
          character(1)
        ),
        collapse = ", "
      )

      table_list[[i]] <- sub_df
      title_list[[i]] <- paste0(
        "# Clinical Significance Evidence Matrix (",
        grp_label,
        ")"
      )
    }

    out <- insight::export_table(
      table_list,
      title = title_list,
      ...
    )
  }

  cat(out, "\n", sep = "")
  invisible(x)
}

#' Summary method for `pnp_individual_certainty` objects
#'
#' @param object An object of class `"pnp_individual_certainty"`.
#' @param format Character string passed to [print.pnp_individual_certainty()].
#' @param digits Integer passed to [print.pnp_individual_certainty()].
#' @param mid_select Optional numeric scalar or vector passed to
#'   [print.pnp_individual_certainty()].
#' @param ... Additional arguments passed to [print.pnp_individual_certainty()].
#'
#' @return Invisibly returns `object`.
#' @export
summary.pnp_individual_certainty <- function(
  object,
  format = c("both", "count", "prop"),
  digits = 1,
  mid_select = NULL,
  ...
) {
  print.pnp_individual_certainty(
    object,
    format = format,
    digits = digits,
    mid_select = mid_select,
    ...
  )
}

#' Plot a Bayesian Clinical Volcano or Sensitivity Plot for `pnp_individual_certainty` objects
#'
#' @description
#' Visualizes individual-level Bayesian clinical certainty as either:
#' * **`type = "volcano"` (default):** Individual posterior medians and credible
#'   intervals on the x-axis against `directed_certainty` (0% to 100%) on the
#'   y-axis, with vertical ROPE shading via [geom_rope()] and automatic
#'   directional annotations derived from `x$meta$better_is`.
#' * **`type = "sensitivity"`:** A 2D clinical certainty sensitivity heatmap
#'   (`geom_tile()`) across `mid` values on the x-axis and a fine sequence of
#'   probability `thresholds` on the y-axis, where tile fill represents the
#'   proportion of individuals meeting clinical certainty. By default
#'   (`outcome = "both"`), renders a mirrored butterfly heatmap with
#'   improvement thresholds extending upward (`50%` to `100%`) and deterioration
#'   thresholds extending downward (`50%` to `100%`), or a single-sided heatmap
#'   when `outcome = "improvement"` or `"deterioration"`.
#'
#' When `.by` grouping variables were supplied to [individual_certainty()], both
#' plot types are automatically faceted across those groups (unless
#' `facet = FALSE`).
#'
#' @param x An object of class `"pnp_individual_certainty"`.
#' @param type Character string selecting the visualization type: `"volcano"`
#'   (default) or `"sensitivity"`.
#' @param mid_select Optional numeric value(s) selecting which `mid`
#'   threshold(s) from `x$meta$mid` to plot. For `type = "volcano"`, must be a
#'   single value (defaults to the first value in `x$meta$mid`). For
#'   `type = "sensitivity"`, can be a numeric vector to subset the plotted `mid`
#'   sequence (defaults to `NULL`, using all `x$meta$mid` values).
#' @param color_by Variable used to color both points and error bars when
#'   `type = "volcano"`. Can be supplied as a character string or unquoted name.
#'   Must be one of `"signed_certainty"` (default, ranging from `-1` for 100%
#'   deterioration to `+1` for 100% improvement), `"pd"`, `"makowski_tier"`,
#'   `"rope_pct"`, `"improvement_pct"`, `"deterioration_pct"`, or
#'   `"directed_certainty"`.
#' @param thresholds Numeric vector of probability thresholds in `(0, 1]` mapped
#'   to the y-axis when `type = "sensitivity"`. Defaults to
#'   `seq(0.50, 0.99, by = 0.01)`.
#' @param outcome Character string specifying which clinical direction(s) to
#'   display when `type = "sensitivity"`: `"both"` (default, mirrored butterfly
#'   heatmap with improvement above the center line and deterioration below),
#'   `"improvement"`, or `"deterioration"`.
#' @param color_palette Character string specifying the fill palette when
#'   `type = "sensitivity"`. Defaults to `"Temps"` (via `{rcartocolor}`). Also
#'   supports any `{rcartocolor}` palette, any `{viridisLite}` option
#'   (`"viridis"`, `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`, `"rocket"`,
#'   `"mako"`, `"turbo"`), or any palette name in [custom_palettes].
#' @param show_ci Logical flag (`TRUE` by default). For `type = "volcano"`,
#'   draws horizontal credible interval error bars (`ci_low` to `ci_high`) when
#'   `TRUE`.
#' @param point_size Numeric controlling point size in the volcano plot (default
#'   `2`).
#' @param point_alpha Numeric in `[0, 1]` controlling point/tile opacity
#'   (default `1`).
#' @param interval_alpha Numeric in `[0, 1]` controlling horizontal error bar
#'   opacity in the volcano plot (default `0.5`).
#' @param x_lab Optional string controlling the x-axis label. When `NULL`
#'   (default), uses `"Posterior Median (with Credible Interval)"` (or
#'   `"Posterior Median"` when `show_ci = FALSE`) for `type = "volcano"`, and
#'   `"Minimal Important Difference (MID)"` for `type = "sensitivity"`.
#' @param y_lab Optional string controlling the y-axis label. When `NULL`
#'   (default), uses `"Directed Clinical Certainty"` for `type = "volcano"`, and
#'   `"Certainty Threshold"` for `type = "sensitivity"`.
#' @param facet Logical flag (`TRUE` by default). When `TRUE` and `x` contains
#'   `.by` grouping columns in `x$meta$group_cols`, facets the plot via
#'   [ggplot2::facet_wrap()].
#' @param ... Additional arguments passed to [ggplot2::geom_point()] (for
#'   `type = "volcano"`) or [ggplot2::geom_tile()] (for `type = "sensitivity"`).
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.pnp_individual_certainty <- function(
  x,
  type = c("volcano", "sensitivity"),
  mid_select = NULL,
  color_by = c(
    "signed_certainty",
    "pd",
    "makowski_tier",
    "rope_pct",
    "improvement_pct",
    "deterioration_pct",
    "directed_certainty"
  ),
  thresholds = seq(0.50, 0.99, by = 0.01),
  outcome = c("both", "improvement", "deterioration"),
  color_palette = "Temps",
  show_ci = TRUE,
  point_size = 2,
  point_alpha = 1,
  interval_alpha = 0.5,
  x_lab = NULL,
  y_lab = NULL,
  facet = TRUE,
  ...
) {
  if (identical(type, c("volcano", "sensitivity"))) {
    type <- "volcano"
  }

  if (
    !is.character(type) ||
      length(type) != 1L ||
      is.na(type) ||
      !type %in% c("volcano", "sensitivity")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg type} must be either {.val volcano} or {.val sensitivity}.",
        "i" = "Received {.val {type}}."
      ),
      class = "peoplenotp_error_bad_type"
    )
  }

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

  if (type == "sensitivity") {
    return(
      plot_individual_sensitivity(
        x = x,
        mid_select = mid_select,
        thresholds = thresholds,
        outcome = outcome,
        color_palette = color_palette,
        point_size = point_size,
        point_alpha = point_alpha,
        x_lab = x_lab,
        y_lab = y_lab,
        facet = facet,
        ...
      )
    )
  }

  allowed_colors <- c(
    "signed_certainty",
    "pd",
    "makowski_tier",
    "rope_pct",
    "improvement_pct",
    "deterioration_pct",
    "directed_certainty"
  )

  color_quo <- rlang::enquo(color_by)
  color_expr <- rlang::quo_get_expr(color_quo)

  if (missing(color_by)) {
    color_col <- "signed_certainty"
  } else if (is.symbol(color_expr)) {
    color_col <- rlang::as_string(color_expr)
  } else if (
    is.character(color_by) &&
      identical(color_by, allowed_colors)
  ) {
    color_col <- "signed_certainty"
  } else if (
    is.character(color_by) &&
      length(color_by) == 1L &&
      !is.na(color_by)
  ) {
    color_col <- color_by
  } else {
    color_col <- ""
  }

  if (!color_col %in% allowed_colors) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg color_by} must be one of {.or {.val {allowed_colors}}}.",
        "i" = "Received {.val {rlang::as_label(color_quo)}}."
      ),
      class = "peoplenotp_error_bad_color_by"
    )
  }

  if (!is.logical(show_ci) || length(show_ci) != 1L || is.na(show_ci)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg show_ci} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {show_ci}}."
      ),
      class = "peoplenotp_error_bad_flag"
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

  if (is.null(x_lab)) {
    x_lab <- if (isTRUE(show_ci)) {
      "Posterior Median (with Credible Interval)"
    } else {
      "Posterior Median"
    }
  }

  if (is.null(y_lab)) {
    y_lab <- "Directed Clinical Certainty"
  }

  selected_mid <- resolve_print_mids(
    available_mids = x$meta$mid,
    mid_select = mid_select,
    mid_bounds = x$meta$mid_bounds,
    single = TRUE
  )

  plot_df <- dplyr::filter(x$individual_summary, mid == selected_mid)
  step_idx <- match(selected_mid, x$meta$mid)
  if (!is.null(x$meta$mid_bounds) && !is.na(step_idx)) {
    r_low <- x$meta$mid_bounds$rope_low[[step_idx]]
    r_high <- x$meta$mid_bounds$rope_high[[step_idx]]
  } else if (all(c("rope_low", "rope_high") %in% names(plot_df))) {
    r_low <- plot_df$rope_low[[1]]
    r_high <- plot_df$rope_high[[1]]
  } else {
    r_low <- -selected_mid
    r_high <- selected_mid
  }

  if (x$meta$better_is == "lower") {
    left_label <- "\u2190 Favors Improvement"
    right_label <- "Favors Deterioration \u2192"
  } else {
    left_label <- "\u2190 Favors Deterioration"
    right_label <- "Favors Improvement \u2192"
  }

  color_labels <- c(
    signed_certainty = "Signed Certainty",
    pd = "PD",
    makowski_tier = "Evidence Tier",
    rope_pct = "% in ROPE",
    improvement_pct = "% Improvement",
    deterioration_pct = "% Deterioration",
    directed_certainty = "Directed Certainty"
  )

  color_scale <- if (color_col == "makowski_tier") {
    rcartocolor::scale_color_carto_d(palette = "Temps", direction = -1)
  } else if (color_col == "signed_certainty") {
    rcartocolor::scale_color_carto_c(
      palette = "Temps",
      direction = -1,
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 1),
      labels = scales::label_percent()
    )
  } else {
    rcartocolor::scale_color_carto_c(
      palette = "Temps",
      direction = -1,
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::label_percent()
    )
  }

  x_ref_low <- if (isTRUE(show_ci)) plot_df$ci_low else plot_df$median
  x_ref_high <- if (isTRUE(show_ci)) plot_df$ci_high else plot_df$median

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = median,
      y = directed_certainty,
      color = .data[[color_col]]
    )
  ) +
    geom_rope(x_range = c(r_low, r_high)) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dotted"
    )

  if (isTRUE(show_ci)) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = ci_low, xmax = ci_high),
        orientation = "y",
        width = 0,
        alpha = interval_alpha
      )
  }

  p <- p +
    ggplot2::geom_point(
      size = point_size,
      alpha = point_alpha,
      ...
    ) +
    ggplot2::annotate(
      "label",
      x = 0,
      y = 1.07,
      label = left_label,
      label.r = ggplot2::unit(0, "lines"),
      hjust = 1.1,
      size = 3.5
    ) +
    ggplot2::annotate(
      "label",
      x = 0,
      y = 1.07,
      label = right_label,
      label.r = ggplot2::unit(0, "lines"),
      hjust = -0.1,
      size = 3.5
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      limits = c(0, 1.12),
      breaks = seq(0, 1, by = 0.25)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    color_scale +
    ggplot2::labs(
      x = x_lab,
      y = y_lab,
      color = unname(color_labels[[color_col]])
    )

  group_cols <- x$meta$group_cols %||% character(0)
  if (isTRUE(facet) && length(group_cols) > 0L) {
    p <- p + ggplot2::facet_wrap(stats::reformulate(group_cols))
  }

  p
}

#' Resolve and normalize `mid` input for `individual_certainty()`
#'
#' @param mid User-supplied `mid` argument.
#'
#' @return A list with `mid_vals` and `mid_bounds` tibble.
#' @noRd
resolve_individual_mid <- function(mid) {
  abort_bad_mid <- function(val) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid} must be a positive numeric scalar/vector, an asymmetric range `c(low, high)` with `low < 0 < high`, or a list of such thresholds.",
        "i" = "Received {.val {val}}."
      ),
      class = "peoplenotp_error_bad_mid"
    )
  }

  if (missing(mid) || is.null(mid)) {
    abort_bad_mid(NULL)
  }

  format_range_lbl <- function(low, high) {
    paste0(
      "[",
      format(low, trim = TRUE, nsmall = 1),
      ", ",
      format(high, trim = TRUE, nsmall = 1),
      "]"
    )
  }

  format_scalar_lbl <- function(m) {
    format(m, trim = TRUE, nsmall = 1)
  }

  if (is.numeric(mid)) {
    if (length(mid) == 0L || any(is.na(mid)) || any(is.infinite(mid))) {
      abort_bad_mid(mid)
    }

    if (length(mid) == 2L && mid[[1]] < 0 && mid[[2]] > 0) {
      low_val <- as.numeric(mid[[1]])
      high_val <- as.numeric(mid[[2]])
      lbl <- format_range_lbl(low_val, high_val)
      bounds <- tibble::tibble(
        mid = lbl,
        mid_label = lbl,
        rope_low = low_val,
        rope_high = high_val,
        is_asymmetric = TRUE
      )
      return(list(mid_vals = lbl, mid_bounds = bounds))
    }

    if (all(mid > 0)) {
      mid_num <- unique(as.numeric(mid))
      lbls <- vapply(mid_num, format_scalar_lbl, character(1))
      bounds <- tibble::tibble(
        mid = mid_num,
        mid_label = lbls,
        rope_low = -mid_num,
        rope_high = mid_num,
        is_asymmetric = FALSE
      )
      return(list(mid_vals = mid_num, mid_bounds = bounds))
    }

    abort_bad_mid(mid)
  }

  if (is.list(mid) && !is.data.frame(mid) && length(mid) >= 1L) {
    parsed_rows <- lapply(mid, function(el) {
      if (
        !is.numeric(el) ||
          length(el) == 0L ||
          any(is.na(el)) ||
          any(is.infinite(el))
      ) {
        abort_bad_mid(mid)
      }
      if (length(el) == 1L && el > 0) {
        m_val <- as.numeric(el)
        tibble::tibble(
          num_val = m_val,
          mid_label = format_scalar_lbl(m_val),
          rope_low = -m_val,
          rope_high = m_val,
          is_asymmetric = FALSE
        )
      } else if (length(el) == 2L && el[[1]] < 0 && el[[2]] > 0) {
        tibble::tibble(
          num_val = NA_real_,
          mid_label = format_range_lbl(el[[1]], el[[2]]),
          rope_low = as.numeric(el[[1]]),
          rope_high = as.numeric(el[[2]]),
          is_asymmetric = TRUE
        )
      } else {
        abort_bad_mid(mid)
      }
    })

    bounds <- dplyr::bind_rows(parsed_rows)
    bounds <- dplyr::distinct(bounds, rope_low, rope_high, .keep_all = TRUE)

    if (any(bounds$is_asymmetric)) {
      mid_vals <- bounds$mid_label
      bounds$mid <- mid_vals
    } else {
      mid_vals <- bounds$num_val
      bounds$mid <- mid_vals
    }
    bounds$num_val <- NULL
    bounds <- dplyr::select(
      bounds,
      mid,
      mid_label,
      rope_low,
      rope_high,
      is_asymmetric
    )
    return(list(mid_vals = mid_vals, mid_bounds = bounds))
  }

  abort_bad_mid(mid)
}

#' Resolve MID selection for `print.pnp_individual_certainty()` and `plot.pnp_individual_certainty()`
#'
#' @param available_mids Vector of MIDs stored in `x$meta$mid`.
#' @param mid_select User-supplied `mid_select`.
#' @param mid_bounds Optional `mid_bounds` tibble from `x$meta$mid_bounds`.
#' @param single Logical; whether only a single MID selection is allowed.
#' @param default_all Logical; whether `NULL` selects all MIDs.
#'
#' @return Vector of selected MIDs matching `available_mids` type.
#' @noRd
resolve_print_mids <- function(
  available_mids,
  mid_select = NULL,
  mid_bounds = NULL,
  single = FALSE,
  default_all = FALSE
) {
  n_mids <- length(available_mids)
  if (is.null(mid_select)) {
    if (isTRUE(single)) {
      return(available_mids[[1]])
    }
    if (isTRUE(default_all) || n_mids <= 2L) {
      return(available_mids)
    }
    mid_idx <- ceiling(n_mids / 2)
    return(unique(available_mids[c(1L, mid_idx, n_mids)]))
  }

  resolved <- NULL

  if (
    (is.numeric(mid_select) || is.character(mid_select)) &&
      length(mid_select) >= 1L &&
      !any(is.na(mid_select)) &&
      all(mid_select %in% available_mids)
  ) {
    resolved <- unique(mid_select)
  } else if (
    is.numeric(mid_select) &&
      length(mid_select) == 2L &&
      !any(is.na(mid_select)) &&
      mid_select[[1]] < 0 &&
      mid_select[[2]] > 0 &&
      !is.null(mid_bounds)
  ) {
    idx <- which(
      mid_bounds$rope_low == mid_select[[1]] &
        mid_bounds$rope_high == mid_select[[2]]
    )
    if (length(idx) >= 1L) {
      resolved <- unique(available_mids[idx])
    }
  } else if (
    is.list(mid_select) &&
      length(mid_select) >= 1L &&
      !is.null(mid_bounds)
  ) {
    matched_idx <- integer(0)
    ok <- TRUE
    for (el in mid_select) {
      if (is.numeric(el) && length(el) == 2L && !any(is.na(el))) {
        idx <- which(
          mid_bounds$rope_low == el[[1]] & mid_bounds$rope_high == el[[2]]
        )
      } else if (is.numeric(el) && length(el) == 1L && !is.na(el)) {
        idx <- which(
          mid_bounds$rope_low == -el & mid_bounds$rope_high == el
        )
      } else {
        idx <- integer(0)
      }
      if (length(idx) == 0L) {
        ok <- FALSE
        break
      }
      matched_idx <- c(matched_idx, idx)
    }
    if (ok && length(matched_idx) >= 1L) {
      resolved <- unique(available_mids[matched_idx])
    }
  } else if (
    is.character(available_mids) &&
      is.numeric(mid_select) &&
      length(mid_select) >= 1L &&
      !any(is.na(mid_select)) &&
      all(mid_select == trunc(mid_select)) &&
      all(mid_select >= 1L & mid_select <= n_mids)
  ) {
    resolved <- unique(available_mids[as.integer(mid_select)])
  }

  if (is.null(resolved) || (isTRUE(single) && length(resolved) != 1L)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid_select} must {if (single) 'be a single value' else 'contain values'} present in `x$meta$mid` ({.val {available_mids}}).",
        "i" = "Received {.val {mid_select}}."
      ),
      class = "peoplenotp_error_bad_mid_select"
    )
  }

  resolved
}

#' Build and plot the clinical certainty sensitivity heatmap for `pnp_individual_certainty`
#'
#' @param x A `pnp_individual_certainty` object.
#' @param mid_select Optional subset of `x$meta$mid`.
#' @param thresholds Numeric vector of probability thresholds in `(0, 1]`.
#' @param outcome `"both"`, `"improvement"`, or `"deterioration"`.
#' @param color_palette Fill palette name.
#' @param point_size Unused for tile plot (kept for internal signature consistency).
#' @param point_alpha Tile opacity in `[0, 1]`.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param facet Logical flag for `.by` group faceting.
#' @param ... Additional arguments passed to `geom_tile()`.
#'
#' @return A `ggplot` object.
#' @noRd
plot_individual_sensitivity <- function(
  x,
  mid_select = NULL,
  thresholds = seq(0.50, 0.99, by = 0.01),
  outcome = c("both", "improvement", "deterioration"),
  color_palette = "Temps",
  point_size = 2,
  point_alpha = 1,
  x_lab = NULL,
  y_lab = NULL,
  facet = TRUE,
  ...
) {
  if (identical(outcome, c("both", "improvement", "deterioration"))) {
    outcome <- "both"
  }

  if (
    !is.character(outcome) ||
      length(outcome) != 1L ||
      is.na(outcome) ||
      !outcome %in% c("both", "improvement", "deterioration")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg outcome} must be one of {.val both}, {.val improvement}, or {.val deterioration}.",
        "i" = "Received {.val {outcome}}."
      ),
      class = "peoplenotp_error_bad_outcome"
    )
  }

  if (
    !is.numeric(thresholds) ||
      length(thresholds) == 0L ||
      any(is.na(thresholds)) ||
      any(thresholds <= 0) ||
      any(thresholds > 1)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg thresholds} must be a numeric vector with values strictly between 0 and 1 (`(0, 1]`).",
        "i" = "Received {.val {thresholds}}."
      ),
      class = "peoplenotp_error_bad_thresholds"
    )
  }

  selected_mids <- resolve_print_mids(
    available_mids = x$meta$mid,
    mid_select = mid_select,
    mid_bounds = x$meta$mid_bounds,
    single = FALSE,
    default_all = TRUE
  )

  fill_scale <- resolve_sensitivity_fill_scale(
    color_palette = color_palette,
    outcome = outcome
  )

  if (is.null(x_lab)) {
    x_lab <- "Minimal Important Difference (MID)"
  }
  if (is.null(y_lab)) {
    y_lab <- "Certainty Threshold"
  }

  group_cols <- x$meta$group_cols %||% character(0)
  thresh_vals <- sort(unique(thresholds))
  t_base <- min(thresh_vals)
  t_step <- if (length(thresh_vals) > 1L) min(diff(thresh_vals)) else 0.01
  t_span <- (max(thresh_vals) - t_base) + t_step

  grid_df <- build_sensitivity_grid(
    individual_summary = x$individual_summary,
    selected_mids = selected_mids,
    thresholds = thresh_vals,
    group_cols = group_cols,
    outcome = outcome,
    t_base = t_base,
    t_step = t_step
  )
  if (is.character(selected_mids)) {
    grid_df$mid <- factor(grid_df$mid, levels = selected_mids)
  }

  if (outcome == "both") {
    x_annot <- if (is.character(selected_mids)) {
      (length(selected_mids) + 1) / 2
    } else {
      mean(range(selected_mids))
    }
    y_annot <- t_span * 1.12
    y_lim <- t_span * 1.24
    y_breaks <- seq(-t_span, t_span, length.out = 5)

    p <- ggplot2::ggplot(
      grid_df,
      ggplot2::aes(
        x = mid,
        y = y_tile,
        fill = prop_signed
      )
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(height = tile_height),
        alpha = point_alpha,
        ...
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted"
      ) +
      ggplot2::annotate(
        "label",
        x = x_annot,
        y = y_annot,
        label = "\u2191 Improvement",
        label.r = ggplot2::unit(0, "lines"),
        size = 3.5
      ) +
      ggplot2::annotate(
        "label",
        x = x_annot,
        y = -y_annot,
        label = "\u2193 Deterioration",
        label.r = ggplot2::unit(0, "lines"),
        size = 3.5
      ) +
      ggplot2::scale_y_continuous(
        labels = function(b) scales::label_percent()(t_base + abs(b)),
        limits = c(-y_lim, y_lim),
        breaks = y_breaks
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      fill_scale +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        fill = "Proportion Certain"
      )
  } else {
    p <- ggplot2::ggplot(
      grid_df,
      ggplot2::aes(
        x = mid,
        y = y_tile,
        fill = prop_certain
      )
    ) +
      ggplot2::geom_tile(
        ggplot2::aes(height = tile_height),
        alpha = point_alpha,
        ...
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(),
        limits = c(t_base - t_step / 2, max(thresh_vals) + t_step / 2)
      ) +
      fill_scale +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        fill = "Proportion Certain"
      )
  }

  if (isTRUE(facet) && length(group_cols) > 0L) {
    p <- p + ggplot2::facet_wrap(stats::reformulate(group_cols))
  }

  p
}

#' Build summarized sensitivity grid from `$individual_summary`
#'
#' @param individual_summary Tibble from `x$individual_summary`.
#' @param selected_mids Numeric vector of MID values to include.
#' @param thresholds Sorted numeric vector of probability thresholds.
#' @param group_cols Character vector of grouping column names.
#' @param outcome `"both"`, `"improvement"`, or `"deterioration"`.
#' @param t_base Minimum threshold value used as mirror baseline.
#' @param t_step Step height of each threshold tile.
#'
#' @return A summarized tibble for plotting with `geom_tile()`.
#' @noRd
build_sensitivity_grid <- function(
  individual_summary,
  selected_mids,
  thresholds,
  group_cols = character(0),
  outcome = "both",
  t_base = min(thresholds),
  t_step = 0.01
) {
  sub_df <- dplyr::filter(individual_summary, mid %in% selected_mids)
  thresh_vals <- sort(unique(thresholds))

  summary_by <- c(group_cols, "mid")

  thresh_list <- lapply(thresh_vals, function(t_val) {
    agg <- dplyr::summarise(
      sub_df,
      n_total = dplyr::n(),
      prop_imp = sum(improvement_pct >= t_val, na.rm = TRUE) / n_total,
      prop_det = sum(deterioration_pct >= t_val, na.rm = TRUE) / n_total,
      .by = dplyr::all_of(summary_by)
    )
    agg$threshold <- t_val
    agg$tile_height <- t_step
    agg
  })

  combined <- dplyr::bind_rows(thresh_list)

  if (outcome == "improvement") {
    return(
      dplyr::mutate(
        combined,
        direction = "Improvement",
        y_tile = threshold,
        prop_certain = prop_imp,
        prop_signed = prop_imp
      )
    )
  }

  if (outcome == "deterioration") {
    return(
      dplyr::mutate(
        combined,
        direction = "Deterioration",
        y_tile = threshold,
        prop_certain = prop_det,
        prop_signed = prop_det
      )
    )
  }

  imp_rows <- dplyr::mutate(
    combined,
    direction = "Improvement",
    y_tile = (threshold - t_base) + t_step / 2,
    prop_certain = prop_imp,
    prop_signed = prop_imp
  )
  det_rows <- dplyr::mutate(
    combined,
    direction = "Deterioration",
    y_tile = -(threshold - t_base) - t_step / 2,
    prop_certain = prop_det,
    prop_signed = -prop_det
  )

  dplyr::bind_rows(imp_rows, det_rows)
}

#' Resolve continuous fill scale for sensitivity heatmap tiles
#'
#' @param color_palette Palette name.
#' @param outcome `"both"`, `"improvement"`, or `"deterioration"`.
#'
#' @return A `ggplot2` continuous fill scale.
#' @noRd
resolve_sensitivity_fill_scale <- function(color_palette, outcome = "both") {
  carto_options <- rcartocolor::cartocolors$Name
  viridis_options <- c(
    "viridis",
    "magma",
    "inferno",
    "plasma",
    "cividis",
    "rocket",
    "mako",
    "turbo",
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H"
  )
  custom_options <- names(peoplenotp::custom_palettes)

  if (
    !is.character(color_palette) ||
      length(color_palette) != 1L ||
      is.na(color_palette) ||
      !color_palette %in% c(carto_options, viridis_options, custom_options)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg color_palette} must be a valid {.pkg rcartocolor} palette, {.pkg viridisLite} option, or a palette in {.code custom_palettes}.",
        "i" = "Received {.val {color_palette}}."
      ),
      class = "peoplenotp_error_bad_color_palette"
    )
  }

  scale_limits <- if (outcome == "both") c(-1, 1) else c(0, 1)
  scale_breaks <- if (outcome == "both") {
    seq(-1, 1, by = 1)
  } else {
    seq(0, 1, by = 0.25)
  }
  scale_labels <- function(b) scales::label_percent()(abs(b))

  if (color_palette %in% carto_options) {
    rcartocolor::scale_fill_carto_c(
      palette = color_palette,
      direction = -1,
      limits = scale_limits,
      breaks = scale_breaks,
      labels = scale_labels
    )
  } else if (color_palette %in% viridis_options) {
    ggplot2::scale_fill_viridis_c(
      option = color_palette,
      limits = scale_limits,
      breaks = scale_breaks,
      labels = scale_labels
    )
  } else {
    ggplot2::scale_fill_gradientn(
      colours = unname(peoplenotp::custom_palettes[[color_palette]]),
      limits = scale_limits,
      breaks = scale_breaks,
      labels = scale_labels
    )
  }
}

#' Ordered directed Makowski tier factor levels (6-Tier symmetric)
#'
#' @return Character vector of length 6.
#' @noRd
makowski_tier_levels <- function() {
  c(
    "Definite Deterioration (>99%)",
    "Probable Deterioration (>97.5%)",
    "Undecided",
    "Practically Equivalent (>97.5% ROPE)",
    "Probable Improvement (>97.5%)",
    "Definite Improvement (>99%)"
  )
}

#' Classify directional and ROPE probabilities into ordered Makowski evidence tiers
#'
#' @param rope_pct Numeric vector of ROPE probabilities in `[0, 1]`.
#' @param imp_pct Numeric vector of improvement probabilities in `[0, 1]`.
#' @param det_pct Numeric vector of deterioration probabilities in `[0, 1]`.
#'
#' @return Ordered factor.
#' @noRd
assign_makowski_tier <- function(rope_pct, imp_pct, det_pct) {
  tiers <- dplyr::case_when(
    det_pct > 0.99 ~ "Definite Deterioration (>99%)",
    det_pct > 0.975 ~ "Probable Deterioration (>97.5%)",
    rope_pct > 0.975 ~ "Practically Equivalent (>97.5% ROPE)",
    imp_pct > 0.99 ~ "Definite Improvement (>99%)",
    imp_pct > 0.975 ~ "Probable Improvement (>97.5%)",
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
resolve_id_column <- function(data, id_quo, id_missing = FALSE) {
  if (rlang::quo_is_null(id_quo)) {
    return(NULL)
  }
  if (id_missing && !"id" %in% names(data)) {
    return(NULL)
  }
  resolve_column_name(data, id_quo)
}

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

#' Resolve `.by` grouping columns or existing `dplyr::group_by()` variables
#'
#' @param data Original input data frame.
#' @param ungrouped_data Ungrouped input data frame.
#' @param by_quo Quosure for `.by`.
#' @param id_name Resolved individual ID column name.
#'
#' @return Character vector of additional grouping column names.
#' @noRd
resolve_group_columns <- function(data, ungrouped_data, by_quo, id_name) {
  if (!rlang::quo_is_null(by_quo)) {
    selected <- tryCatch(
      names(dplyr::select(ungrouped_data[0, , drop = FALSE], !!by_quo)),
      error = function(cnd) {
        cli::cli_abort(
          c(
            "x" = "Grouping column(s) {.code {rlang::as_label(by_quo)}} not found in {.arg data}.",
            "i" = "Available columns: {.val {names(ungrouped_data)}}."
          ),
          class = "peoplenotp_error_missing_column",
          parent = cnd
        )
      }
    )
    return(setdiff(selected, c(id_name, ".row")))
  }

  if (dplyr::is_grouped_df(data)) {
    return(setdiff(dplyr::group_vars(data), c(id_name, ".row")))
  }

  character(0)
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
#' @param m MID step value or label.
#' @param rope_low Lower ROPE boundary.
#' @param rope_high Upper ROPE boundary.
#' @param better_is `"lower"` or `"higher"`.
#'
#' @return A single-row tibble.
#' @noRd
compute_individual_mid <- function(
  draws,
  m,
  rope_low = -m,
  rope_high = m,
  better_is
) {
  rope_pct <- mean(draws >= rope_low & draws <= rope_high, na.rm = TRUE)

  if (better_is == "lower") {
    imp_pct <- mean(draws < rope_low, na.rm = TRUE)
    det_pct <- mean(draws > rope_high, na.rm = TRUE)
  } else {
    imp_pct <- mean(draws > rope_high, na.rm = TRUE)
    det_pct <- mean(draws < rope_low, na.rm = TRUE)
  }

  dir_cert <- pmax(imp_pct, det_pct)
  signed_cert <- imp_pct - det_pct

  tibble::tibble(
    mid = m,
    rope_low = rope_low,
    rope_high = rope_high,
    rope_pct = rope_pct,
    improvement_pct = imp_pct,
    deterioration_pct = det_pct,
    directed_certainty = dir_cert,
    signed_certainty = signed_cert,
    makowski_tier = assign_makowski_tier(rope_pct, imp_pct, det_pct)
  )
}

#' Build long evidence matrix counts and proportions (optionally stratified)
#'
#' @param individual_summary Individual summary tibble.
#' @param mid_vals Numeric vector of MID values.
#' @param group_cols Character vector of additional grouping column names.
#'
#' @return A tibble with `group_cols` (if any), `makowski_tier`, `mid`, `n`,
#'   `prop`, `pct`.
#' @noRd
build_evidence_matrix_long <- function(
  individual_summary,
  mid_vals,
  group_cols = character(0)
) {
  all_tiers <- factor(
    makowski_tier_levels(),
    levels = makowski_tier_levels(),
    ordered = TRUE
  )

  if (length(group_cols) == 0L) {
    n_individuals <- nrow(dplyr::filter(
      individual_summary,
      mid == mid_vals[[1]]
    ))
    grid <- tidyr::expand_grid(
      makowski_tier = all_tiers,
      mid = mid_vals
    )

    counts <- dplyr::summarise(
      individual_summary,
      n = dplyr::n(),
      .by = c(makowski_tier, mid)
    )

    return(
      dplyr::left_join(grid, counts, by = c("makowski_tier", "mid")) |>
        dplyr::mutate(
          n = dplyr::coalesce(n, 0L),
          prop = if (n_individuals > 0L) n / n_individuals else 0,
          pct = prop * 100
        )
    )
  }

  first_mid_df <- dplyr::filter(individual_summary, mid == mid_vals[[1]])
  grp_totals <- dplyr::summarise(
    first_mid_df,
    n_group = dplyr::n(),
    .by = dplyr::all_of(group_cols)
  )

  grid <- tidyr::expand_grid(
    grp_totals,
    makowski_tier = all_tiers,
    mid = mid_vals
  )

  counts <- dplyr::summarise(
    individual_summary,
    n = dplyr::n(),
    .by = c(dplyr::all_of(group_cols), makowski_tier, mid)
  )

  dplyr::left_join(
    grid,
    counts,
    by = c(group_cols, "makowski_tier", "mid")
  ) |>
    dplyr::mutate(
      n = dplyr::coalesce(n, 0L),
      prop = ifelse(n_group > 0L, n / n_group, 0),
      pct = prop * 100
    ) |>
    dplyr::select(-n_group)
}

#' Pivot long evidence matrix into wide inverted matrix format
#'
#' @param matrix_long Long evidence matrix tibble.
#' @param mid_vals Numeric vector of MID values.
#' @param group_cols Character vector of additional grouping column names.
#' @param format `"both"`, `"count"`, or `"prop"`.
#' @param digits Number of decimal digits for percentages.
#'
#' @return A wide tibble with `group_cols` (if any), `makowski_tier` rows, and
#'   `MID = ...` columns.
#' @noRd
format_evidence_matrix_wide <- function(
  matrix_long,
  mid_vals,
  group_cols = character(0),
  format = "both",
  digits = 1
) {
  formatted_mids <- if (is.numeric(mid_vals)) {
    format(mid_vals, trim = TRUE, nsmall = 1)
  } else {
    as.character(mid_vals)
  }
  mid_labels <- stats::setNames(
    paste0("MID = ", formatted_mids),
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
    dplyr::select(
      formatted_long,
      dplyr::all_of(group_cols),
      makowski_tier,
      mid_col,
      cell
    ),
    names_from = mid_col,
    values_from = cell
  )
}
