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
#'   uniquely identifying individuals (default `id`).
#' @param value <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name containing the individual target estimate per draw (e.g., `.epred`,
#'   `net_change`). Defaults to `.epred`.
#' @param mid Numeric scalar (e.g., `2.5`) or numeric vector of positive values
#'   (e.g., `seq(1.5, 3.5, by = 0.5)`) defining the Minimal Important
#'   Difference threshold(s).
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
#'     and `mid` value containing `<id>`, any `.by` grouping columns, `mid`,
#'     `median`, `ci_low`, `ci_high`, `pd`, `rope_pct`, `improvement_pct`,
#'     `deterioration_pct`, `directed_certainty`, `signed_certainty`, and
#'     `makowski_tier` (an ordered factor).
#'   * `$matrix_data`: A [tibble::tibble()] representing the inverted Makowski
#'     Evidence Matrix crossing `makowski_tier` (rows) by `mid` thresholds
#'     (columns), stratified by `.by` if supplied.
#'   * `$meta`: A list of metadata (`mid`, `better_is`, `ci`, `id_col`,
#'     `group_cols`, `value_col`, `draw_col`, `n_individuals`, and
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
  by_quo <- rlang::enquo(.by)
  draw_quo <- rlang::enquo(draw_col)

  id_name <- resolve_column_name(ungrouped_data, id_quo)
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

  mid_vals <- unique(mid)
  summary_by <- c(id_name, group_cols)

  base_summary <- dplyr::summarise(
    ungrouped_data,
    compute_individual_base(
      draws = .data[[value_name]],
      ci = ci
    ),
    .by = dplyr::all_of(summary_by)
  )

  mid_summaries <- lapply(mid_vals, function(m) {
    rope_df <- dplyr::summarise(
      ungrouped_data,
      compute_individual_mid(
        draws = .data[[value_name]],
        m = m,
        better_is = better_is
      ),
      .by = dplyr::all_of(summary_by)
    )

    merged <- dplyr::left_join(base_summary, rope_df, by = summary_by)
    dplyr::select(
      merged,
      dplyr::all_of(summary_by),
      mid,
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
    mid_select = mid_select
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
#' * **`type = "sensitivity"`:** A clinical certainty sensitivity curve across
#'   `mid` values on the x-axis and discrete probability `thresholds` as colored
#'   lines. By default (`outcome = "both"`), renders a mirrored butterfly plot
#'   showing the proportion of individuals achieving certainty for improvement
#'   upward (`0%` to `+100%`) and for deterioration downward (`0%` to `100%`),
#'   or single-sided curves when `outcome = "improvement"` or
#'   `"deterioration"`.
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
#' @param thresholds Numeric vector of probability thresholds in `(0, 1]` used
#'   when `type = "sensitivity"`. Defaults to `c(0.5, 0.7, 0.8, 0.9, 0.95)`.
#' @param outcome Character string specifying which clinical direction(s) to
#'   display when `type = "sensitivity"`: `"both"` (default, mirrored butterfly
#'   plot with improvement above 0 and deterioration below 0), `"improvement"`,
#'   or `"deterioration"`.
#' @param color_palette Character string specifying the discrete color palette
#'   for `thresholds` when `type = "sensitivity"`. Defaults to `"Temps"` (via
#'   `{rcartocolor}`). Also supports any `{rcartocolor}` palette, any
#'   `{viridisLite}` option (`"viridis"`, `"magma"`, `"inferno"`, `"plasma"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`), or any palette name in
#'   [custom_palettes].
#' @param show_ci Logical flag (`TRUE` by default). For `type = "volcano"`,
#'   draws horizontal credible interval error bars (`ci_low` to `ci_high`) when
#'   `TRUE`.
#' @param point_size Numeric controlling point size (default `2`).
#' @param point_alpha Numeric in `[0, 1]` controlling point opacity (default
#'   `1`).
#' @param interval_alpha Numeric in `[0, 1]` controlling horizontal error bar
#'   opacity in the volcano plot (default `0.5`).
#' @param x_lab Optional string controlling the x-axis label. When `NULL`
#'   (default), uses `"Posterior Median (with Credible Interval)"` (or
#'   `"Posterior Median"` when `show_ci = FALSE`) for `type = "volcano"`, and
#'   `"Minimal Important Difference (MID)"` for `type = "sensitivity"`.
#' @param y_lab Optional string controlling the y-axis label. When `NULL`
#'   (default), uses `"Directed Clinical Certainty"` for `type = "volcano"`, and
#'   `"Proportion Meeting Certainty"` for `type = "sensitivity"`.
#' @param facet Logical flag (`TRUE` by default). When `TRUE` and `x` contains
#'   `.by` grouping columns in `x$meta$group_cols`, facets the plot via
#'   [ggplot2::facet_wrap()].
#' @param ... Additional arguments passed to [ggplot2::geom_point()] (and
#'   [ggplot2::geom_line()] when `type = "sensitivity"`).
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
  thresholds = c(0.5, 0.7, 0.8, 0.9, 0.95),
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
    geom_rope(x_range = c(-selected_mid, selected_mid)) +
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

#' Resolve MID selection for `print.pnp_individual_certainty()`
#'
#' @param available_mids Numeric vector of MIDs stored in `x$meta$mid`.
#' @param mid_select User-supplied `mid_select` (`NULL` or numeric vector).
#'
#' @return Numeric vector of selected MIDs.
#' @noRd
resolve_print_mids <- function(available_mids, mid_select = NULL) {
  if (is.null(mid_select)) {
    n_mids <- length(available_mids)
    if (n_mids <= 2L) {
      return(available_mids)
    }
    mid_idx <- ceiling(n_mids / 2)
    return(unique(available_mids[c(1L, mid_idx, n_mids)]))
  }

  if (
    !is.numeric(mid_select) ||
      length(mid_select) == 0L ||
      any(is.na(mid_select)) ||
      !all(mid_select %in% available_mids)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid_select} must contain values present in `x$meta$mid` ({.val {available_mids}}).",
        "i" = "Received {.val {mid_select}}."
      ),
      class = "peoplenotp_error_bad_mid_select"
    )
  }

  unique(mid_select)
}

#' Build and plot the clinical certainty sensitivity curve for `pnp_individual_certainty`
#'
#' @param x A `pnp_individual_certainty` object.
#' @param mid_select Optional subset of `x$meta$mid`.
#' @param thresholds Numeric vector of probability thresholds in `(0, 1]`.
#' @param outcome `"both"`, `"improvement"`, or `"deterioration"`.
#' @param color_palette Discrete palette name.
#' @param point_size Point size.
#' @param point_alpha Point opacity.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param facet Logical flag for `.by` group faceting.
#' @param ... Additional arguments passed to `geom_line()` and `geom_point()`.
#'
#' @return A `ggplot` object.
#' @noRd
plot_individual_sensitivity <- function(
  x,
  mid_select = NULL,
  thresholds = c(0.5, 0.7, 0.8, 0.9, 0.95),
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

  available_mids <- x$meta$mid
  if (is.null(mid_select)) {
    selected_mids <- available_mids
  } else {
    if (
      !is.numeric(mid_select) ||
        length(mid_select) == 0L ||
        any(is.na(mid_select)) ||
        !all(mid_select %in% available_mids)
    ) {
      cli::cli_abort(
        c(
          "x" = "Argument {.arg mid_select} must contain values present in `x$meta$mid` ({.val {available_mids}}).",
          "i" = "Received {.val {mid_select}}."
        ),
        class = "peoplenotp_error_bad_mid_select"
      )
    }
    selected_mids <- unique(mid_select)
  }

  color_scale <- resolve_sensitivity_color_scale(color_palette)

  if (is.null(x_lab)) {
    x_lab <- "Minimal Important Difference (MID)"
  }
  if (is.null(y_lab)) {
    y_lab <- "Proportion Meeting Certainty"
  }

  group_cols <- x$meta$group_cols %||% character(0)
  grid_df <- build_sensitivity_grid(
    individual_summary = x$individual_summary,
    selected_mids = selected_mids,
    thresholds = thresholds,
    group_cols = group_cols,
    outcome = outcome
  )

  if (outcome == "both") {
    x_annot <- mean(range(selected_mids))

    p <- ggplot2::ggplot(
      grid_df,
      ggplot2::aes(
        x = mid,
        y = prop_signed,
        color = threshold_label,
        group = interaction(threshold_label, direction)
      )
    ) +
      ggplot2::geom_hline(
        yintercept = 0,
        linetype = "dotted"
      ) +
      ggplot2::geom_line(...) +
      ggplot2::geom_point(
        size = point_size,
        alpha = point_alpha,
        ...
      ) +
      ggplot2::annotate(
        "label",
        x = x_annot,
        y = 1.07,
        label = "\u2191 Improvement",
        label.r = ggplot2::unit(0, "lines"),
        size = 3.5
      ) +
      ggplot2::annotate(
        "label",
        x = x_annot,
        y = -1.07,
        label = "\u2193 Deterioration",
        label.r = ggplot2::unit(0, "lines"),
        size = 3.5
      ) +
      ggplot2::scale_y_continuous(
        labels = function(b) scales::label_percent()(abs(b)),
        limits = c(-1.14, 1.14),
        breaks = seq(-1, 1, by = 0.25)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      color_scale +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        color = "Certainty Threshold"
      )
  } else {
    p <- ggplot2::ggplot(
      grid_df,
      ggplot2::aes(
        x = mid,
        y = prop_certain,
        color = threshold_label,
        group = threshold_label
      )
    ) +
      ggplot2::geom_line(...) +
      ggplot2::geom_point(
        size = point_size,
        alpha = point_alpha,
        ...
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25)
      ) +
      color_scale +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        color = "Certainty Threshold"
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
#' @param thresholds Numeric vector of probability thresholds.
#' @param group_cols Character vector of grouping column names.
#' @param outcome `"both"`, `"improvement"`, or `"deterioration"`.
#'
#' @return A summarized tibble for plotting.
#' @noRd
build_sensitivity_grid <- function(
  individual_summary,
  selected_mids,
  thresholds,
  group_cols = character(0),
  outcome = "both"
) {
  sub_df <- dplyr::filter(individual_summary, mid %in% selected_mids)
  thresh_vals <- sort(unique(thresholds))
  thresh_labels <- scales::label_percent()(thresh_vals)
  label_map <- stats::setNames(thresh_labels, as.character(thresh_vals))

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
    agg$threshold_label <- factor(
      unname(label_map[as.character(t_val)]),
      levels = thresh_labels,
      ordered = TRUE
    )
    agg
  })

  combined <- dplyr::bind_rows(thresh_list)

  if (outcome == "improvement") {
    return(
      dplyr::mutate(
        combined,
        direction = "Improvement",
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
        prop_certain = prop_det,
        prop_signed = prop_det
      )
    )
  }

  imp_rows <- dplyr::mutate(
    combined,
    direction = "Improvement",
    prop_certain = prop_imp,
    prop_signed = prop_imp
  )
  det_rows <- dplyr::mutate(
    combined,
    direction = "Deterioration",
    prop_certain = prop_det,
    prop_signed = -prop_det
  )

  dplyr::bind_rows(imp_rows, det_rows)
}

#' Resolve discrete color scale for sensitivity threshold curves
#'
#' @param color_palette Palette name.
#'
#' @return A `ggplot2` discrete color scale.
#' @noRd
resolve_sensitivity_color_scale <- function(color_palette) {
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

  if (color_palette %in% carto_options) {
    rcartocolor::scale_color_carto_d(palette = color_palette, direction = -1)
  } else if (color_palette %in% viridis_options) {
    ggplot2::scale_color_viridis_d(option = color_palette)
  } else {
    ggplot2::scale_color_manual(
      values = unname(peoplenotp::custom_palettes[[color_palette]])
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
  signed_cert <- imp_pct - det_pct

  tibble::tibble(
    mid = m,
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
