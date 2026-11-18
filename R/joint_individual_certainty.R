#' Evaluate joint individual-level Bayesian clinical certainty across multiple outcomes
#'
#' @description
#' Summarizes individual-level posterior draws simultaneously across two or more
#' co-primary outcome variables (`values`) and one or more Minimal Important
#' Difference (`mid`) thresholds. By evaluating directional and ROPE conditions
#' row-wise per MCMC draw, it naturally accounts for the posterior correlation
#' among outcomes within each individual.
#'
#' Computes outcome-specific posterior medians, credible intervals, and
#' marginal Probabilities of Direction (`pd_<outcome>`), alongside joint
#' posterior probabilities (`joint_pd`, `joint_rope_pct`,
#' `joint_improvement_pct`, `joint_deterioration_pct`), `directed_certainty`,
#' `signed_certainty`, ordered Makowski evidence tiers (`makowski_tier`), and an
#' inverted Joint Clinical Significance Evidence Matrix (optionally stratified
#' by `.by` grouping variables).
#'
#' @param data A tidy data frame or tibble containing posterior draws per
#'   individual.
#' @param id <[`data-masking`][rlang::args_data_masking]> Unquoted column name
#'   uniquely identifying individuals (default `id`), or `NULL` for single-case
#'   ($N = 1$) datasets without an ID column.
#' @param values <[`tidy-select`][dplyr::dplyr_tidy_select]> Two or more numeric
#'   columns containing the individual posterior draws for each outcome (e.g.,
#'   `values = c(strength, anxiety)`).
#' @param mid Minimal Important Difference threshold(s) for the selected
#'   `values`. Can be supplied as:
#'   * A positive numeric scalar (e.g., `2.5`), recycled across all $K$ outcomes
#'     for a single joint MID evaluation.
#'   * A positive numeric vector of length $K$ (e.g., `c(5, 2)` or named
#'     `c(strength = 5, anxiety = 2)`), defining one outcome-specific joint MID
#'     tuple.
#'   * A list of length $K$ of positive numeric vectors (e.g.,
#'     `list(strength = c(4, 5, 6), anxiety = c(1.5, 2.0, 2.5))`), defining $M$
#'     paired joint MID steps across outcomes for sensitivity analyses and
#'     multi-column Evidence Matrices. Elements of length 1 are recycled to $M$.
#' @param .by <[`tidy-select`][dplyr::dplyr_tidy_select]> Optional additional
#'   grouping column(s) such as therapy condition (e.g., `.by = group`). When
#'   provided (or when `data` is grouped via [dplyr::group_by()]),
#'   `$individual_summary`, `$matrix_data`, `print()`, and `plot()` are
#'   stratified across these groups. Defaults to `NULL`.
#' @param better_is Character scalar (`"lower"` or `"higher"`, recycled across
#'   all outcomes) or character vector of length $K$ (e.g.,
#'   `c("higher", "lower")`, optionally named by outcome) indicating clinical
#'   polarity for each outcome. Defaults to `"lower"`.
#' @param draw_col <[`data-masking`][rlang::args_data_masking]> Unquoted column
#'   name identifying the MCMC draw index (default `.draw`).
#' @param ci Credible interval width for individual outcome summaries (default
#'   `0.95`).
#' @param include_marginals Logical flag (`FALSE` by default). When `TRUE`,
#'   also includes outcome-specific marginal ROPE probabilities
#'   (`rope_pct_<outcome>`, `improvement_pct_<outcome>`,
#'   `deterioration_pct_<outcome>`) in `$individual_summary`.
#'
#' @return An S3 list object of class `"pnp_joint_individual_certainty"`
#'   containing:
#'   * `$individual_summary`: A [tibble::tibble()] with one row per individual
#'     and joint `mid` step containing `<id>`, any `.by` grouping columns,
#'     `mid` (character label of the joint MID tuple), `mid_<outcome>` columns,
#'     `median_<outcome>`, `ci_low_<outcome>`, `ci_high_<outcome>`, and
#'     `pd_<outcome>` for each outcome (plus marginal ROPE columns when
#'     `include_marginals = TRUE`), followed by `joint_pd`, `joint_rope_pct`,
#'     `joint_improvement_pct`, `joint_deterioration_pct`, `directed_certainty`,
#'     `signed_certainty`, and `makowski_tier` (an ordered factor).
#'   * `$matrix_data`: A [tibble::tibble()] representing the inverted Makowski
#'     Evidence Matrix crossing `makowski_tier` (rows) by joint `mid` steps
#'     (columns), stratified by `.by` if supplied.
#'   * `$meta`: A list of metadata (`mid`, `mid_grid`, `better_is`, `ci`,
#'     `id_col`, `group_cols`, `value_cols`, `draw_col`, `include_marginals`,
#'     `n_individuals`, and `matrix_long`).
#'
#' @examples
#' set.seed(123)
#' draws_df <- tibble::tibble(
#'   id = rep(paste0("P", 1:10), each = 200),
#'   group = rep(rep(c("Waitlist", "Intervention"), each = 5), each = 200),
#'   .draw = rep(1:200, times = 10),
#'   strength = stats::rnorm(
#'     2000,
#'     mean = rep(seq(-2, 8, length.out = 10), each = 200),
#'     sd = 1.5
#'   ),
#'   anxiety = stats::rnorm(
#'     2000,
#'     mean = rep(seq(2, -5, length.out = 10), each = 200),
#'     sd = 1.0
#'   )
#' )
#'
#' res <- joint_individual_certainty(
#'   data = draws_df,
#'   id = id,
#'   values = c(strength, anxiety),
#'   mid = list(strength = c(4, 5), anxiety = c(1.5, 2.0)),
#'   .by = group,
#'   better_is = c("higher", "lower")
#' )
#'
#' print(res)
#' plot(res)
#' plot(res, type = "sensitivity", outcome = "improvement")
#'
#' @export
joint_individual_certainty <- function(
  data,
  id = id,
  values,
  mid,
  .by = NULL,
  better_is = "lower",
  draw_col = .draw,
  ci = 0.95,
  include_marginals = FALSE
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
      class = c(
        "peoplenotp_error_bad_include_marginals",
        "peoplenotp_error_bad_flag"
      )
    )
  }

  ungrouped_data <- dplyr::ungroup(data)

  id_quo <- rlang::enquo(id)
  values_quo <- rlang::enquo(values)
  by_quo <- rlang::enquo(.by)
  draw_quo <- rlang::enquo(draw_col)

  if (missing(values) || rlang::quo_is_null(values_quo)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg values} must select at least two numeric outcome columns.",
        "i" = "Example: {.code values = c(outcome_1, outcome_2)}."
      ),
      class = "peoplenotp_error_bad_values"
    )
  }

  id_name <- resolve_id_column(
    data = ungrouped_data,
    id_quo = id_quo,
    id_missing = missing(id)
  )

  value_names <- resolve_value_columns(ungrouped_data, values_quo)
  n_outcomes <- length(value_names)

  for (v_name in value_names) {
    val_vec <- ungrouped_data[[v_name]]
    if (!is.numeric(val_vec)) {
      cli::cli_abort(
        c(
          "x" = "Column {.code {v_name}} must be numeric.",
          "i" = "Column has class {.cls {class(val_vec)}}."
        ),
        class = "peoplenotp_error_bad_column_type"
      )
    }
  }

  better_is_vec <- resolve_joint_better_is(better_is, value_names)
  mid_info <- resolve_joint_mid(mid, value_names)
  mid_grid <- mid_info$mid_grid
  mid_labels <- mid_info$mid_labels

  group_cols <- resolve_group_columns(
    data = data,
    ungrouped_data = ungrouped_data,
    by_quo = by_quo,
    id_name = id_name
  )
  group_cols <- setdiff(group_cols, value_names)

  if (!missing(draw_col)) {
    draw_name <- resolve_column_name(ungrouped_data, draw_quo)
  } else {
    draw_label <- rlang::as_label(draw_quo)
    draw_name <- if (draw_label %in% names(ungrouped_data)) draw_label else NULL
  }

  summary_by <- c(id_name, group_cols)

  base_summary <- dplyr::summarise(
    ungrouped_data,
    compute_joint_base(
      df_sub = dplyr::pick(dplyr::all_of(value_names)),
      value_names = value_names,
      ci = ci
    ),
    .by = dplyr::all_of(summary_by)
  )

  mid_summaries <- lapply(seq_len(nrow(mid_grid)), function(step_idx) {
    step_mids <- as.numeric(mid_grid[step_idx, ])
    names(step_mids) <- value_names
    step_label <- mid_labels[[step_idx]]

    rope_df <- dplyr::summarise(
      ungrouped_data,
      compute_joint_mid_step(
        df_sub = dplyr::pick(dplyr::all_of(value_names)),
        value_names = value_names,
        step_mids = step_mids,
        step_label = step_label,
        better_is_vec = better_is_vec,
        include_marginals = include_marginals
      ),
      .by = dplyr::all_of(summary_by)
    )

    merged <- if (length(summary_by) == 0L) {
      dplyr::bind_cols(base_summary, rope_df)
    } else {
      dplyr::left_join(base_summary, rope_df, by = summary_by)
    }

    mid_col_names <- paste0("mid_", value_names)
    base_col_names <- unlist(
      lapply(value_names, function(vn) {
        cols <- paste0(c("median_", "ci_low_", "ci_high_", "pd_"), vn)
        if (isTRUE(include_marginals)) {
          cols <- c(
            cols,
            paste0(c("rope_pct_", "improvement_pct_", "deterioration_pct_"), vn)
          )
        }
        cols
      }),
      use.names = FALSE
    )

    dplyr::select(
      merged,
      dplyr::all_of(summary_by),
      mid,
      dplyr::all_of(mid_col_names),
      dplyr::all_of(base_col_names),
      joint_pd,
      joint_rope_pct,
      joint_improvement_pct,
      joint_deterioration_pct,
      directed_certainty,
      signed_certainty,
      makowski_tier
    )
  })

  individual_summary <- tibble::as_tibble(dplyr::bind_rows(mid_summaries))
  n_individuals <- nrow(base_summary)

  matrix_long <- build_evidence_matrix_long(
    individual_summary = individual_summary,
    mid_vals = mid_labels,
    group_cols = group_cols
  )

  matrix_data <- format_joint_evidence_matrix_wide(
    matrix_long = matrix_long,
    mid_labels = mid_labels,
    group_cols = group_cols,
    format = "both",
    digits = 1
  )
  attr(matrix_data, "long") <- matrix_long

  meta <- list(
    mid = mid_labels,
    mid_grid = mid_grid,
    better_is = better_is_vec,
    ci = ci,
    id_col = id_name,
    group_cols = group_cols,
    value_cols = value_names,
    draw_col = draw_name,
    include_marginals = include_marginals,
    n_individuals = n_individuals,
    matrix_long = matrix_long
  )

  structure(
    list(
      individual_summary = individual_summary,
      matrix_data = matrix_data,
      meta = meta
    ),
    class = "pnp_joint_individual_certainty"
  )
}

#' Print method for `pnp_joint_individual_certainty` objects
#'
#' @param x An object of class `"pnp_joint_individual_certainty"`.
#' @param format Character string specifying cell format: `"both"` (default,
#'   `n (%)`), `"count"` (`n`), or `"prop"` (`%`).
#' @param digits Number of decimal places for percentages (default `1`).
#' @param mid_select Optional character vector of joint `mid` labels (from
#'   `x$meta$mid`), integer step indices, or numeric vector of length $K$
#'   selecting which joint `mid` threshold(s) to display as columns in the
#'   printed Evidence Matrix. When `NULL` (default), displays up to three
#'   representative steps (`min`, `middle`, and `max`).
#' @param ... Additional arguments passed to [insight::export_table()].
#'
#' @return Invisibly returns `x`.
#' @export
print.pnp_joint_individual_certainty <- function(
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

  selected_mids <- resolve_joint_mid_select(
    meta = x$meta,
    mid_select = mid_select,
    single = FALSE,
    default_all = FALSE
  )

  group_cols <- x$meta$group_cols %||% character(0)
  filtered_long <- dplyr::filter(
    x$meta$matrix_long,
    mid %in% selected_mids
  )

  wide_df <- format_joint_evidence_matrix_wide(
    matrix_long = filtered_long,
    mid_labels = selected_mids,
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
      title = "# Joint Clinical Significance Evidence Matrix",
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
        "# Joint Clinical Significance Evidence Matrix (",
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

#' Summary method for `pnp_joint_individual_certainty` objects
#'
#' @param object An object of class `"pnp_joint_individual_certainty"`.
#' @param format Character string passed to
#'   [print.pnp_joint_individual_certainty()].
#' @param digits Integer passed to [print.pnp_joint_individual_certainty()].
#' @param mid_select Optional selector passed to
#'   [print.pnp_joint_individual_certainty()].
#' @param ... Additional arguments passed to
#'   [print.pnp_joint_individual_certainty()].
#'
#' @return Invisibly returns `object`.
#' @export
summary.pnp_joint_individual_certainty <- function(
  object,
  format = c("both", "count", "prop"),
  digits = 1,
  mid_select = NULL,
  ...
) {
  print.pnp_joint_individual_certainty(
    object,
    format = format,
    digits = digits,
    mid_select = mid_select,
    ...
  )
}

#' Plot a 2D Joint Clinical Certainty or Sensitivity Plot for `pnp_joint_individual_certainty` objects
#'
#' @description
#' Visualizes joint individual-level Bayesian clinical certainty as either:
#' * **`type = "joint"` (default):** A 2D joint clinical scatter plot displaying
#'   posterior medians and credible intervals for two outcomes on the x- and
#'   y-axes, overlaid with a 2D ROPE rectangle via [geom_rope()] and an optional
#'   2D Region of Practical Significance target quadrant via [geom_rops()].
#' * **`type = "sensitivity"`:** A 2D clinical certainty sensitivity heatmap
#'   (`geom_tile()`) across joint `mid` steps on the x-axis and probability
#'   `thresholds` on the y-axis, where tile fill represents the proportion of
#'   individuals meeting joint clinical certainty.
#'
#' @param x An object of class `"pnp_joint_individual_certainty"`.
#' @param type Character string selecting the visualization type: `"joint"`
#'   (default) or `"sensitivity"`.
#' @param mid_select Optional selector choosing which joint `mid` step(s) from
#'   `x$meta$mid` to plot. Can be a character label (e.g., `"5.0 / 2.0"`), an
#'   integer step index (e.g., `1`), or a numeric vector of length $K$ matching
#'   a row of `x$meta$mid_grid`. For `type = "joint"`, must select a single step
#'   (defaults to the first step). For `type = "sensitivity"`, can select
#'   multiple steps (defaults to all steps).
#' @param outcomes_select Optional character or integer vector of length 2
#'   selecting which two outcomes from `x$meta$value_cols` to display on the x-
#'   and y-axes when `type = "joint"` and more than two outcomes were analyzed.
#'   Defaults to the first two outcomes in `x$meta$value_cols`.
#' @param color_by Variable used to color points and error bars when
#'   `type = "joint"`. Can be supplied as a character string or unquoted name.
#'   Must be one of `"signed_certainty"` (default), `"joint_pd"` (or `"pd"`),
#'   `"makowski_tier"`, `"joint_rope_pct"` (or `"rope_pct"`),
#'   `"joint_improvement_pct"` (or `"improvement_pct"`),
#'   `"joint_deterioration_pct"` (or `"deterioration_pct"`), or
#'   `"directed_certainty"`.
#' @param thresholds Numeric vector of probability thresholds in `(0, 1]` mapped
#'   to the y-axis when `type = "sensitivity"`. Defaults to
#'   `seq(0.50, 0.99, by = 0.01)`.
#' @param outcome Character string specifying which clinical direction(s) to
#'   display when `type = "sensitivity"`: `"both"` (default, mirrored butterfly
#'   heatmap), `"improvement"`, or `"deterioration"`.
#' @param color_palette Character string specifying the fill palette when
#'   `type = "sensitivity"`. Defaults to `"Temps"`.
#' @param show_ci Logical flag (`TRUE` by default). For `type = "joint"`, draws
#'   horizontal and vertical credible interval error bars when `TRUE`.
#' @param show_rops Logical flag (`TRUE` by default). For `type = "joint"`,
#'   highlights the 2D target quadrant of joint improvement via [geom_rops()].
#' @param point_size Numeric controlling point size in the joint plot (default
#'   `2`).
#' @param point_alpha Numeric in `[0, 1]` controlling point/tile opacity
#'   (default `1`).
#' @param interval_alpha Numeric in `[0, 1]` controlling error bar opacity in
#'   the joint plot (default `0.4`).
#' @param x_lab Optional string controlling the x-axis label.
#' @param y_lab Optional string controlling the y-axis label.
#' @param facet Logical flag (`TRUE` by default). When `TRUE` and `x` contains
#'   `.by` grouping columns in `x$meta$group_cols`, facets the plot via
#'   [ggplot2::facet_wrap()].
#' @param ... Additional arguments passed to [ggplot2::geom_point()] (for
#'   `type = "joint"`) or [ggplot2::geom_tile()] (for `type = "sensitivity"`).
#'
#' @return A [ggplot2::ggplot()] object.
#' @export
plot.pnp_joint_individual_certainty <- function(
  x,
  type = c("joint", "sensitivity"),
  mid_select = NULL,
  outcomes_select = NULL,
  color_by = c(
    "signed_certainty",
    "joint_pd",
    "makowski_tier",
    "joint_rope_pct",
    "joint_improvement_pct",
    "joint_deterioration_pct",
    "directed_certainty"
  ),
  thresholds = seq(0.50, 0.99, by = 0.01),
  outcome = c("both", "improvement", "deterioration"),
  color_palette = "Temps",
  show_ci = TRUE,
  show_rops = TRUE,
  point_size = 2,
  point_alpha = 1,
  interval_alpha = 0.4,
  x_lab = NULL,
  y_lab = NULL,
  facet = TRUE,
  ...
) {
  if (identical(type, c("joint", "sensitivity"))) {
    type <- "joint"
  }

  if (
    !is.character(type) ||
      length(type) != 1L ||
      is.na(type) ||
      !type %in% c("joint", "sensitivity")
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg type} must be either {.val joint} or {.val sensitivity}.",
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
      plot_joint_sensitivity(
        x = x,
        mid_select = mid_select,
        thresholds = thresholds,
        outcome = outcome,
        color_palette = color_palette,
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
    "joint_pd",
    "pd",
    "makowski_tier",
    "joint_rope_pct",
    "rope_pct",
    "joint_improvement_pct",
    "improvement_pct",
    "joint_deterioration_pct",
    "deterioration_pct",
    "directed_certainty"
  )

  color_quo <- rlang::enquo(color_by)
  color_expr <- rlang::quo_get_expr(color_quo)

  if (missing(color_by)) {
    color_col <- "signed_certainty"
  } else if (is.symbol(color_expr)) {
    color_col <- rlang::as_string(color_expr)
  } else if (is.character(color_by) && length(color_by) > 1L) {
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
        "x" = "Argument {.arg color_by} must be one of {.or {.val {unique(allowed_colors)}}}.",
        "i" = "Received {.val {rlang::as_label(color_quo)}}."
      ),
      class = "peoplenotp_error_bad_color_by"
    )
  }

  alias_map <- c(
    pd = "joint_pd",
    rope_pct = "joint_rope_pct",
    improvement_pct = "joint_improvement_pct",
    deterioration_pct = "joint_deterioration_pct"
  )
  if (color_col %in% names(alias_map)) {
    color_col <- unname(alias_map[[color_col]])
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

  if (!is.logical(show_rops) || length(show_rops) != 1L || is.na(show_rops)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg show_rops} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {show_rops}}."
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

  selected_mid <- resolve_joint_mid_select(
    meta = x$meta,
    mid_select = mid_select,
    single = TRUE,
    default_all = FALSE
  )

  pair_cols <- resolve_outcomes_select(x$meta$value_cols, outcomes_select)
  col_x <- pair_cols[[1]]
  col_y <- pair_cols[[2]]

  step_idx <- match(selected_mid, x$meta$mid)
  mid_x <- as.numeric(x$meta$mid_grid[step_idx, col_x])
  mid_y <- as.numeric(x$meta$mid_grid[step_idx, col_y])
  dir_x <- unname(x$meta$better_is[[col_x]])
  dir_y <- unname(x$meta$better_is[[col_y]])

  plot_df <- dplyr::filter(x$individual_summary, mid == selected_mid)

  med_x_col <- paste0("median_", col_x)
  low_x_col <- paste0("ci_low_", col_x)
  high_x_col <- paste0("ci_high_", col_x)

  med_y_col <- paste0("median_", col_y)
  low_y_col <- paste0("ci_low_", col_y)
  high_y_col <- paste0("ci_high_", col_y)

  if (is.null(x_lab)) {
    x_lab <- paste0(
      col_x,
      " (Posterior Median",
      if (isTRUE(show_ci)) " & CI)" else ")"
    )
  }
  if (is.null(y_lab)) {
    y_lab <- paste0(
      col_y,
      " (Posterior Median",
      if (isTRUE(show_ci)) " & CI)" else ")"
    )
  }

  color_labels <- c(
    signed_certainty = "Signed Certainty",
    joint_pd = "Joint PD",
    makowski_tier = "Evidence Tier",
    joint_rope_pct = "Joint % in ROPE",
    joint_improvement_pct = "Joint % Improvement",
    joint_deterioration_pct = "Joint % Deterioration",
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

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[med_x_col]],
      y = .data[[med_y_col]],
      color = .data[[color_col]]
    )
  )

  if (isTRUE(show_rops)) {
    p <- p +
      geom_rops(
        x_range = c(-mid_x, mid_x),
        y_range = c(-mid_y, mid_y),
        better_is = c(dir_x, dir_y)
      )
  }

  p <- p +
    geom_rope(
      x_range = c(-mid_x, mid_x),
      y_range = c(-mid_y, mid_y),
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted")

  if (isTRUE(show_ci)) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data[[low_x_col]], xmax = .data[[high_x_col]]),
        orientation = "y",
        width = 0,
        alpha = interval_alpha
      ) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[low_y_col]], ymax = .data[[high_y_col]]),
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

resolve_value_columns <- function(ungrouped_data, values_quo) {
  selected <- tryCatch(
    names(dplyr::select(ungrouped_data[0, , drop = FALSE], !!values_quo)),
    error = function(cnd) {
      cli::cli_abort(
        c(
          "x" = "Outcome column(s) {.code {rlang::as_label(values_quo)}} not found in {.arg data}.",
          "i" = "Available columns: {.val {names(ungrouped_data)}}."
        ),
        class = "peoplenotp_error_missing_column",
        parent = cnd
      )
    }
  )

  if (length(selected) < 2L) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg values} must select at least 2 outcome columns (selected {length(selected)}).",
        "i" = "Use {.fn individual_certainty} for single-outcome analyses."
      ),
      class = "peoplenotp_error_bad_values"
    )
  }

  selected
}

resolve_joint_better_is <- function(better_is, value_names) {
  k <- length(value_names)
  if (
    !is.character(better_is) ||
      !(length(better_is) %in% c(1L, k)) ||
      any(is.na(better_is)) ||
      !all(better_is %in% c("lower", "higher"))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg better_is} must be {.val lower} or {.val higher} (length 1 or length {k} matching {.arg values}).",
        "i" = "Received {.val {better_is}}."
      ),
      class = "peoplenotp_error_bad_better_is"
    )
  }

  if (length(better_is) == 1L) {
    out <- rep(better_is, k)
    names(out) <- value_names
    return(out)
  }

  if (!is.null(names(better_is)) && all(value_names %in% names(better_is))) {
    better_is <- better_is[value_names]
  } else {
    names(better_is) <- value_names
  }

  better_is
}

resolve_joint_mid <- function(mid, value_names) {
  k <- length(value_names)

  if (missing(mid) || is.null(mid)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid} must be a positive numeric scalar, a numeric vector of length {k}, or a list of length {k} of positive numeric vectors.",
        "i" = "Received {.val {if (missing(mid)) NULL else mid}}."
      ),
      class = "peoplenotp_error_bad_mid"
    )
  }

  is_valid_num <- function(v) {
    is.numeric(v) &&
      length(v) >= 1L &&
      !any(is.na(v)) &&
      !any(is.infinite(v)) &&
      all(v > 0)
  }

  if (is.numeric(mid)) {
    if (!is_valid_num(mid) || !(length(mid) %in% c(1L, k))) {
      cli::cli_abort(
        c(
          "x" = "Numeric {.arg mid} must contain positive finite values of length 1 or length {k} (matching the number of outcomes).",
          "i" = "To supply multiple MID steps across outcomes, pass a list of length {k} (e.g., {.code list(c(4, 5), c(1.5, 2.0))})."
        ),
        class = "peoplenotp_error_bad_mid"
      )
    }

    if (length(mid) == 1L) {
      mid_vec <- rep(mid, k)
    } else if (!is.null(names(mid)) && all(value_names %in% names(mid))) {
      mid_vec <- unname(mid[value_names])
    } else {
      mid_vec <- unname(mid)
    }

    mid_grid <- as.data.frame(
      matrix(mid_vec, nrow = 1L, ncol = k),
      stringsAsFactors = FALSE
    )
    names(mid_grid) <- value_names
  } else if (is.list(mid) && !is.data.frame(mid)) {
    if (length(mid) != k || !all(vapply(mid, is_valid_num, logical(1)))) {
      cli::cli_abort(
        c(
          "x" = "List {.arg mid} must have length {k} (one positive numeric vector per outcome).",
          "i" = "Received a list of length {length(mid)}."
        ),
        class = "peoplenotp_error_bad_mid"
      )
    }

    if (!is.null(names(mid)) && all(value_names %in% names(mid))) {
      mid <- mid[value_names]
    }

    lens <- lengths(mid)
    max_len <- max(lens)
    if (!all(lens %in% c(1L, max_len))) {
      cli::cli_abort(
        c(
          "x" = "All vectors in list {.arg mid} must have the same length (or length 1 to be recycled).",
          "i" = "Received element lengths: {.val {lens}}."
        ),
        class = "peoplenotp_error_bad_mid"
      )
    }

    cols_expanded <- lapply(mid, function(v) rep(v, length.out = max_len))
    names(cols_expanded) <- value_names
    mid_grid <- as.data.frame(cols_expanded, stringsAsFactors = FALSE)
    mid_grid <- dplyr::distinct(mid_grid)
  } else {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid} must be a positive numeric vector or a list of positive numeric vectors.",
        "i" = "Received an object of class {.cls {class(mid)}}."
      ),
      class = "peoplenotp_error_bad_mid"
    )
  }

  formatted_cols <- lapply(mid_grid, function(col_vals) {
    format(col_vals, trim = TRUE, nsmall = 1)
  })
  mid_labels <- do.call(
    paste,
    c(formatted_cols, list(sep = " / "))
  )

  list(
    mid_grid = tibble::as_tibble(mid_grid),
    mid_labels = mid_labels
  )
}

compute_joint_base <- function(df_sub, value_names, ci) {
  out_list <- list()
  for (vn in value_names) {
    draws <- df_sub[[vn]]
    qi <- tidybayes::median_qi(draws, .width = ci)
    pd_val <- max(
      mean(draws > 0, na.rm = TRUE),
      mean(draws < 0, na.rm = TRUE)
    )
    out_list[[paste0("median_", vn)]] <- qi$y
    out_list[[paste0("ci_low_", vn)]] <- qi$ymin
    out_list[[paste0("ci_high_", vn)]] <- qi$ymax
    out_list[[paste0("pd_", vn)]] <- pd_val
  }
  tibble::as_tibble(out_list)
}

compute_joint_mid_step <- function(
  df_sub,
  value_names,
  step_mids,
  step_label,
  better_is_vec,
  include_marginals = FALSE
) {
  n_draws <- nrow(df_sub)
  all_pd_dom <- rep(TRUE, n_draws)
  all_rope <- rep(TRUE, n_draws)
  all_imp <- rep(TRUE, n_draws)
  all_det <- rep(TRUE, n_draws)

  out_list <- list(mid = step_label)
  for (vn in value_names) {
    out_list[[paste0("mid_", vn)]] <- unname(step_mids[[vn]])
  }

  for (vn in value_names) {
    draws <- df_sub[[vn]]
    m <- step_mids[[vn]]
    dir <- better_is_vec[[vn]]

    pos_prop <- mean(draws > 0, na.rm = TRUE)
    neg_prop <- mean(draws < 0, na.rm = TRUE)
    dom_mask <- if (pos_prop >= neg_prop) draws > 0 else draws < 0

    rope_mask <- draws >= -m & draws <= m
    if (dir == "lower") {
      imp_mask <- draws < -m
      det_mask <- draws > m
    } else {
      imp_mask <- draws > m
      det_mask <- draws < -m
    }

    all_pd_dom <- all_pd_dom & dom_mask
    all_rope <- all_rope & rope_mask
    all_imp <- all_imp & imp_mask
    all_det <- all_det & det_mask

    if (isTRUE(include_marginals)) {
      out_list[[paste0("rope_pct_", vn)]] <- mean(rope_mask, na.rm = TRUE)
      out_list[[paste0("improvement_pct_", vn)]] <- mean(imp_mask, na.rm = TRUE)
      out_list[[paste0("deterioration_pct_", vn)]] <- mean(
        det_mask,
        na.rm = TRUE
      )
    }
  }

  joint_pd <- mean(all_pd_dom, na.rm = TRUE)
  joint_rope_pct <- mean(all_rope, na.rm = TRUE)
  joint_imp_pct <- mean(all_imp, na.rm = TRUE)
  joint_det_pct <- mean(all_det, na.rm = TRUE)

  dir_cert <- pmax(joint_imp_pct, joint_det_pct)
  signed_cert <- joint_imp_pct - joint_det_pct

  out_list$joint_pd <- joint_pd
  out_list$joint_rope_pct <- joint_rope_pct
  out_list$joint_improvement_pct <- joint_imp_pct
  out_list$joint_deterioration_pct <- joint_det_pct
  out_list$directed_certainty <- dir_cert
  out_list$signed_certainty <- signed_cert
  out_list$makowski_tier <- assign_makowski_tier(
    joint_rope_pct,
    joint_imp_pct,
    joint_det_pct
  )

  tibble::as_tibble(out_list)
}

format_joint_evidence_matrix_wide <- function(
  matrix_long,
  mid_labels,
  group_cols = character(0),
  format = "both",
  digits = 1
) {
  col_labels <- stats::setNames(
    paste0("MID = ", mid_labels),
    mid_labels
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
    mid_col = unname(col_labels[as.character(mid)]),
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

resolve_joint_mid_select <- function(
  meta,
  mid_select = NULL,
  single = FALSE,
  default_all = FALSE
) {
  available_mids <- meta$mid
  n_mids <- length(available_mids)
  k <- length(meta$value_cols)

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
    is.character(mid_select) &&
      length(mid_select) >= 1L &&
      !any(is.na(mid_select))
  ) {
    if (all(mid_select %in% available_mids)) {
      resolved <- unique(mid_select)
    }
  } else if (
    is.numeric(mid_select) &&
      !any(is.na(mid_select)) &&
      length(mid_select) >= 1L
  ) {
    if (length(mid_select) == k) {
      candidate_label <- paste(
        format(as.numeric(mid_select), trim = TRUE, nsmall = 1),
        collapse = " / "
      )
      if (candidate_label %in% available_mids) {
        resolved <- candidate_label
      }
    }
    if (
      is.null(resolved) &&
        all(mid_select == trunc(mid_select)) &&
        all(mid_select >= 1L & mid_select <= n_mids)
    ) {
      resolved <- unique(available_mids[as.integer(mid_select)])
    }
  }

  if (is.null(resolved) || (isTRUE(single) && length(resolved) != 1L)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg mid_select} must select {if (single) 'a single value' else 'value(s)'} present in `x$meta$mid` ({.val {available_mids}}).",
        "i" = "Received {.val {mid_select}}."
      ),
      class = "peoplenotp_error_bad_mid_select"
    )
  }

  resolved
}

resolve_outcomes_select <- function(value_cols, outcomes_select = NULL) {
  if (is.null(outcomes_select)) {
    return(value_cols[1:2])
  }

  if (
    is.character(outcomes_select) &&
      length(outcomes_select) == 2L &&
      !any(is.na(outcomes_select)) &&
      all(outcomes_select %in% value_cols) &&
      outcomes_select[[1]] != outcomes_select[[2]]
  ) {
    return(outcomes_select)
  }

  if (
    is.numeric(outcomes_select) &&
      length(outcomes_select) == 2L &&
      !any(is.na(outcomes_select)) &&
      all(outcomes_select %in% seq_along(value_cols)) &&
      outcomes_select[[1]] != outcomes_select[[2]]
  ) {
    return(value_cols[as.integer(outcomes_select)])
  }

  cli::cli_abort(
    c(
      "x" = "Argument {.arg outcomes_select} must specify 2 distinct outcomes from {.val {value_cols}}.",
      "i" = "Received {.val {outcomes_select}}."
    ),
    class = "peoplenotp_error_bad_outcomes_select"
  )
}

plot_joint_sensitivity <- function(
  x,
  mid_select = NULL,
  thresholds = seq(0.50, 0.99, by = 0.01),
  outcome = c("both", "improvement", "deterioration"),
  color_palette = "Temps",
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

  selected_mids <- resolve_joint_mid_select(
    meta = x$meta,
    mid_select = mid_select,
    single = FALSE,
    default_all = TRUE
  )

  fill_scale <- resolve_sensitivity_fill_scale(
    color_palette = color_palette,
    outcome = outcome
  )

  if (is.null(x_lab)) {
    x_lab <- paste0(
      "Joint MID (",
      paste(x$meta$value_cols, collapse = " / "),
      ")"
    )
  }
  if (is.null(y_lab)) {
    y_lab <- "Certainty Threshold"
  }

  group_cols <- x$meta$group_cols %||% character(0)
  thresh_vals <- sort(unique(thresholds))
  t_base <- min(thresh_vals)
  t_step <- if (length(thresh_vals) > 1L) min(diff(thresh_vals)) else 0.01
  t_span <- (max(thresh_vals) - t_base) + t_step

  adapted_summary <- dplyr::mutate(
    x$individual_summary,
    improvement_pct = joint_improvement_pct,
    deterioration_pct = joint_deterioration_pct
  )

  grid_df <- build_sensitivity_grid(
    individual_summary = adapted_summary,
    selected_mids = selected_mids,
    thresholds = thresh_vals,
    group_cols = group_cols,
    outcome = outcome,
    t_base = t_base,
    t_step = t_step
  )
  grid_df$mid <- factor(grid_df$mid, levels = selected_mids)

  if (outcome == "both") {
    x_annot <- (length(selected_mids) + 1) / 2
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
        label = "\u2191 Joint Improvement",
        label.r = ggplot2::unit(0, "lines"),
        size = 3.5
      ) +
      ggplot2::annotate(
        "label",
        x = x_annot,
        y = -y_annot,
        label = "\u2193 Joint Deterioration",
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
