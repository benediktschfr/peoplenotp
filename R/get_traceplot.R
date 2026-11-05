#' Plot MCMC trace plots for key structural parameters of a Bayesian model
#'
#' @description
#' Generates faceted MCMC trace plots for a fitted Bayesian model (such as a
#' `brmsfit` object) while omitting individual-level random effect trajectories
#' (`r_*`) by default for visual clarity.
#'
#' When `regex = NULL` (the default), parameters are organized into distinct
#' structural groups controlled by logical flags:
#' - `fixed`: Population-level (fixed) effects (`^b_`)
#' - `sd`: Group-level standard deviations (`^sd_`)
#' - `sigma`: Residual or distributional scale parameters (`^sigma`)
#' - `ar`: Autoregressive parameters (`^(ar(\\[|$)|cor_ar)`)
#'
#' Any enabled group that is present in `model` receives its own faceted
#' `ggplot2` trace plot section, and all active sections are combined into a
#' unified layout with a shared chain legend via [patchwork::wrap_plots()].
#' Groups that are not present in `model` (e.g., `sd` in a non-hierarchical
#' model or `ar` in a model without autocorrelation) are silently skipped.
#'
#' Supplying a custom regular expression via `regex` overrides the four logical
#' flags and creates a single faceted trace plot for all matching parameters.
#'
#' @param model A fitted Bayesian model object supported by `{insight}` and
#'   `{tidybayes}` (typically a `brmsfit` object).
#' @param fixed A single logical flag (`TRUE` or `FALSE`). If `TRUE` (the
#'   default), includes population-level effects matching `"^b_"`.
#' @param sd A single logical flag (`TRUE` or `FALSE`). If `TRUE` (the default),
#'   includes group-level standard deviations matching `"^sd_"`.
#' @param sigma A single logical flag (`TRUE` or `FALSE`). If `TRUE` (the
#'   default), includes residual scale parameters matching `"^sigma"`.
#' @param ar A single logical flag (`TRUE` or `FALSE`). If `TRUE` (the default),
#'   includes autoregressive parameters matching `"^(ar(\\[|$)|cor_ar)"`.
#' @param regex An optional single character string specifying a custom regular
#'   expression to select model parameters, or `NULL` (the default) to use the
#'   `fixed`, `sd`, `sigma`, and `ar` flags. When provided, `regex` overrides
#'   the individual logical flags.
#' @param alpha A single numeric value between `0` and `1` controlling the
#'   opacity of the trace lines. Defaults to `1`.
#' @param linewidth A single positive numeric value specifying the line width of
#'   the trace lines. Defaults to `0.5`.
#' @param scales A single character string passed to [ggplot2::facet_wrap()]
#'   controlling axis scaling across facets. Must be one of `"free_y"` (the
#'   default), `"fixed"`, `"free_x"`, or `"free"`.
#' @param facet_ncol An optional single positive integer passed to the `ncol`
#'   argument of [ggplot2::facet_wrap()] within each parameter section, or
#'   `NULL` (the default) for automatic facet grid dimensions.
#' @param patchwork_ncol An optional single positive integer passed to the
#'   `ncol` argument of [patchwork::wrap_plots()] to control how parameter
#'   sections are arranged, or `NULL`. Defaults to `1` (vertical stacking).
#'
#' @return A `patchwork` / `ggplot` object displaying the MCMC trace plots for
#'   the selected model parameters.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' data(misophonia)
#'
#' fit <- brm(
#'   severity_score ~ day * phase + ar(time = day),
#'   data = misophonia
#' )
#'
#' # Default: fixed effects, SDs, sigma, and AR stacked in sections
#' get_traceplot(fit)
#'
#' # Only fixed effects and sigma, with custom styling
#' get_traceplot(fit, sd = FALSE, ar = FALSE, alpha = 0.7, facet_ncol = 2)
#'
#' # Custom regex override
#' get_traceplot(fit, regex = "^b_.*")
#' }
#'
#' @export
get_traceplot <- function(
  model,
  fixed = TRUE,
  sd = TRUE,
  sigma = TRUE,
  ar = TRUE,
  regex = NULL,
  alpha = 1,
  linewidth = 0.5,
  scales = "free_y",
  facet_ncol = NULL,
  patchwork_ncol = 1
) {
  if (!insight::is_model(model)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg model} must be a valid fitted model object supported by {.pkg insight}.",
        "i" = "Received an object of class {.cls {class(model)}}."
      ),
      class = "peoplenotp_error_bad_model"
    )
  }

  validate_logical_flag(fixed, "fixed")
  validate_logical_flag(sd, "sd")
  validate_logical_flag(sigma, "sigma")
  validate_logical_flag(ar, "ar")

  if (
    !is.null(regex) &&
      (!is.character(regex) ||
        length(regex) != 1L ||
        is.na(regex) ||
        !nzchar(regex))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg regex} must be a single non-empty character string or {.code NULL}.",
        "i" = "Received {.val {regex}}."
      ),
      class = "peoplenotp_error_bad_regex"
    )
  }

  if (
    !is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      alpha < 0 ||
      alpha > 1
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg alpha} must be a single numeric value between 0 and 1.",
        "i" = "Received {.val {alpha}}."
      ),
      class = "peoplenotp_error_bad_alpha"
    )
  }

  if (
    !is.numeric(linewidth) ||
      length(linewidth) != 1L ||
      is.na(linewidth) ||
      is.infinite(linewidth) ||
      linewidth <= 0
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg linewidth} must be a single positive numeric value.",
        "i" = "Received {.val {linewidth}}."
      ),
      class = "peoplenotp_error_bad_linewidth"
    )
  }

  allowed_scales <- c("free_y", "fixed", "free_x", "free")
  if (
    !is.character(scales) ||
      length(scales) != 1L ||
      is.na(scales) ||
      !scales %in% allowed_scales
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg scales} must be one of {.or {.val {allowed_scales}}}.",
        "i" = "Received {.val {scales}}."
      ),
      class = "peoplenotp_error_bad_scales"
    )
  }

  validate_positive_int_or_null(facet_ncol, "facet_ncol")
  validate_positive_int_or_null(patchwork_ncol, "patchwork_ncol")

  if (is.null(regex) && !any(c(fixed, sd, sigma, ar))) {
    cli::cli_abort(
      c(
        "x" = "At least one parameter group ({.arg fixed}, {.arg sd}, {.arg sigma}, or {.arg ar}) must be {.code TRUE} when {.arg regex} is {.code NULL}.",
        "i" = "Enable at least one parameter flag or provide a custom {.arg regex}."
      ),
      class = "peoplenotp_error_no_parameter_groups"
    )
  }

  all_vars <- tidybayes::get_variables(model)

  if (!is.null(regex)) {
    active_patterns <- stats::setNames(
      regex,
      "custom"
    )
  } else {
    group_patterns <- c(
      fixed = "^b_",
      sd = "^sd_",
      sigma = "^sigma",
      ar = "^(ar(\\[|$)|cor_ar)"
    )
    enabled_flags <- c(
      fixed = fixed,
      sd = sd,
      sigma = sigma,
      ar = ar
    )
    active_patterns <- group_patterns[enabled_flags]
  }

  matched_patterns <- purrr::keep(active_patterns, function(pat) {
    any(grepl(pat, all_vars))
  })

  if (length(matched_patterns) == 0L) {
    cli::cli_abort(
      c(
        "x" = "No variables in {.arg model} matched the requested parameter selection ({.val {unname(active_patterns)}}).",
        "i" = "Available variables include: {.val {utils::head(all_vars, 10)}}."
      ),
      class = "peoplenotp_error_no_matching_variables"
    )
  }

  plot_list <- purrr::map(unname(matched_patterns), function(pat) {
    build_single_traceplot(
      model = model,
      pattern = pat,
      alpha = alpha,
      linewidth = linewidth,
      scales = scales,
      facet_ncol = facet_ncol
    )
  })

  if (length(plot_list) == 1L) {
    return(
      patchwork::wrap_plots(
        plot_list,
        ncol = patchwork_ncol,
        guides = "collect"
      )
    )
  }

  patchwork::wrap_plots(
    plot_list,
    ncol = patchwork_ncol,
    guides = "collect"
  )
}

#' Build a single faceted trace plot for one regex pattern
#'
#' @param model A fitted Bayesian model object.
#' @param pattern A single regex character string.
#' @param alpha Line opacity.
#' @param linewidth Line width.
#' @param scales Facet scales option.
#' @param facet_ncol Facet column count or `NULL`.
#' @return A `ggplot` object.
#' @noRd
build_single_traceplot <- function(
  model,
  pattern,
  alpha,
  linewidth,
  scales,
  facet_ncol
) {
  pattern_sym <- rlang::sym(pattern)

  draws <- tidybayes::gather_draws(
    model,
    !!pattern_sym,
    regex = TRUE
  )

  ggplot2::ggplot(draws, ggplot2::aes(x = .iteration, y = .value)) +
    ggplot2::geom_line(
      ggplot2::aes(color = factor(.chain), group = .chain),
      alpha = alpha,
      linewidth = linewidth
    ) +
    ggplot2::facet_wrap(
      ~.variable,
      scales = scales,
      ncol = facet_ncol
    ) +
    ggplot2::labs(color = "Chain")
}

#' Validate a single non-missing logical flag
#'
#' @param x Value to check.
#' @param arg_name Name of the argument.
#' @noRd
validate_logical_flag <- function(x, arg_name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg {arg_name}} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {x}}."
      ),
      class = "peoplenotp_error_bad_flag"
    )
  }
}

#' Validate a single positive integer or NULL
#'
#' @param x Value to check.
#' @param arg_name Name of the argument.
#' @noRd
validate_positive_int_or_null <- function(x, arg_name) {
  if (
    !is.null(x) &&
      (!is.numeric(x) ||
        length(x) != 1L ||
        is.na(x) ||
        is.infinite(x) ||
        x < 1 ||
        x != trunc(x))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg {arg_name}} must be a single positive integer or {.code NULL}.",
        "i" = "Received {.val {x}}."
      ),
      class = paste0("peoplenotp_error_bad_", arg_name)
    )
  }
}
