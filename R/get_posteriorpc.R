#' Plot a posterior predictive check for a Bayesian model
#'
#' @description
#' Generates a density-based posterior predictive check (`ggplot2` object) for a
#' fitted Bayesian model (such as a `brmsfit` object). Extracts the underlying
#' model data and response variable via [insight::get_data()] and
#' [insight::find_response()], draws simulated outcome predictions from the
#' posterior predictive distribution via [tidybayes::add_predicted_draws()],
#' plots each draw as a semi-transparent density curve, and overlays the
#' empirical density of the observed outcome variable.
#'
#' If the model formula includes an autoregressive structure (e.g., `ar()` or
#' `arma()`), the outcome variable is automatically removed from the prediction
#' input data before calling [tidybayes::add_predicted_draws()] (accompanied by
#' an informational note), unless `keep_outcome = TRUE`.
#'
#' @param model A fitted Bayesian model object supported by `{insight}` and
#'   `{tidybayes}` (typically a `brmsfit` object).
#' @param ndraws A single positive integer specifying the number of posterior
#'   predictive draws to sample and plot, or `NULL` to use all available draws.
#'   Defaults to `300`.
#' @param alpha A single numeric value between `0` and `1` controlling the
#'   opacity of the individual posterior predictive draw density lines.
#'   Defaults to `0.2`.
#' @param accent_color A single character string specifying the color (hex code
#'   or R color name) of the empirical observed outcome density curve.
#'   Defaults to `"#E76F51"`.
#' @param linewidth A single positive numeric value specifying the line width of
#'   the empirical observed outcome density curve. Defaults to `2`.
#' @param keep_outcome A single logical flag (`TRUE` or `FALSE`). When `FALSE`
#'   (the default) and an autoregressive term (`ar()` or `arma()`) is detected
#'   in `model`, the observed outcome column is excluded from the model data
#'   before generating predicted draws, and an informational message is emitted.
#'   Set to `TRUE` to force keeping the outcome column in the prediction data.
#'
#' @return A [ggplot2::ggplot()] object displaying the posterior predictive
#'   density draws alongside the observed outcome density.
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
#' get_posteriorpc(fit, ndraws = 300)
#' }
#'
#' @export
get_posteriorpc <- function(
  model,
  ndraws = 300,
  alpha = 0.2,
  accent_color = "#E76F51",
  linewidth = 2,
  keep_outcome = FALSE
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

  if (
    !is.null(ndraws) &&
      (!is.numeric(ndraws) ||
        length(ndraws) != 1L ||
        is.na(ndraws) ||
        is.infinite(ndraws) ||
        ndraws < 1 ||
        ndraws != trunc(ndraws))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg ndraws} must be a single positive integer or {.code NULL}.",
        "i" = "Received {.val {ndraws}}."
      ),
      class = "peoplenotp_error_bad_ndraws"
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
    !is.character(accent_color) ||
      length(accent_color) != 1L ||
      is.na(accent_color) ||
      !nzchar(accent_color)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg accent_color} must be a single non-empty character string.",
        "i" = "Received {.val {accent_color}}."
      ),
      class = "peoplenotp_error_bad_accent_color"
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

  if (
    !is.logical(keep_outcome) ||
      length(keep_outcome) != 1L ||
      is.na(keep_outcome)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg keep_outcome} must be a single non-missing logical value (`TRUE` or `FALSE`).",
        "i" = "Received {.val {keep_outcome}}."
      ),
      class = "peoplenotp_error_bad_keep_outcome"
    )
  }

  model_data <- insight::get_data(model, verbose = FALSE)
  outcome_col <- insight::find_response(model)

  if (
    is.null(outcome_col) ||
      !is.data.frame(model_data) ||
      !outcome_col %in% names(model_data)
  ) {
    cli::cli_abort(
      c(
        "x" = "Could not extract the outcome variable or data from {.arg model}."
      ),
      class = "peoplenotp_error_missing_outcome"
    )
  }

  pred_data <- model_data
  if (has_ar_structure(model) && !keep_outcome) {
    cli::cli_inform(
      c(
        "i" = "Autoregressive structure detected in {.arg model}: removing outcome column {.field {outcome_col}} before drawing predictions.",
        "i" = "Set {.code keep_outcome = TRUE} to include the outcome column when predicting."
      ),
      class = "peoplenotp_message_ar_outcome_removed"
    )
    pred_data <- dplyr::select(pred_data, -dplyr::all_of(outcome_col))
  }

  pred_draws <- pred_data |>
    tidybayes::add_predicted_draws(model, ndraws = ndraws)

  ggplot2::ggplot(pred_draws, ggplot2::aes(x = .prediction)) +
    ggplot2::geom_line(
      stat = "density",
      ggplot2::aes(group = .draw),
      alpha = alpha
    ) +
    ggplot2::geom_density(
      data = model_data,
      ggplot2::aes(x = .data[[outcome_col]]),
      color = accent_color,
      linewidth = linewidth
    )
}

#' Detect whether a model contains an autoregressive structure via `{insight}`
#'
#' @param model A fitted model object.
#' @return A single logical (`TRUE` or `FALSE`).
#' @noRd
has_ar_structure <- function(model) {
  formulas <- insight::find_formula(model, verbose = FALSE)
  if (is.null(formulas)) {
    return(FALSE)
  }

  ar_fns <- c("ar", "arma", "cor_ar", "cor_arma")
  formula_names <- unlist(
    lapply(formulas, function(f) {
      if (inherits(f, "formula")) {
        c(all.names(f), all.names(attr(f, "autocor")))
      } else {
        character(0L)
      }
    }),
    use.names = FALSE
  )

  any(ar_fns %in% formula_names)
}
