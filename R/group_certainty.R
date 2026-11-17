#' Evaluate group-level Bayesian clinical certainty and heterogeneity
#'
#' @description
#' Summarizes population-level (fixed) effects and inter-individual
#' random-effects heterogeneity (`sd__*` and `cor__*`) from a fitted Bayesian
#' model (such as a `brmsfit` object). For each fixed effect, computes the
#' Posterior Median, Credible Interval (`CI_low`, `CI_high`), Probability of
#' Direction (`pd`), and the proportion of posterior draws inside the Region of
#' Practical Equivalence (`rope_pct`).
#'
#' @param model A fitted Bayesian model object supported by `{insight}` and
#'   `{tidybayes}` (typically a `brmsfit` object), or a data frame of posterior
#'   draws.
#' @param rope Numeric vector of length 2 specifying the Region of Practical
#'   Equivalence, e.g., `c(-2.5, 2.5)`.
#' @param parameters Optional character vector or regular expression pattern to
#'   filter target fixed-effect parameters. Defaults to `NULL` (all fixed
#'   effects).
#' @param ci Credible interval level for parameter point-and-interval summaries
#'   (default `0.95`).
#' @param ... Additional arguments (currently unused).
#'
#' @return An S3 object of class `c("pnp_group_certainty", "data.frame")`
#'   containing columns `parameter`, `effect`, `median`, `CI_low`, `CI_high`,
#'   `pd`, and `rope_pct`.
#'
#' @examples
#' \dontrun{
#' library(brms)
#' data(burnout)
#'
#' fit <- brm(eei ~ time_weeks * group + (1 + time_weeks | id), data = burnout)
#' group_certainty(fit, rope = c(-10, 10))
#' }
#'
#' @export
group_certainty <- function(
  model,
  rope,
  parameters = NULL,
  ci = 0.95,
  ...
) {
  is_valid_model <- insight::is_model(model) ||
    inherits(model, "brmsfit") ||
    is.data.frame(model)

  if (!is_valid_model) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg model} must be a valid fitted model object (e.g., {.cls brmsfit}) or a data frame of posterior draws.",
        "i" = "Received an object of class {.cls {class(model)}}."
      ),
      class = "peoplenotp_error_bad_model"
    )
  }

  if (
    missing(rope) ||
      !is.numeric(rope) ||
      length(rope) != 2L ||
      any(is.na(rope)) ||
      rope[1] > rope[2]
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg rope} must be a numeric vector of length 2 with `rope[1] <= rope[2]`.",
        "i" = "Received {.val {if (missing(rope)) NULL else rope}}."
      ),
      class = c("peoplenotp_error_bad_rope", "peoplenotp_error_bad_range")
    )
  }

  if (
    !is.null(parameters) &&
      (!is.character(parameters) ||
        length(parameters) == 0L ||
        any(is.na(parameters)) ||
        any(!nzchar(parameters)))
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg parameters} must be a non-empty character vector or {.code NULL}.",
        "i" = "Received {.val {parameters}}."
      ),
      class = "peoplenotp_error_bad_parameters"
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

  all_vars <- if (is.data.frame(model)) {
    setdiff(names(model), c(".chain", ".iteration", ".draw"))
  } else {
    tidybayes::get_variables(model)
  }

  fixed_vars <- grep("^b_", all_vars, value = TRUE)
  if (length(fixed_vars) == 0L && is.data.frame(model)) {
    fixed_vars <- grep(
      "^(?!sd_|cor_|r_|lp__|lprior)",
      all_vars,
      perl = TRUE,
      value = TRUE
    )
  }

  if (!is.null(parameters)) {
    clean_fixed <- sub("^b_", "", fixed_vars)
    matched_mask <- fixed_vars %in%
      parameters |
      clean_fixed %in% parameters |
      Reduce(
        `|`,
        lapply(parameters, function(pat) {
          grepl(pat, fixed_vars) | grepl(pat, clean_fixed)
        }),
        init = rep(FALSE, length(fixed_vars))
      )
    fixed_vars <- fixed_vars[matched_mask]
  }

  if (length(fixed_vars) == 0L) {
    cli::cli_abort(
      c(
        "x" = "No fixed-effect parameters in {.arg model} matched the requested selection.",
        "i" = "Available variables include: {.val {utils::head(all_vars, 10)}}."
      ),
      class = "peoplenotp_error_no_matching_variables"
    )
  }

  re_vars <- grep("^(sd_|cor_)", all_vars, value = TRUE)

  fixed_draws_df <- extract_group_draws(
    model = model,
    vars = fixed_vars,
    default_pattern = "^b_.*"
  )

  fixed_rows <- lapply(fixed_vars, function(var_name) {
    draws <- fixed_draws_df$.value[fixed_draws_df$.variable == var_name]
    summarize_group_fixed_draws(
      draws = draws,
      parameter = var_name,
      rope = rope,
      ci = ci
    )
  })
  fixed_df <- dplyr::bind_rows(fixed_rows)

  if (length(re_vars) > 0L) {
    re_draws_df <- extract_group_draws(
      model = model,
      vars = re_vars,
      default_pattern = "^(sd_|cor_).*"
    )
    re_rows <- lapply(re_vars, function(var_name) {
      draws <- re_draws_df$.value[re_draws_df$.variable == var_name]
      summarize_group_random_draws(
        draws = draws,
        parameter = var_name,
        ci = ci
      )
    })
    random_df <- dplyr::bind_rows(re_rows)
    combined_df <- dplyr::bind_rows(fixed_df, random_df)
  } else {
    random_df <- data.frame(
      parameter = character(0),
      effect = character(0),
      median = numeric(0),
      CI_low = numeric(0),
      CI_high = numeric(0),
      pd = numeric(0),
      rope_pct = numeric(0),
      stringsAsFactors = FALSE
    )
    combined_df <- fixed_df
  }

  out <- as.data.frame(combined_df, stringsAsFactors = FALSE)
  rownames(out) <- NULL

  structure(
    out,
    rope = rope,
    ci = ci,
    random_effects = as.data.frame(
      random_df[, c("parameter", "median", "CI_low", "CI_high"), drop = FALSE],
      stringsAsFactors = FALSE
    ),
    class = c("pnp_group_certainty", "data.frame")
  )
}

#' Print method for `pnp_group_certainty` objects
#'
#' @param x An object of class `"pnp_group_certainty"`.
#' @param digits Number of decimal places for numeric values (default `2`).
#' @param ... Additional arguments passed to [insight::export_table()].
#'
#' @return Invisibly returns `x`.
#' @export
print.pnp_group_certainty <- function(x, digits = 2, ...) {
  fixed_cols <- c(
    "parameter",
    "median",
    "CI_low",
    "CI_high",
    "pd",
    "rope_pct"
  )
  random_cols <- c("parameter", "median", "CI_low", "CI_high")

  df_base <- as.data.frame(x, stringsAsFactors = FALSE)
  fixed_df <- df_base[df_base$effect == "fixed", fixed_cols, drop = FALSE]
  random_df <- df_base[df_base$effect == "random", random_cols, drop = FALSE]
  rownames(fixed_df) <- NULL
  rownames(random_df) <- NULL

  if (nrow(random_df) > 0L) {
    formatted <- insight::export_table(
      list(fixed_df, random_df),
      digits = digits,
      protect_integers = FALSE,
      title = list(
        "# Group-Level Clinical Certainty (Fixed Effects)",
        "# Individual Heterogeneity (Random Effects)"
      ),
      ...
    )
  } else {
    formatted <- insight::export_table(
      fixed_df,
      digits = digits,
      protect_integers = FALSE,
      title = "# Group-Level Clinical Certainty (Fixed Effects)",
      ...
    )
  }

  cat(formatted, "\n", sep = "")
  invisible(x)
}

#' @export
`$.pnp_group_certainty` <- function(x, name) {
  val <- .subset2(x, name)
  if (!is.null(val)) {
    return(val)
  }
  alias_map <- c(
    ci_low = "CI_low",
    ci_high = "CI_high",
    lower = "CI_low",
    upper = "CI_high",
    Parameter = "parameter",
    Median = "median"
  )
  if (name %in% names(alias_map)) {
    return(.subset2(x, alias_map[[name]]))
  }
  NULL
}

#' Extract posterior draws for a set of model variables
#'
#' @param model A model object or data frame of draws.
#' @param vars Character vector of variable names to extract.
#' @param default_pattern Regex pattern used with `tidybayes::gather_draws()`.
#'
#' @return A data frame with `.variable` and `.value` columns.
#' @noRd
extract_group_draws <- function(model, vars, default_pattern) {
  if (is.data.frame(model)) {
    draws_long <- tidyr::pivot_longer(
      dplyr::select(model, dplyr::all_of(vars)),
      cols = dplyr::everything(),
      names_to = ".variable",
      values_to = ".value"
    )
    return(draws_long)
  }

  pattern_sym <- rlang::sym(default_pattern)
  draws_long <- tidybayes::gather_draws(
    model,
    !!pattern_sym,
    regex = TRUE
  )
  dplyr::filter(draws_long, .variable %in% vars)
}

#' Summarize fixed-effect posterior draws for `group_certainty()`
#'
#' @param draws Numeric vector of posterior draws.
#' @param parameter Parameter name.
#' @param rope Numeric vector of length 2.
#' @param ci Credible interval level.
#'
#' @return A single-row tibble.
#' @noRd
summarize_group_fixed_draws <- function(
  draws,
  parameter,
  rope,
  ci
) {
  qi <- tidybayes::median_qi(draws, .width = ci)
  pd_val <- max(
    mean(draws > 0, na.rm = TRUE),
    mean(draws < 0, na.rm = TRUE)
  )
  rope_pct <- mean(draws >= rope[1] & draws <= rope[2], na.rm = TRUE)

  tibble::tibble(
    parameter = parameter,
    effect = "fixed",
    median = qi$y,
    CI_low = qi$ymin,
    CI_high = qi$ymax,
    pd = pd_val,
    rope_pct = rope_pct
  )
}

#' Summarize random-effect posterior draws for `group_certainty()`
#'
#' @param draws Numeric vector of posterior draws.
#' @param parameter Parameter name.
#' @param ci Credible interval level.
#'
#' @return A single-row tibble.
#' @noRd
summarize_group_random_draws <- function(draws, parameter, ci) {
  qi <- tidybayes::median_qi(draws, .width = ci)

  tibble::tibble(
    parameter = parameter,
    effect = "random",
    median = qi$y,
    CI_low = qi$ymin,
    CI_high = qi$ymax,
    pd = NA_real_,
    rope_pct = NA_real_
  )
}
