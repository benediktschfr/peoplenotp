#' Calculate Clinical Certainty Sensitivity
#'
#' @description
#' Performs a sensitivity analysis on posterior draws over a grid of ROPE
#' (Region of Practical Equivalence) values and certainty thresholds.
#' Evaluates the robustness of clinical outcomes for either improvement or deterioration.
#' Requires a data frame of posterior draws (e.g., generated via `tidybayes::add_epred_draws()`).
#'
#' @param data A data frame containing the posterior draws.
#' @param col The column containing the posterior values to analyze (e.g., `.epred` or `difference`).
#'   Can be passed unquoted. Defaults to `.epred`.
#' @param id_col The column identifying individuals. Can be passed unquoted. Defaults to `id`.
#' @param rope_range A numeric vector of length 2 specifying the min and max values
#'   for the ROPE sequence to test. Defaults to `c(10, 15)`.
#' @param rope_step Numeric step size for the ROPE sequence. Defaults to `1`.
#' @param threshold_range A numeric vector of length 2 specifying the min and max
#'   certainty thresholds. Defaults to `c(0.50, 0.99)`.
#' @param threshold_step Numeric step size for the thresholds. Defaults to `0.1`.
#' @param outcome Character. Either `"improvement"` or `"deterioration"`. Determines
#'   which percentage is evaluated. Defaults to `"improvement"`.
#' @param better_is Character. Either `"higher"` or `"lower"`. Passed to `get_clinical_certainty()`. Defaults to `"higher"`.
#'
#' @return A summarized data frame (tibble) ready for visualization.
#'
#' @importFrom dplyr ungroup mutate select left_join summarise n_distinct
#' @importFrom tidyr nest crossing unnest
#' @importFrom purrr map2
#' @export
get_clinical_sensitivity_grid <- function(
  data,
  col = .epred,
  id_col = id,
  rope_range = c(10, 15),
  rope_step = 1,
  threshold_range = c(0.50, 0.99),
  threshold_step = 0.1,
  outcome = c("improvement", "deterioration"),
  better_is = "higher"
) {
  # 1. Argument-Check
  outcome <- match.arg(outcome)

  # 2. Sequenzen generieren
  rope_seq <- seq(rope_range[1], rope_range[2], by = rope_step)
  threshold_seq <- seq(
    threshold_range[1],
    threshold_range[2],
    by = threshold_step
  )

  # 3. Den übergebenen Datensatz (der schon Draws enthält!) gruppieren/verschachteln.
  # Wir nutzen {{ id_col }}, um die Spalte dynamisch auszuwerten.
  nested_draws <- data |>
    dplyr::ungroup() |>
    tidyr::nest(.by = {{ id_col }})

  # 4. Interne Helferfunktion
  # Wir leiten die Zielspalte mit {{ col }} direkt an get_clinical_certainty weiter
  calc_cc_helper <- function(df, current_rope) {
    get_clinical_certainty(
      df,
      col = {{ col }},
      better_is = better_is,
      range = c(-current_rope, current_rope)
    )
  }

  # 5. Berechnung über alle ROPE-Werte (purrr::map2)
  sensitivity_results <- nested_draws |>
    tidyr::crossing(rope = rope_seq) |>
    dplyr::mutate(
      results = purrr::map2(
        data,
        rope,
        \(a, b) calc_cc_helper(a, b),
        .progress = TRUE
      )
    ) |>
    tidyr::unnest(results) |>
    dplyr::select(-data)

  # 6. Ziel-Spalte basierend auf dem Argument 'outcome' wählen
  target_col <- if (outcome == "improvement") {
    "improvement_pct"
  } else {
    "deterioration_pct"
  }

  # 7. Aggregation über alle Thresholds
  sensitivity_grid <- tidyr::crossing(
    rope = rope_seq,
    threshold = threshold_seq
  ) |>
    dplyr::left_join(
      sensitivity_results,
      by = "rope",
      relationship = "many-to-many"
    ) |>
    dplyr::summarise(
      n_total = dplyr::n_distinct({{ id_col }}),
      n_success = sum(.data[[target_col]] >= threshold, na.rm = TRUE),
      prop_certain = n_success / n_total,
      .by = c(rope, threshold)
    )

  sensitivity_grid
}
