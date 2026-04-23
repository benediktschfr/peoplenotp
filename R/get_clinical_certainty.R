#' Calculate Clinical Certainty Based on ROPE
#'
#' @description
#' This function evaluates posterior draws (typically the difference between
#' two measurement points) and categorizes the effect into clinical statuses.
#' It calculates the median, credible interval, and Probability of Direction (pd).
#' It also compares the entire posterior distribution against a defined
#' *Region of Practical Equivalence (ROPE)*. Depending on the desired direction
#' of the effect (`better_is`) and predefined probability thresholds
#' (`certainty_boundaries`), a final `clinical_status` is assigned.
#'
#' @param data A data frame or tibble containing the posterior draws.
#' @param col The unquoted name of the column in `data` containing the numeric posterior draws.
#' @param better_is Character. Indicates the desired direction of the effect. Allowed values are `"higher"` (higher values are better, e.g., well-being) or `"lower"` (lower values are better, e.g., depression). Default is `"higher"`.
#' @param range Numeric vector of length 2. Defines the boundaries of the ROPE (Minimum Clinically Important Difference). Default is `c(-1, 1)`.
#' @param ci_width Numeric value between 0 and 1. Defines the width of the credible interval returned for the median effect. Default is `0.95`.
#' @param certainty_boundaries Numeric vector of length 2. Defines the probability thresholds for assigning clinical categories in the order `c(Likely_Boundary, Significant_Boundary)`. Default is `c(0.975, 0.99)`.
#' @param include_rope Logical. Should the ROPE boundaries and the percentage of draws within the ROPE be included in the output? Default is `TRUE`.
#'
#' @return A \code{tibble} with a single row containing the median, credible intervals, probability of direction (pd), ROPE percentages (optional), and the final \code{clinical_status}.
#'
#' @export
#'
#' @examples
#' # Simulate a dataframe of posterior draws
#' set.seed(123)
#' df_patient <- tibble::tibble(
#'   iteration = 1:4000,
#'   score_difference = rnorm(4000, mean = 5, sd = 1.5)
#' )
#'
#' # Use the function in a pipe workflow
#' df_patient |>
#'   get_clinical_certainty(
#'     col = score_difference,
#'     better_is = "higher",
#'     range = c(-2, 2)
#'   )
get_clinical_certainty <- function(
  data,
  col,
  better_is = "higher",
  range = c(-1, 1),
  ci_width = 0.95,
  certainty_boundaries = c(0.975, 0.99),
  include_rope = TRUE
) {
  # 1. Spalte ziehen
  posterior_draws <- dplyr::pull(data, {{ col }})

  # 2. Probability of Direction (pd) berechnen
  # Entspricht dem Maximum aus der Wahrscheinlichkeit für Werte > 0 und Werte < 0
  pd_val <- max(mean(posterior_draws > 0), mean(posterior_draws < 0))

  # 3. Wahrscheinlichkeiten manuell berechnen (Base R, extrem schnell)
  if (better_is == "higher") {
    imp_pct <- mean(posterior_draws > range[2])
    det_pct <- mean(posterior_draws < range[1])
  } else {
    imp_pct <- mean(posterior_draws < range[1])
    det_pct <- mean(posterior_draws > range[2])
  }

  # 4. Status ermitteln (via utils.R Helper)
  status <- assign_clinical_status(
    imp_pct = imp_pct,
    det_pct = det_pct,
    bound_likely = certainty_boundaries[[1]],
    bound_signif = certainty_boundaries[[2]]
  )

  # 5. Credible Interval berechnen (tidybayes)
  effect <- tidybayes::median_qi(posterior_draws, .width = ci_width)

  # 6. Basis-Output generieren
  res <- tibble::tibble(
    median = effect$y,
    lower = effect$ymin,
    upper = effect$ymax,
    pd = pd_val,
    improvement_pct = imp_pct,
    deterioration_pct = det_pct,
    clinical_status = status
  )

  # 7. ROPE-Spalten bedingt hinzufügen
  if (include_rope) {
    rope_pct <- mean(posterior_draws >= range[1] & posterior_draws <= range[2])

    # Fügt die ROPE-Spalten direkt vor improvement_pct ein,
    # damit die Reihenfolge logisch bleibt
    res <- dplyr::mutate(
      res,
      rope_low = range[1],
      rope_high = range[2],
      rope_pct = rope_pct,
      .before = improvement_pct
    )
  }

  res
}
