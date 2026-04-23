#' Calculate Joint Clinical Certainty for Two Outcomes
#'
#' @description
#' This function calculates the joint probability of clinical success (or deterioration)
#' across two co-primary endpoints. By evaluating the condition row-wise (per MCMC draw),
#' it naturally accounts for the correlation between the two posterior distributions.
#' It calculates the marginal Probability of Direction (pd) for each outcome, as well
#' as the `joint_pd` (the probability that both effects are simultaneously in their
#' respective dominant directions). Optionally, marginal clinical certainties for each
#' individual outcome can be returned.
#'
#' @param data A data frame or tibble containing the posterior draws.
#' @param col_1 The unquoted name of the first numeric column.
#' @param col_2 The unquoted name of the second numeric column.
#' @param better_is Character vector of length 2. Indicates the desired direction for `col_1` and `col_2`. Allowed values are `"higher"` or `"lower"`. Default is `c("higher", "higher")`.
#' @param ranges List of two numeric vectors. Defines the ROPE boundaries for `col_1` and `col_2`. Default is `list(c(-1, 1), c(-1, 1))`.
#' @param ci_width Numeric. Width of the credible interval returned for the medians. Default is `0.95`.
#' @param certainty_boundaries Numeric vector of length 2. Thresholds for "Likely" and "Significant". Default is `c(0.975, 0.99)`.
#' @param include_marginals Logical. If \code{TRUE}, the output will also include the marginal probabilities and clinical statuses for each individual column. Default is \code{FALSE}.
#'
#' @return A \code{tibble} with a single row containing the medians, CIs, marginal and joint Probability of Direction (pd), joint ROPE probabilities, and the final \code{joint_clinical_status}. If \code{include_marginals = TRUE}, individual ROPE probabilities and statuses are appended.
#'
#' @export
#'
#' @examples
#' # Simulate correlated posterior draws
#' set.seed(123)
#' df_patient <- tibble::tibble(
#'   wellbeing = rnorm(4000, mean = 5, sd = 2),       # higher is better
#'   depression = wellbeing * -0.8 + rnorm(4000, 0, 1) # lower is better
#' )
#'
#' # Default: Only joint certainties, pds, and base descriptives
#' df_patient |>
#'   get_joint_clinical_certainty(
#'     col_1 = wellbeing,
#'     col_2 = depression,
#'     better_is = c("higher", "lower"),
#'     ranges = list(c(-2, 2), c(-1, 1))
#'   )
#'
#' # Extended: Include marginal certainties for plotting/summaries
#' df_patient |>
#'   get_joint_clinical_certainty(
#'     col_1 = wellbeing,
#'     col_2 = depression,
#'     better_is = c("higher", "lower"),
#'     ranges = list(c(-2, 2), c(-1, 1)),
#'     include_marginals = TRUE
#'   )
get_joint_clinical_certainty <- function(
  data,
  col_1,
  col_2,
  better_is = c("higher", "higher"),
  ranges = list(c(-1, 1), c(-1, 1)),
  ci_width = 0.95,
  certainty_boundaries = c(0.975, 0.99),
  include_marginals = FALSE
) {
  # 1. Spalten ziehen (Tidy Eval)
  draws1 <- dplyr::pull(data, {{ col_1 }})
  draws2 <- dplyr::pull(data, {{ col_2 }})

  # 2. Probability of Direction (pd) marginal und joint berechnen
  # pd_1
  pd_1_pos <- mean(draws1 > 0)
  pd_1_neg <- mean(draws1 < 0)
  pd_1 <- max(pd_1_pos, pd_1_neg)
  mask_pd_1 <- if (pd_1_pos >= pd_1_neg) draws1 > 0 else draws1 < 0

  # pd_2
  pd_2_pos <- mean(draws2 > 0)
  pd_2_neg <- mean(draws2 < 0)
  pd_2 <- max(pd_2_pos, pd_2_neg)
  mask_pd_2 <- if (pd_2_pos >= pd_2_neg) draws2 > 0 else draws2 < 0

  # joint_pd (Draw-weise Evaluierung, ob BEIDE auf ihrer dominanten Seite liegen)
  joint_pd <- mean(mask_pd_1 & mask_pd_2)

  # 3. Hilfsfunktion zur Bestimmung von Verbesserung/Verschlechterung pro Draw
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

  # 4. Bedingungen für alle Draws auswerten
  flags1 <- get_flags(draws1, better_is[1], ranges[[1]])
  flags2 <- get_flags(draws2, better_is[2], ranges[[2]])

  # 5. Wahrscheinlichkeiten berechnen (Marginal & Joint)
  imp_pct_1 <- mean(flags1$imp)
  det_pct_1 <- mean(flags1$det)

  imp_pct_2 <- mean(flags2$imp)
  det_pct_2 <- mean(flags2$det)

  joint_imp_pct <- mean(flags1$imp & flags2$imp)
  joint_det_pct <- mean(flags1$det & flags2$det)

  # 6. Status Zuweisungen
  bound_likely <- certainty_boundaries[[1]]
  bound_signif <- certainty_boundaries[[2]]

  status_1 <- assign_clinical_status(
    imp_pct_1,
    det_pct_1,
    bound_likely,
    bound_signif
  )
  status_2 <- assign_clinical_status(
    imp_pct_2,
    det_pct_2,
    bound_likely,
    bound_signif
  )
  joint_clinical_status <- assign_clinical_status(
    joint_imp_pct,
    joint_det_pct,
    bound_likely,
    bound_signif
  )

  # 7. Basis-Statistiken generieren
  e1 <- tidybayes::median_qi(draws1, .width = ci_width)
  e2 <- tidybayes::median_qi(draws2, .width = ci_width)

  # 8. Finalen Datensatz zusammenbauen
  res <- tibble::tibble(
    # Spalte 1 (Basis deskriptiv)
    median_1 = e1$y,
    lower_1 = e1$ymin,
    upper_1 = e1$ymax,
    pd_1 = pd_1,
    # Spalte 1 (Marginal Certainties)
    improvement_pct_1 = imp_pct_1,
    deterioration_pct_1 = det_pct_1,
    clinical_status_1 = status_1,

    # Spalte 2 (Basis deskriptiv)
    median_2 = e2$y,
    lower_2 = e2$ymin,
    upper_2 = e2$ymax,
    pd_2 = pd_2,
    # Spalte 2 (Marginal Certainties)
    improvement_pct_2 = imp_pct_2,
    deterioration_pct_2 = det_pct_2,
    clinical_status_2 = status_2,

    # Joint Certainties
    joint_pd = joint_pd,
    joint_improvement_pct = joint_imp_pct,
    joint_deterioration_pct = joint_det_pct,
    joint_clinical_status = joint_clinical_status
  )

  # 9. Output filtern, falls Marginals nicht gewünscht sind
  if (!include_marginals) {
    res <- res |>
      dplyr::select(
        median_1,
        lower_1,
        upper_1,
        pd_1,
        median_2,
        lower_2,
        upper_2,
        pd_2,
        dplyr::starts_with("joint_")
      )
  }

  res
}
