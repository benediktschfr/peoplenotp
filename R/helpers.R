#' Weist den klinischen Status basierend auf Wahrscheinlichkeiten zu
#'
#' @description
#' Interne Hilfsfunktion, die von den Hauptfunktionen genutzt wird, um
#' Code-Wiederholungen bei der Zuweisung der Kategorien zu vermeiden.
#'
#' @param imp_pct Wahrscheinlichkeit für Verbesserung
#' @param det_pct Wahrscheinlichkeit für Verschlechterung
#' @param bound_likely Schwellenwert für "Likely"
#' @param bound_signif Schwellenwert für "Significant"
#'
#' @return Ein Character-Vektor mit der klinischen Kategorie.
#' @noRd
assign_clinical_status <- function(
  imp_pct,
  det_pct,
  bound_likely,
  bound_signif
) {
  dplyr::case_when(
    imp_pct > bound_signif ~ "Significant Improvement",
    imp_pct > bound_likely ~ "Likely Improvement",
    det_pct > bound_signif ~ "Significant Deterioration",
    det_pct > bound_likely ~ "Likely Deterioration",
    TRUE ~ "Undecided"
  )
}
