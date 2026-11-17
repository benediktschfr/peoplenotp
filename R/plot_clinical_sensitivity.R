#' Plot clinical certainty sensitivity across ROPEs and probability thresholds
#'
#' @description
#' Visualizes the summarized sensitivity grid produced by
#' [get_clinical_sensitivity_grid()], showing how the proportion of the cohort
#' achieving clinical certainty (`prop_certain`) decays as the Minimal Important
#' Difference (`rope`) widens across varying certainty `threshold` levels.
#'
#' @param grid_data A data frame or tibble containing `rope`, `threshold`, and
#'   `prop_certain` (as returned by [get_clinical_sensitivity_grid()]).
#' @param color_palette Character string specifying the discrete color palette.
#'   Supports any `{viridisLite}` option (`"viridis"`, `"magma"`, `"inferno"`,
#'   `"plasma"`, `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`) or any palette
#'   name in [custom_palettes] (e.g., `"oasis"`, `"maritime"`). Defaults to
#'   `"viridis"`.
#' @param ... Additional arguments passed to [ggplot2::geom_line()] and
#'   [ggplot2::geom_point()].
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @examples
#' set.seed(123)
#' draws_df <- tibble::tibble(
#'   id = rep(1:20, each = 100),
#'   .epred = stats::rnorm(2000, mean = 12, sd = 3)
#' )
#'
#' grid <- get_clinical_sensitivity_grid(
#'   draws_df,
#'   rope_range = c(8, 14),
#'   rope_step = 2,
#'   threshold_range = c(0.8, 0.95),
#'   threshold_step = 0.05
#' )
#'
#' plot_clinical_sensitivity(grid)
#'
#' @export
plot_clinical_sensitivity <- function(
  grid_data,
  color_palette = "viridis",
  ...
) {
  if (!is.data.frame(grid_data)) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg grid_data} must be a data frame or tibble.",
        "i" = "Received an object of class {.cls {class(grid_data)}}."
      ),
      class = "peoplenotp_error_bad_data"
    )
  }

  required_cols <- c("rope", "threshold", "prop_certain")
  missing_cols <- setdiff(required_cols, names(grid_data))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg grid_data} is missing required column{?s}: {.val {missing_cols}}.",
        "i" = "Required columns are {.val {required_cols}}."
      ),
      class = "peoplenotp_error_missing_column"
    )
  }

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
      !color_palette %in% c(viridis_options, custom_options)
  ) {
    cli::cli_abort(
      c(
        "x" = "Argument {.arg color_palette} must be a valid viridis option or a palette in {.code custom_palettes}.",
        "i" = "Received {.val {color_palette}}."
      ),
      class = "peoplenotp_error_bad_color_palette"
    )
  }

  color_scale <- if (color_palette %in% viridis_options) {
    ggplot2::scale_color_viridis_d(option = color_palette)
  } else {
    ggplot2::scale_color_manual(
      values = unname(peoplenotp::custom_palettes[[color_palette]])
    )
  }

  ggplot2::ggplot(
    grid_data,
    ggplot2::aes(
      x = rope,
      y = prop_certain,
      color = as.factor(threshold),
      group = as.factor(threshold)
    )
  ) +
    ggplot2::geom_line(...) +
    ggplot2::geom_point(...) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      limits = c(0, 1)
    ) +
    color_scale +
    ggplot2::labs(
      x = "ROPE / Minimal Important Difference (MID)",
      y = "Proportion Meeting Certainty",
      color = "Certainty Threshold"
    ) +
    ggplot2::theme_minimal()
}
