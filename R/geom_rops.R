#' Add a ROPS (Region of Practical Significance) Rectangle to a Plot
#'
#' @description
#' Draws a ROPS in a `ggplot2` plot. This region highlights the area where
#' an effect is considered practically significant, based on the boundaries of
#' a ROPE (Region of Practical Equivalence) and the desired direction of the effect.
#' Supports 1D (vertical/horizontal bands) and 2D (target quadrant) regions.
#'
#' @param x_range Numeric vector of length 2. Boundaries of the ROPE on the x-axis.
#' @param y_range Numeric vector of length 2. Boundaries of the ROPE on the y-axis.
#' @param better_is Character vector indicating the desired direction of the effect for the x and y axes. Allowed values are `"higher"` or `"lower"`. Default is `c("higher", "higher")`. If a single value is provided, it is applied to the available axis.
#' @param alpha Numeric. Transparency of the fill color (0 to 1). Default is `0.2`.
#' @param fill Character. Fill color of the rectangle. Default is `"#a6d96a"` (light green).
#' @param ... Additional arguments passed to `ggplot2::annotate()` (e.g., `color = NA`).
#'
#' @return A `ggplot2::annotate` layer that can be added to a plot via `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' # Base plot
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#'
#' # 1D ROPS on the x-axis (higher is better)
#' p + geom_rops(x_range = c(-1, 1), better_is = "higher")
#'
#' # 2D ROPS: x should be higher, y should be lower (Bottom-Right Quadrant)
#' p +
#'   geom_rope(x_range = c(-1, 1), y_range = c(-2, 2)) + # Add ROPE for context
#'   geom_rops(
#'     x_range = c(-1, 1),
#'     y_range = c(-2, 2),
#'     better_is = c("higher", "lower")
#'   )
geom_rops <- function(
  x_range = NULL,
  y_range = NULL,
  better_is = c("higher", "higher"),
  alpha = 0.2,
  fill = "#a6d96a",
  ...
) {
  # 1. Error handling for missing ranges
  if (is.null(x_range) && is.null(y_range)) {
    cli::cli_abort(c(
      "x" = "Neither {.arg x_range} nor {.arg y_range} was supplied.",
      "i" = "Please provide at least one of the two ranges to calculate the ROPS."
    ))
  }

  # 2. Error handling for better_is input
  if (!all(better_is %in% c("higher", "lower"))) {
    cli::cli_abort(c(
      "x" = "{.arg better_is} contains invalid values.",
      "i" = "Only {.val higher} or {.val lower} are allowed."
    ))
  }

  # 3. Set default boundaries to "infinity"
  x_min <- -Inf
  x_max <- Inf
  y_min <- -Inf
  y_max <- Inf

  # 4. Calculate X boundaries based on direction
  if (!is.null(x_range)) {
    if (better_is[1] == "higher") {
      x_min <- x_range[2] # Starts at upper ROPE bound, goes to Inf
    } else {
      x_max <- x_range[1] # Starts at -Inf, ends at lower ROPE bound
    }
  }

  # 5. Calculate Y boundaries based on direction
  if (!is.null(y_range)) {
    # Match the correct direction for y (fallback to index 1 if only a single value was passed)
    y_dir <- if (length(better_is) > 1) better_is[2] else better_is[1]

    if (y_dir == "higher") {
      y_min <- y_range[2] # Starts at upper ROPE bound, goes to Inf
    } else {
      y_max <- y_range[1] # Starts at -Inf, ends at lower ROPE bound
    }
  }

  # 6. Return the complete ggplot annotation object
  ggplot2::annotate(
    "rect",
    xmin = x_min,
    xmax = x_max,
    ymin = y_min,
    ymax = y_max,
    fill = fill,
    alpha = alpha,
    ...
  )
}
