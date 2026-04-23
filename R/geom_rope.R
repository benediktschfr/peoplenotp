#' Add a ROPE (Region of Practical Equivalence) Rectangle to a Plot
#'
#' @description
#' Draws a ROPE in a `ggplot2` plot. Supports ROPEs on the x-axis
#' (vertical band), on the y-axis (horizontal band), or on both axes
#' simultaneously (2D rectangle, e.g., for joint analyses).
#'
#' @param x_range Numeric vector of length 2 (e.g., `c(-1, 1)`). Boundaries of the ROPE on the x-axis.
#' @param y_range Numeric vector of length 2. Boundaries of the ROPE on the y-axis.
#' @param alpha Numeric. Transparency of the fill color (0 to 1). Default is `0.2`.
#' @param fill Character. Fill color of the rectangle. Default is `"grey70"`.
#' @param ... Additional arguments passed to `ggplot2::annotate()` (e.g., `color = "black"`, `linetype = "dashed"`).
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
#' # 1D ROPE on the x-axis
#' p + geom_rope(x_range = c(2.5, 3.5))
#'
#' # 1D ROPE on the y-axis
#' p + geom_rope(y_range = c(15, 20), fill = "blue", alpha = 0.1)
#'
#' # 2D ROPE (Joint ROPE Box) with dashed border
#' p + geom_rope(
#'   x_range = c(2.5, 3.5),
#'   y_range = c(15, 20),
#'   color = "black",
#'   linetype = "dashed"
#' )
geom_rope <- function(
  x_range = NULL,
  y_range = NULL,
  alpha = 0.2,
  fill = "grey70",
  ...
) {
  # 1. Error handling with cli for a clean console output
  if (is.null(x_range) && is.null(y_range)) {
    cli::cli_abort(c(
      "x" = "Neither {.arg x_range} nor {.arg y_range} was supplied.",
      "i" = "Please provide at least one of the two ranges to draw a ROPE."
    ))
  }

  # 2. Set default boundaries to "infinity"
  x_min <- -Inf
  x_max <- Inf
  y_min <- -Inf
  y_max <- Inf

  # 3. Overwrite boundaries if specific ranges were provided
  if (!is.null(x_range)) {
    x_min <- x_range[[1]]
    x_max <- x_range[[2]]
  }

  if (!is.null(y_range)) {
    y_min <- y_range[[1]]
    y_max <- y_range[[2]]
  }

  # 4. Return the complete ggplot annotation object
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
