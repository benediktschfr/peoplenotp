#' Highlight Specific Individuals using Tidy Selection
#'
#' @description
#' A flexible `ggplot2` helper to highlight specific individuals. Supports tidy evaluation,
#' allowing you to pass the ID column as a bare name, a string, or a tidyselect helper.
#'
#' @param target_id Vector of IDs to highlight.
#' @param id_col The column identifying individuals. Can be a bare name (e.g., `id`)
#'   or a string (`"id"`). Defaults to `id`.
#' @param geom Character. One of `"point"`, `"line"`, or `"both"`. Defaults to `"point"`.
#' @param highlight_col Color for highlighted data. Defaults to the oasis palette.
#' @param base_col Color for background data. Defaults to `"grey30"`.
#' @param size Numeric. Point size for the top layer. Defaults to `5`.
#' @param linewidth Numeric. Line width for the top layer. Defaults to `1.5`.
#' @param ... Additional arguments passed on to \code{geom_point()} and/or \code{geom_line()}
#'   (e.g., \code{alpha}, \code{shape}, \code{linetype}, or \code{position}). Note that if
#'   \code{geom = "both"}, these arguments are passed to both geometries.
#'
#' @return A list of `ggplot2` components (aesthetics, scales, and layers) that
#'   can be added to a `ggplot` object.
#'
#' @importFrom ggplot2 aes scale_color_manual geom_line geom_point
#' @importFrom rlang ensym as_string
#' @export
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' # Create a dummy dataset
#' set.seed(42)
#' df <- data.frame(
#'   proband_nr = rep(1:4, each = 10),
#'   time = rep(1:10, 4),
#'   value = as.vector(replicate(4, cumsum(rnorm(10))))
#' )
#'
#' # Example 1: Use a bare name for id_col and default points
#' ggplot(df, aes(x = time, y = value)) +
#'   geom_point(color = "grey70") +
#'   highlight_individual(target_id = 2, id_col = proband_nr)
#'
#' # Example 2: Use a string for id_col, highlight multiple IDs, and draw both
#' # Note: Remember to include group = proband_nr in the main ggplot call for lines!
#' ggplot(df, aes(x = time, y = value, group = proband_nr)) +
#'   geom_line(color = "grey70") +
#'   highlight_individual(
#'     target_id = c(1, 3),
#'     id_col = "proband_nr",
#'     geom = "both",
#'     highlight_col = "blue"
#'   )
#'
#' # Example 3: Using the `...` argument to pass position_jitter and alpha
#' ggplot(df, aes(x = time, y = value)) +
#'   geom_point(color = "grey70") +
#'   highlight_individual(
#'     target_id = 4,
#'     id_col = proband_nr,
#'     geom = "point",
#'     alpha = 0.8,
#'     position = position_jitter(width = 0.2, height = 0)
#'   )
#' }
highlight_individual <- function(
  target_id,
  id_col = id,
  geom = c("point", "line", "both"),
  highlight_col = peoplenotp::custom_palettes$oasis[[2]],
  base_col = "grey30",
  size = 5,
  linewidth = 1.5,
  ...
) {
  geom <- match.arg(geom)
  id_sym <- rlang::ensym(id_col)
  id_name <- rlang::as_string(id_sym)

  layers <- list(
    ggplot2::aes(color = !!id_sym %in% target_id),
    ggplot2::scale_color_manual(
      values = c("TRUE" = highlight_col, "FALSE" = base_col),
      guide = "none"
    )
  )

  top_layer_data <- function(plot_data) {
    plot_data[plot_data[[id_name]] %in% target_id, ]
  }

  if (geom %in% c("line", "both")) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_line(
          data = top_layer_data,
          color = highlight_col,
          linewidth = linewidth,
          ...
        )
      )
    )
  }

  if (geom %in% c("point", "both")) {
    layers <- append(
      layers,
      list(
        ggplot2::geom_point(
          data = top_layer_data,
          color = highlight_col,
          size = size,
          ...
        )
      )
    )
  }

  layers
}
