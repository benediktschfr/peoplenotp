#' Custom ggplot2 theme for the peoplenotp project
#'
#' @description
#' A highly customized `ggplot2` theme tailored for publication-ready figures.
#' It removes distracting elements like facet background boxes, properly justifies
#' titles, and sets default palettes and geom colors natively.
#'
#' @param base_family Character. The font family to use. Default is `"Rubik"`.
#' @param base_size Numeric. The base font size in pt. Default is `11`.
#' @param continuous_scale Character vector of hex codes. The palette used for continuous data. Default is the `rcartocolor` "Geyser" palette.
#' @param discrete_scale Character vector of hex codes. The palette used for discrete data. Default is the custom "oasis" palette.
#' @param ink Character. The default foreground color for geoms (lines, points). Default is `"grey30"`.
#' @param ... Additional arguments passed to `ggplot2::theme_light()`.
#'
#' @return A `ggplot2::theme` object.
#' @export
#'
#' @examples
#' # Only run this example if the 'Rubik' font is installed on the system
#' if ("Rubik" %in% systemfonts::system_fonts()$family) {
#'   ggplot(mtcars, aes(wt, mpg, color = factor(cyl))) +
#'     geom_point() +
#'     theme_peoplenotp()
#' }
theme_peoplenotp <- function(
  base_family = "Rubik",
  base_size = 11,
  continuous_scale = rev(rcartocolor::carto_pal(100, "Geyser")),
  discrete_scale = custom_palettes[["oasis"]],
  ink = "grey30",
  ...
) {
  systemfonts::require_font(base_family, verbose = FALSE)
  continuous_scale <- continuous_scale
  accent_color <- discrete_scale[[1]]

  ggplot2::theme_light(
    base_family = base_family,
    base_size = base_size,
    accent = accent_color,
    ...
  ) +
    ggplot2::theme(
      # Title tuning
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.2)
      ),
      plot.title.position = "plot",

      # Axis title justification
      axis.title.x = ggplot2::element_text(hjust = 0),
      axis.title.y = ggplot2::element_text(hjust = 0),

      # Facet text
      strip.text = ggplot2::element_text(
        face = "bold",
        hjust = 0,
        colour = "grey10",
        size = ggplot2::rel(1.05)
      ),
      strip.background = ggplot2::element_blank(),

      # Guides
      legend.position = "bottom",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),

      # Default palettes
      palette.colour.discrete = unname(discrete_scale),
      palette.fill.discrete = unname(discrete_scale),
      palette.colour.continuous = continuous_scale,
      palette.fill.continuous = continuous_scale,

      geom = ggplot2::element_geom(
        ink = ink
      )
    )
}

#' Set the global theme to peoplenotp
#'
#' @description
#' A convenience function that sets `theme_peoplenotp()` as the global default
#' for all subsequent `ggplot2` plots in the current session.
#'
#' @export
#'
#' @examples
#' use_peoplenotp_theme()
use_peoplenotp_theme <- function() {
  ggplot2::set_theme(theme_peoplenotp())
}

#' Reset the global theme to default
#'
#' @description
#' Resets the global `ggplot2` theme back to the standard `theme_gray()`.
#'
#' @export
#'
#' @examples
#' use_default_theme()
use_default_theme <- function() {
  ggplot2::set_theme(ggplot2::theme_gray())
}
