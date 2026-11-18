test_that("geom_rope() and geom_rops() create ggplot2 rect annotations and validate inputs", {
  p_base <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()

  p_rope <- p_base + geom_rope(x_range = c(2, 4), y_range = c(15, 25))
  expect_s3_class(p_rope, "ggplot")
  expect_length(p_rope$layers, 2L)

  p_rops <- p_base +
    geom_rops(
      x_range = c(2, 4),
      y_range = c(15, 25),
      better_is = c("higher", "lower")
    )
  expect_s3_class(p_rops, "ggplot")
  expect_length(p_rops$layers, 2L)

  expect_error(geom_rope())
  expect_error(geom_rops())
  expect_error(geom_rops(x_range = c(-1, 1), better_is = "invalid"))
})

test_that("highlight_individual() adds highlight layers and theme_peoplenotp() returns a theme", {
  df <- tibble::tibble(
    id = rep(1:3, each = 4),
    time = rep(1:4, 3),
    val = seq_len(12)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(time, val, group = id)) +
    ggplot2::geom_line() +
    highlight_individual(target_id = 2, id_col = id, geom = "both")

  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 3L)

  thm <- theme_peoplenotp(base_family = "sans")
  expect_s3_class(thm, "theme")
})
