#' Display the edges in the ICELLNET network.
#'
#' @description Display the edges in the ICELLNET network.
#'
#' @param mapping mapping
#' @param data  data
#' @param position position
#' @param arrow arrow
#' @param curvature curvature
#' @param angle angle
#' @param ncp ncp
#' @param na.rm na.rm
#' @param show.legend to show the legend or not
#' @param inherit.aes  inherit.aes
#'
#' @examples
#' \dontrun{geom_edges (data = my_df_net_2, aes(size = Score), arrow = arrow( length = unit(7, "pt"), angle = 35, type = "open", ends = "first") }
#'
geom_edges = function (mapping = NULL,
data = NULL,
position = "identity",
arrow = NULL,
curvature = 0,
angle = 90,
ncp = 5,
na.rm = FALSE,
show.legend = NA,
inherit.aes = TRUE,
...) {
  if (!curvature) {
    geom = ggplot2::GeomSegment
    params = list(arrow = arrow, na.rm = na.rm, ...)
  }
  else {
    geom = ggplot2::GeomCurve
    params = list(
      arrow = arrow,
      curvature = curvature,
      angle = angle,
      ncp = ncp,
      na.rm = na.rm,
      ...
    )
  }
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatEdges,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = params
  )
}
