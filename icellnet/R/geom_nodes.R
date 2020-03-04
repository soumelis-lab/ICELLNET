#' Display the nodes in the ICELLNET network.
#'
#' @description Display the nodes in the ICELLNET network.
#'
#' @param mapping mapping
#' @param data  data
#' @param position position
#' @param na.rm na.rm
#' @param show.legend to show the legend or not
#' @param inherit.aes  inherit.aes
#'
#' @examples
#' \dontrun{geom_nodes(aes(color = Cell_type), size = 30) }

geom_nodes = function (mapping = NULL,
                       data = NULL,
                       position = "identity",
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE,
                       ...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatNodes,
    geom = ggplot2::GeomPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm,
                  ...)
  )
}
StatNodes <-
  ggplot2::ggproto(
    "StatNodes",
    ggplot2::Stat,
    compute_layer = function(data, scales, params) {
      if (all(c("xend", "yend") %in% names(data))) {
        unique(subset(data, select = c(-xend,-yend)))
      } else {
        unique(data)
      }
    }
  )
StatEdges <-
  ggplot2::ggproto(
    "StatEdges",
    ggplot2::Stat,
    compute_layer = function(data, scales, params) {
      unique(subset(data,!(x == xend & y == yend)))
    }
  )
