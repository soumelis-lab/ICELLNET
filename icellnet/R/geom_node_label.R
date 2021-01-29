#' Label the node with the name of the cell type in the ICELLNET network.
#'
#' @description Label the node with the name of the cell type in the ICELLNET network.
#'
#' @param mapping mapping
#' @param data  data
#' @param position position
#' @param parse parse
#' @param nudge_x nudge_x
#' @param nudge_y nudge_y
#' @param label.padding label.padding
#' @param label.r label.r
#' @param label.size label.size
#' @param na.rm na.rm
#' @param show.legend to show the legend or not
#' @param inherit.aes  inherit.aes
#' @param ...
#'
#' @examples
#' \dontrun{geom_node_label(aes(label = Cell_type, fill = Cell_type), size = rel(6), fontface = "bold") }

geom_node_label = function (mapping = NULL,
                            data = NULL,
                            position = "identity",
                            ...,
                            parse = FALSE,
                            nudge_x = 0,
                            nudge_y = 0,
                            label.padding = unit(0.25, "lines"),
                            label.r = unit(0.15, "lines"),
                            label.size = 0.25,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`",
           call. = FALSE)
    }
    position <- ggplot2::position_nudge(nudge_x, nudge_y)
  }
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatNodes,
    geom = GeomLabel,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      parse = parse,
      label.padding = label.padding,
      label.r = label.r,
      label.size = label.size,
      na.rm = na.rm,
      ...
    )
  )
}
