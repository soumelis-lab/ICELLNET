#' Display the communication network with ICELLNET method
#'
#' @description Display the network from global communication scores matrix between a central cell (CC) and partner cell types (PC).
#' @details The matrix of global communication score is obtained with the icellnet.score function.
#' - out == display the communication score from the CC to PC
#' - in == display the communication score from the PC to CC
#' If you want to compare the communication between different conditions, we suggest to bind the different global
#' communication score matrix into a single one, so that all the scores are scaled together to display the
#' network.
#'
#' @param icn.score Matrix of global communication scores
#' @param scale  vector that contains the minimum and the maximum of the global communication score matrix.
#' @param direction Direction of the communication (either "out" or "in")
#' @param PC.col Color vector for the partner cells c("cell1"= "color1", "cell2"="color2").
#'
#' @export
#' @examples
#' \dontrun{network.create(icn.score = Scores.norm, scale = c(round(min(Scores.norm)),...
#' round(max(Scores.norm))), direction = "out", PC.col)}
#'
network.create = function(icn.score, scale = c(min(icn.score), max(icn.score)), direction = c("in", "out"), PC.col=NULL) {
  # Color palette
  if (is.null(PC.col)){
    my_pal = c(rep ("white", dim(icn.score)[1]+1))
    names(my_pal)= c(rownames(icn.score), colnames(icn.score))
  }else {
    CC.col = "grey"
    names(CC.col) = colnames(icn.score)
    my_pal = c(PC.col, CC.col)
  }

  # PC Nodes
  mtx.coord = data.frame(row.names = rownames(icn.score))
  mtx.coord$x = 0 + 0.2 * cos(seq(0, 2 * pi, length.out = (dim(icn.score)[1] +
                                                             1)))[-(dim(icn.score)[1] + 1)]
  mtx.coord$y = 0 + 0.2 * sin(seq(0, 2 * pi, length.out = (dim(icn.score)[1] +
                                                             1)))[-(dim(icn.score)[1] + 1)]
  mtx.coord$xend = mtx.coord$x
  mtx.coord$yend = mtx.coord$y
  mtx.coord$Cell_type = rownames(icn.score)
  # CC Node
  mtx.coord = rbind(mtx.coord, c(0, 0, 0, 0, colnames(icn.score)))
  rownames(mtx.coord)[dim(icn.score)[1] + 1] = colnames(icn.score)
  mtx.coord$Score = "NA"
  # Edges
  mtx.edges = data.frame(
    "x" = rep(0, dim(icn.score)[1]),
    "y" = rep(0, dim(icn.score)[1]),
    "xend" = 0 + 0.14 * cos(seq(0, 2 * pi, length.out =
                                  (
                                    dim(icn.score)[1] + 1
                                  )))[-(dim(icn.score)[1] + 1)],
    "yend" = 0 + 0.14 * sin(seq(0, 2 * pi, length.out =
                                  (
                                    dim(icn.score)[1] + 1
                                  )))[-(dim(icn.score)[1] + 1)],
    "Cell_type" = colnames(icn.score),
    "Score" = icn.score[, 1]
  )

  mtx.edges.2 = data.frame(
    "x" = 0 + 0.06 * cos(seq(0, 2 * pi, length.out = (
      dim(icn.score)[1] + 1
    )))[-(dim(icn.score)[1] + 1)],
    "y" = 0 + 0.06 * sin(seq(0, 2 * pi, length.out =
                               (
                                 dim(icn.score)[1] + 1
                               )))[-(dim(icn.score)[1] + 1)],
    "xend" = 0 + 0.14 * cos(seq(0, 2 * pi, length.out =
                                  (
                                    dim(icn.score)[1] + 1
                                  )))[-(dim(icn.score)[1] + 1)],
    "yend" = 0 + 0.14 * sin(seq(0, 2 * pi, length.out =
                                  (
                                    dim(icn.score)[1] + 1
                                  )))[-(dim(icn.score)[1] + 1)],
    "Cell_type" = colnames(icn.score),
    "Score" = icn.score[, 1]
  )

  # Network dataframe
  my_df_net = rbind(mtx.coord, mtx.edges)
  my_df_net$x = as.numeric(my_df_net$x)
  my_df_net$y = as.numeric(my_df_net$y)
  my_df_net$xend = as.numeric(my_df_net$xend)
  my_df_net$yend = as.numeric(my_df_net$yend)
  my_df_net$Score = as.numeric(my_df_net$Score)

  my_df_net_2 = rbind(mtx.coord, mtx.edges.2)
  my_df_net_2$x = as.numeric(my_df_net_2$x)
  my_df_net_2$y = as.numeric(my_df_net_2$y)
  my_df_net_2$xend = as.numeric(my_df_net_2$xend)
  my_df_net_2$yend = as.numeric(my_df_net_2$yend)
  my_df_net_2$Score = as.numeric(my_df_net_2$Score)

  if (direction =="in"){
    gg.net <-
      ggplot2::ggplot(my_df_net, ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      )) +
      geom_edges(
        data = my_df_net_2,
        ggplot2::aes(size = Score),
        arrow = grid::arrow(
          length = grid::unit(7, "pt"),
          angle = 35,
          type = "open",
          ends = "first"
        )
      ) +
      ggplot2::scale_size(
        range = c(1, 6),
        breaks = c(round(min(scale)), round((max(scale)-min(scale))/2), round(max(scale))),
        labels = c(round(min(scale)), round((max(scale)-min(scale))/2), round(max(scale))),
        limits = scale
      ) +
      ggplot2::geom_point(color = "black",
                 size = 31,
                 shape = 1) +
      geom_nodes(ggplot2::aes(color = Cell_type), size = 30) +
      ggplot2::scale_color_manual(values = my_pal, guide = "none") +
      ggplot2::scale_fill_manual(values = my_pal, guide = "none") +
      geom_node_label(ggplot2::aes(label = Cell_type, fill = Cell_type),
                      size = ggplot2::rel(6),
                      fontface = "bold") +
      ggplot2::expand_limits(x = c(-0.25, 0.25), y = c(-0.25, 0.25)) +
      ggplot2::theme_bw(base_size = 12, base_family = "") +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
      )
  }else if (direction =="out"){
    gg.net <-
      ggplot2::ggplot(my_df_net, ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      )) +
      geom_edges(
        data = my_df_net_2,
        ggplot2::aes(size = Score),
        arrow = grid::arrow(
          length = grid::unit(7, "pt"),
          angle = 35,
          type = "open",
          ends = "last"
        )
      ) +
      ggplot2::scale_size(
        range = c(1, 6),
        breaks = c(round(min(scale)), round((max(scale)-min(scale))/2), round(max(scale))),
        labels = c(round(min(scale)), round((max(scale)-min(scale))/2), round(max(scale))),
        limits = scale
      ) +
      ggplot2::geom_point(color = "black",
                 size = 31,
                 shape = 1) +
      geom_nodes(ggplot2::aes(color = Cell_type), size = 30) +
      ggplot2::scale_color_manual(values = my_pal, guide = "none") +
      ggplot2::scale_fill_manual(values = my_pal, guide = "none") +
      geom_node_label(ggplot2::aes(label = Cell_type, fill = Cell_type),
                      size = ggplot2::rel(6),
                      fontface = "bold") +
      ggplot2::expand_limits(x = c(-0.25, 0.25), y = c(-0.25, 0.25)) +
      ggplot2::theme_bw(base_size = 12, base_family = "") +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
      )
  } else
    stop('Error : Direction of the communication ("in" or "out") must be specified ')
  return(gg.net)
}
