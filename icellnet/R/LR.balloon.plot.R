#'
#' Plot of individual communication scores of contributing L/R pairs.
#'
#' @description Display the individual communication scores of ligand/receptor pairs that contributes with a score superior to a threshold (thresh = 10 as a default) and/or the top n interactions contributing to the score. Return a ggplot object.
#'
#' @param lr Matrix of individual communication scores
#' @param thresh Value set as a threshold to display the L/R pairs contributing with a score superior to this threshold
#' @param topn Value set as top n interactions to display
#' @param db.name.couple output of the name.lr.couple() function. name.lr.couple(db, "Family") is set as a default
#' @param title Title of the balloon plot
#' @param family.col Color vector for the family of molecules c("family1"= "color1", "family2"="color2")
#
#'
#' @export
#' @examples
#' \dontrun{lr <- icellnet.score.results[[2]]
#' LR.balloon.plot(lr = lr, thresh = 20, top_n=NULL, db.name.couple=name.lr.couple(db, "Family"))
#' }
#'
#'
LR.balloon.plot <- function (lr = lr, thresh = 0, topn = NULL, db.name.couple = db.name.couple,
          title = title, family.col = family.col)
{
  interactions = rownames(lr)
  lr = as.data.frame(lr)
  lr$Pair = interactions
  lr = as.data.frame(lr[stats::complete.cases(lr), ])
  lr = lr %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. > 0)) %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. >= thresh))
  if (!is.null(topn)) {
    if (topn > dim(lr)[1]) {
      note(paste0("lr contains only ", dim(lr)[1], " after filtering interaction highest than theshold"))
      lr = dplyr::left_join(lr, as.data.frame(db.name.couple),
                            by = "Pair")

    }
    else {
      lr = lr %>% dplyr::mutate(sum = rowSums(dplyr::across(where(is.numeric)))) %>%
        dplyr::arrange(dplyr::desc(sum)) %>% dplyr::top_n(topn,
                                                          sum)
      lr = dplyr::left_join(lr, as.data.frame(db.name.couple),
                            by = "Pair") %>% dplyr::select(-sum)
    }
  } else {
    lr = dplyr::left_join(lr, as.data.frame(db.name.couple), by = "Pair")
  }
  if (!(is.null(lr$Subfamily))){
    lr=lr %>% dplyr::rename(Family=Subfamily)
  }
  melted <- reshape2::melt(lr, id.vars = c("Pair", "Family"))
  melted$Family[is.na(melted$Family)] <- "NA"
  melted = melted %>% dplyr::arrange(Family, variable)
  melted <- melted %>% dplyr::mutate(row = group_indices_(melted,
                                                          .dots = c("Family", "Pair")))
  melted <- melted %>% dplyr::mutate(col = group_indices_(melted,
                                                          .dots = c("variable")))
  vars_x_axis <- c(melted %>% dplyr::arrange(col) %>% dplyr::select(variable) %>%
                     dplyr::distinct())$variable
  names_y_axis <- c(melted %>% dplyr::arrange(row) %>% dplyr::group_by(Pair) %>%
                      dplyr::distinct(Pair) %>% dplyr::select(Pair))$Pair
  plot <- ggplot2::ggplot(melted, ggplot2::aes(x = factor(col),
                                               y = factor(row), color = factor(Family), size = value)) +
    ggplot2::geom_point() +
    ggplot2::geom_text(ggplot2::aes(label = round(value), x = col + 0.4), alpha = 1, size = 3) +
    ggplot2::scale_size_area(max_size = 8) +
    ggplot2::scale_x_discrete(breaks = 1:length(vars_x_axis),
                              labels = vars_x_axis, position = "top") +
    ggplot2::scale_y_discrete(breaks = 1:length(names_y_axis), labels = names_y_axis) +
    ggplot2::scale_color_manual(values = family.col) +
    ggplot2::theme_bw() + ggplot2::labs(title = title) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                                                                        panel.border = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(),
                                                                        panel.grid.minor.x = ggplot2::element_blank(),
                                                                        axis.text.x = ggplot2::element_text(angle = 0), axis.ticks.x = ggplot2::element_blank(),
                                                                        axis.ticks.y = ggplot2::element_blank())
  return(plot)
}
