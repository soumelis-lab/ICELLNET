
#' Plot of individual communication scores of contributing L/R pairs.
#'
#' @description Display the individual communication scores of ligand/receptor pairs that contributes with a score superior to a threshold and/or the top n interactions contributing to the score. Return a ggplot object.
#'
#' @param lr Matrix of individual communication scores
#' @param thresh Value set as a threshold to display the L/R pairs contributing with a score superior to this threshold
#' @param topn Value set as top n interactions to display
#' @param db.couple output of the name.lr.couple() function. name.lr.couple(db, "Family") is used as default
#' @param sort.by character, "sum" or "var", "sum" as default. In combination with topn parameter, allows to sort L/R pairs according to the most contributing ("sum"), or the most different among conditions ("var")
#' @param title Title of the heatmap
#' @param value_display Numerical value, used as a threshold value to add numerical values to the heatmap above this value.
#
#'
#' @export
#' @examples
#' \dontrun{LR.heatmap(lr = lr1, topn = 30, sort.by = "sum")
#' }
#'


LR.heatmap <- function (lr = lr, thresh = 0, topn = NULL, sort.by = "sum",
                         db.couple = NULL, title = title, value_display=NA) {
  if (is.null(db.couple)){
    db.couple=name.lr.couple(db, "Family")
  }
  interactions = rownames(lr)
  lr = as.data.frame(lr)
  lr$Pair = interactions
  lr = as.data.frame(lr[stats::complete.cases(lr), ])
  lr = lr %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. >
                                                             0)) %>% dplyr::filter_if(is.numeric, dplyr::any_vars(. >=
                                                                                                                    thresh))
  if (!is.null(topn)) {
    if (topn > dim(lr)[1]) {
      note(paste0("lr contains only ", dim(lr)[1], " after filtering interaction highest than theshold"))
      lr = dplyr::left_join(lr, as.data.frame(db.couple),
                            by = "Pair")
    }
    else {
      if (sort.by == "sum") {
        lr = lr %>% dplyr::mutate(sum = rowSums(dplyr::across(where(is.numeric)))) %>%
          dplyr::arrange(dplyr::desc(sum)) %>% dplyr::top_n(topn,
                                                            sum)
        lr = dplyr::left_join(lr, as.data.frame(db.couple),
                              by = "Pair") %>% dplyr::select(-sum)
      }
      else if (sort.by == "var") {
        lr$variance = apply(dplyr::select_if(lr, is.numeric),
                            1, var, na.rm = TRUE)
        lr = lr %>% dplyr::arrange(dplyr::desc(variance)) %>%
          dplyr::top_n(topn, variance)
        lr = dplyr::left_join(lr, as.data.frame(db.couple),
                              by = "Pair") %>% dplyr::select(-variance)
      }
      else stop("sort.by argument should be fixed on var or sum")
    }
  }
  else {
    lr = dplyr::left_join(lr, as.data.frame(db.couple),
                          by = "Pair")
  }
  if (!(is.null(lr$Subfamily))) {
    lr = lr %>% dplyr::rename(Family = Subfamily)
  }
  melted <- reshape2::melt(lr, id.vars = c("Pair", "Family"))
  melted$Family[is.na(melted$Family)] <- "NA"
  melted = melted %>% dplyr::arrange(Family, variable)

  if (is.na(value_display)){
    plot <- ggplot(melted, ggplot2::aes(x = variable,
                                        y = Pair, fill=value)) +  ggplot2::geom_tile(aes(fill = value)) +
      ggplot2::scale_fill_gradient(low="white", high="red")+
      ggplot2::theme(axis.text.x=element_text(angle=45, hjust = 1))  +
      ggplot2::labs(title = title)
  }else{
    plot <- ggplot(melted, ggplot2::aes(x = variable,
                                        y = Pair, fill=value)) +  ggplot2::geom_tile(aes(fill = value)) +
      ggplot2::scale_fill_gradient(low="white", high="red")+
      ggplot2::theme(axis.text.x=element_text(angle=45, hjust = 1)) +
      ggplot2::geom_text( data=melted %>% filter(value>value_display), ggplot2::aes(label = round(value)), size = 2) +
      ggplot2::labs(title = title)
  }
  return(plot)
}

