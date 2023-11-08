#'
#' Filter the lr matrix (output from icellnet.score function) to select only the most contributing / differing interactions of the communications scores.
#'
#' @description Filter the lr matrix (output from icellnet.score function) to select only the most contributing / differing interactions of the communications scores : contributes with a score superior to a threshold and/or the top n interactions contributing to the score. Return a matrix, with interactions in rows and communication scores of communicating cell pairs in columns.
#' @param lr Matrix of individual communication scores
#' @param thresh Value set as a threshold to display the L/R pairs contributing with a score superior to this threshold
#' @param topn Value set as top n interactions to display
#' @param sort.by character, "sum" or "var", "sum" as default. In combination with topn parameter, allows to sort L/R pairs according to the most contributing ("sum"), or the most different among conditions ("var")
#
#' @export
#' @examples
#' \dontrun{LR.selection(lr = lr1, topn = 10, sort.by = "sum"))}
#'
LR.selection <- function (lr = lr, thresh = 0, topn = NULL, sort.by = "sum") {
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
    }
  }
  if (sort.by == "sum") {
    lr = lr %>% dplyr::mutate(sum = rowSums(dplyr::across(where(is.numeric)))) %>%
      dplyr::arrange(dplyr::desc(sum)) %>% dplyr::top_n(topn,
                                                        sum)
    lr = lr %>% dplyr::select(-sum)
  }
  else if (sort.by == "var") {
    lr$variance = apply(dplyr::select_if(lr, is.numeric),
                        1, var, na.rm = TRUE)
    lr = lr %>% dplyr::arrange(dplyr::desc(variance)) %>%
      dplyr::top_n(topn, variance)
    lr = lr %>% dplyr::select(-variance)
  } else stop("sort.by argument should be fixed on var or sum")

  rownames(lr)=lr$Pair
  lr =lr %>% dplyr::select(-Pair)

  return(lr)

}
