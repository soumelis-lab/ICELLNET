#'
#' Plot of individual communication scores of contributing L/R pairs.
#'
#' @description Display the individual communication scores of ligand/receptor pairs that contributes with a score superior to a threshold (thresh = 10 as a default)
#'
#' @param lr Matrix of individual communication scores
#' @param PC Vector of selected peripheral cell types
#' @param thresh Value set as a threshold to display the L/R pairs contributing with a score superior to this threshold
#' @param type "raw" as default, can be set to "percentage" if you wish to consider the threshold as a defined percentage of the global communication score
#' @param db.name.couple output of the name.lr.couple() function. name.lr.couple(db, "Family") is set as a default
#' @param title Title of the balloon plot
#
#'
#' @export
#' @examples
#' \dontrun{lr <- icellnet.score.results[[2]]
#' LR.balloon.plot(lr = lr, PC =PC , thresh = 20,  db.name.couple=name.lr.couple(db, "Family"))
#' }
#'
LR.balloon.plot = function (lr = lr,
                            PC = PC,
                            thresh = thresh,
                            type = c("raw" , "percentage"),
                            db.name.couple = db.name.couple,
                            title = title) {
  lr=lr[complete.cases(lr),]
  lr=lr[is.finite(log2(rowSums(lr))),]
  rank0 = as.data.table(x = character(0))
  i = 2
  if (length(PC) == 1) {
    ranked.lr.cell = as.data.table(sort(lr[, cell], decreasing = TRUE), keep.rownames =
                                     TRUE)
    if (type == "percentage") {
      ranked.lr.cell$V2 = round(ranked.lr.cell$V2 / sum(lr[, cell]) * 100, 1)#as.numeric(score)*100,1)
    }
    if (is.null(which (as.numeric(sum(
      ranked.lr.cell$V2
    )) > as.numeric(thresh)))) {
      i = i + 1
    } else
      rank = as.data.table (ranked.lr.cell[which (as.numeric(sum(ranked.lr.cell$V2)) > as.numeric(thresh))] , keep.rownames =
                              TRUE)
    rank0 = rank
  } else
    for (cell in PC) {
      ranked.lr.cell = as.data.table(sort(lr[, cell], decreasing = TRUE), keep.rownames =
                                       TRUE)
      if (type == "percentage") {
        ranked.lr.cell$V2 = round(ranked.lr.cell$V2 / sum(lr[, cell]) * 100, 1)#as.numeric(score)*100,1)
      }
      if (is.null(which (as.numeric(sum(
        ranked.lr.cell$V2
      )) > as.numeric(thresh)))) {
        i = i + 1
      } else
        rank = as.data.table (ranked.lr.cell[which (as.numeric(ranked.lr.cell$V2) > as.numeric(thresh))] , keep.rownames =
                                TRUE)
      rank0 <- merge(rank0, rank, by = "V1", all = TRUE)
      colnames(rank0)[i] <- cell
      i = i + 1
    }
  vars <- PC
  rank0$family = rep(NA, length(rank0$V1))
  for (i in 1:length(rank0$V1)) {
    rank0$family[i] = as.character(db.name.couple[which(db.name.couple[, 1] ==
                                                          rank0$V1[i]), 2])
  }
  melted <- melt(rank0, id.vars = c("V1", "family"))
  melted$family[is.na(melted$family)] <- "NA"
  melted = melted %>% arrange(family, variable)
  melted <-
    melted %>% mutate(row = group_indices_(melted, .dots = c('family', 'V1')))
  melted <-
    melted %>% mutate(col = group_indices_(melted, .dots = c('variable')))
  vars_x_axis <-
    c(melted %>% arrange(col) %>% select(variable) %>% distinct())$variable
  names_y_axis <-
    c(melted %>% arrange(row) %>%  group_by(V1) %>% distinct(V1) %>%  select(V1))$V1
  plot <-
    ggplot(melted,
           aes(
             x = factor(col),
             y = factor(row),
             color = factor(family),
             size = value #,
             #alpha = value
           )) +
    geom_point() +    # plot as points
    geom_text(aes(label = round(value), x = col + 0.4),
              alpha = 1.0,
              size = 3) +   # display the value next to the "balloons"
    #scale_alpha_continuous(range = c(0.4, 0.8)) +
    scale_size_area(max_size = 8) +
    scale_x_discrete(
      breaks = 1:length(vars_x_axis),
      labels = vars_x_axis,
      position = 'top'
    ) +   # set the labels on the X axis
    scale_y_discrete(breaks = 1:length(names_y_axis), labels = names_y_axis) +                 # set the labels on the Y axis
    scale_color_manual(values = family.col) +
    theme_classic() +
    labs(title = title) +
    theme(
      axis.line = element_blank(),
      # disable axis lines
      axis.title = element_blank(),
      # disable axis titles
      panel.border = element_blank(),
      # disable panel border
      panel.grid.major.x = element_blank(),
      # disable lines in grid on X-axis
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_text(angle = 90),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank()
    )

  return(plot)
}
