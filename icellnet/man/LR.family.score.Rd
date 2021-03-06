% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LR.family.score.R
\name{LR.family.score}
\alias{LR.family.score}
\title{Distribution of class of molecules contribution to the communication scores}
\usage{
LR.family.score(
  lr = lr,
  my.family = my.family,
  db.couple = as.data.frame(db.couple),
  plot = F,
  title = NULL,
  ymax = NULL,
  family.col = NULL
)
}
\arguments{
\item{lr}{Matrix of individual communication scores (sum on each column should be equal to the global communication score)}

\item{my.family}{Vector selecting family or subfamilies of molecules for the analysis}

\item{db.couple}{Output of the name.lr.couple() function. name.lr.couple(db, "Family") is set as a default}

\item{plot}{Logical, by default plot=F. If plot=T, function returns a ggplot object}

\item{title}{Title of the barplot (character)}

\item{ymax}{Value as limit of the yaxis. If ymax=NULL, ymax takes the maximum value of the table}

\item{family.col}{Color vector for the family of molecules c("family1"= "color1", "family2"="color2").}
}
\description{
Computes the contribution of families or subfamilies of molecules (included in my.selection) to the global communication scores.
If plot=FALSE, the function returns a table. If plot=T, returns a ggplot object, barplot representing contribution of families of molecules to the communication scores.
}
\examples{
\dontrun{LR.family.score(lr=lr, my.family=my.family, db.couple=db.name.couple, plot=F)}

}
