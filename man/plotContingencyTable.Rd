% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotContingencyTable.R
\name{plotContingencyTable}
\alias{plotContingencyTable}
\title{plots the contingency table of two clustering results}
\usage{
plotContingencyTable(
  cluster_labels_1 = NULL,
  cluster_labels_2 = NULL,
  automateConsensus = T,
  minClustSize = 10,
  filename = "Contingency_Table.pdf"
)
}
\arguments{
\item{cluster_labels_1}{First vector of cluster labels for each cell}

\item{cluster_labels_2}{Second vector of cluster labels for each cell}

\item{automateConsensus}{Boolean indicating whether automated consensus should be returned}

\item{minClustSize}{Minimum number of cells in a cluster. Default = 10.}

\item{filename}{name of contingency table file}
}
\value{
consensus cluster labels vector
}
\description{
plots the contingency table of two clustering results
}
\author{
ranjanb
}
