% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reclusterDEConsensus.R
\name{reclusterDEConsensus}
\alias{reclusterDEConsensus}
\title{Recluster consensus clusters using DE gene analysis}
\usage{
reclusterDEConsensus(
  dataMatrix,
  consensusClusterLabels,
  method = "Wilcoxon",
  meanScalingFactor = 5,
  qValThrs,
  fcThrs,
  deepSplitValues = 1:4,
  minClusterSize = 10,
  filename = "de_gene_object.rds",
  plotName = "DE_Heatmap"
)
}
\arguments{
\item{dataMatrix}{the log-transformed and normalized scRNAseq genes x cells matrix}

\item{consensusClusterLabels}{consensus cell type labels for each cell}

\item{method}{Method used for DE gene statistical test}

\item{meanScalingFactor}{scale of the mean gene expression across the gene expression matrix to set a minimum threshold of average cluster expression for a gene to be considered a DE gene.}

\item{qValThrs}{maximum q-value threshold for statistical test}

\item{fcThrs}{minimum fold-change threshold for DE gene criterion (natural log following Seurat convention)}

\item{deepSplitValues}{vector of WGCNA tree cutting deepsplit parameters}

\item{minClusterSize}{specifies the type of minimum cluster size}

\item{filename}{name of DE gene object file}

\item{plotName}{name of DE Heatmap file}
}
\value{
list containing vector of DE genes, clustering tree and dynamic color list.
}
\description{
Recluster consensus clusters using DE gene analysis
}
\author{
ranjanb, schmidtf
}
