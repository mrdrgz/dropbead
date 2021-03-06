% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/DGEMatrix.R
\name{removeLowQualityCells}
\alias{removeLowQualityCells}
\title{Filter out low quality cells}
\usage{
removeLowQualityCells(object, min.genes = 1000, ...)
}
\arguments{
\item{object}{A \code{data.frame} representing DGE.}

\item{n}{The minimum number of genes required.}
}
\value{
A \code{data.frame} without the low quality cells.
}
\description{
Remove cells expressing less than a minimum of genes (default value is 2000).
}

