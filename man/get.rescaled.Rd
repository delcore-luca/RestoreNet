% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get.rescaled}
\alias{get.rescaled}
\title{Rescaling a clonal tracking dataset}
\usage{
get.rescaled(Y)
}
\arguments{
\item{Y}{A 3-dimensional array whose dimensions are the time, the cell type and the clone respectively.}
}
\value{
A rescaled clonal tracking dataset.
}
\description{
Rescales a clonal tracking dataset based on the sequencing depth.
}
\details{
This function rescales a clonal tracking dataset Y according to the formula
\deqn{Y_{ijk} \leftarrow Y_{ijk} \cdot \frac{min_{ij}\sum_cY_{ijc}}{\sum_cY_{ijc}}}{Y_{ijk} <- Y_{ijk}*(min_{ij}sum_cY_{ijc})/(sum_cY_{ijc})}
}
