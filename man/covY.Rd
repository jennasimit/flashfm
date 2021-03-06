% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covY.R
\docType{data}
\name{covY}
\alias{covY}
\title{Covariance matrix for two simulated traits}
\format{An object of class \code{matrix} with 2 rows and 2 columns.}
\usage{
covY
}
\description{
Simulated data was generated using the mvrnorm function in R and genotype matrix X. The two traits have correlation 0.4 and 
each has two causal variants, of which one, rs61839660 is shared between traits; trait 1 has second 
causal variant rs62626317 and trait 2 has second causal variant rs11594656. The traits are from a sample
of 2000, with no missing measurements.
}
\details{
The matrix has row and column names coinciding with the trait names.
}
\keyword{datasets}
