% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group.multi.R
\name{best.snps}
\alias{best.snps}
\title{Best SNPs}
\usage{
best.snps(d, mppi.thr = pp.thr, pp.thr = 0.1)
}
\arguments{
\item{d}{snpmod object or a list of snpmod objects}

\item{mppi.thr}{MPPI threshold, SNPs with MPPI>mppi.thr will be shown}

\item{pp.thr}{deprecated, alias for mppi.thr}
}
\value{
subset of \code{snps(d)} data.frame including best SNPs
}
\description{
Display the SNPs with greatest marginal posterior probability of inclusion (MPPI)
}
\author{
Chris Wallace
}
