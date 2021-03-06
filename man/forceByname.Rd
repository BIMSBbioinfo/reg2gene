% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BenchmarkGInteractions.R
\name{forceByname}
\alias{forceByname}
\title{Benchmarks help function - forceByName}
\usage{
forceByname(benchInteractions, interactions)
}
\description{
Forces benchmark input object to use gene coordinates of interactions object 
if
gene names overlap. It is crucial that benchInteractions 
& interactions both contain
meta-data named "name" -> where gene names are stored, 
and that anchor2 in both objects corresponds to gene object.
}
\details{
Read description for \code{benchmarkInteractions}
}
\author{
Inga Patarcic
}
\keyword{internal}
