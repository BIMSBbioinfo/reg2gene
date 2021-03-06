% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataIntegRFunctions.R
\name{normalizeScores}
\alias{normalizeScores}
\title{Function that runs normalization for the quantified activity/expression level}
\usage{
normalizeScores(scoresDF, normalizationMet)
}
\arguments{
\item{scoresDF}{a dataframe of quatified scores per genome window and bigwig 
file}

\item{normalizationMet}{charater. which normalization procedure is wanted.}
}
\value{
a normalized dataframe of quatified scores per genome window and 
bigwig files.
}
\description{
Function that runs normalization for the quantified activity/expression level
}
\author{
Inga Patarcic
}
\keyword{internal}
