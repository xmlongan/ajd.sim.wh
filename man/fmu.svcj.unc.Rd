% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fmu.svcj.unc}
\alias{fmu.svcj.unc}
\title{SVCJ Model Unconditional Central Moment Formulas}
\format{
\subsection{\code{fmu.svcj.unc}}{

A list with 8 elements:
\describe{
\item{\code{fmu1}}{The first central moment, equals 0}
\item{\code{fmu2}}{The second central moment}
\item{...}{The ... central moment}
\item{\code{fmu8}}{The eighth central moment}
}
Each of the formula consists of a list of 14 vectors:
\describe{
\item{\code{e^{kt}}}{Power of this term}
\item{...}{Power of this term}
\item{\code{sigma_s}}{Power of this term}
\item{\code{num}}{Numerator for this row coefficient}
\item{\code{den}}{Denominator for this row coefficient}
}
}
}
\usage{
fmu.svcj.unc
}
\description{
The first eight unconditional central moment formulas for the SVCJ model.
}
\keyword{datasets}
