% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fmu_to_val.R
\name{eval_mom_hest}
\alias{eval_mom_hest}
\title{Compute numerical conditional moments for Heston SV}
\usage{
eval_mom_hest(fmu, par)
}
\arguments{
\item{fmu}{moment formula of the centralized conditional return}

\item{par}{list of parameters}
}
\value{
scalar moment
}
\description{
Compute numerical conditional central moments from their formulas for
the Heston SV model.
}
\examples{
fmu.hest = fmu.hest
v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61; rho = -0.7
r = 0.0319; tau = 1
par = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
moms8 = rep(0, 8) # centralized variable whose first moment = 0
for (i in 2:8) {moms8[i] = eval_mom_hest(fmu.hest[[i]], par)}
}
