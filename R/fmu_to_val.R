#' Compute numerical conditional moments for Heston SV
#'
#' @description
#' Compute numerical conditional central moments from their formulas for
#' the Heston SV model.
#'
#' @param fmu moment formula of the centralized conditional return
#' @param par list of parameters
#'
#' @return scalar moment
#' @export
#'
#' @examples
#' fmu.hest = fmu.hest
#' v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61; rho = -0.7
#' r = 0.0319; tau = 1
#' par = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_hest(fmu.hest[[i]], par)}
eval_mom_hest <- function(fmu, par) {
  v0 = par$v0
  k = par$k
  theta = par$theta
  sigma = par$sigma
  rho = par$rho
  h = par$h
  #
  col1 = exp(-fmu$`e^{-kt}` * k * h)
  col2 = h^fmu$t
  col3 = k^(-fmu$`k^{-}`)
  col4 = (v0-theta)^fmu$`v_0-theta`
  col5 = theta^fmu$theta
  col6 = sigma^fmu$sigma
  col7 = rho^fmu$rho
  coef = fmu$num / fmu$den
  #
  val = sum(col1 * col2 * col3 * col4 * col5 * col6 * col7 * coef)
  #
  return(val)
}
