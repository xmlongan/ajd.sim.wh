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


#' Compute numerical conditional moments for SVCJ
#'
#' @description
#' Compute numerical conditional central moments from their formulas for
#' the SVCJ model.
#'
#' @param fmu moment formula of the centralized conditional return
#' @param par list of parameters
#'
#' @return scalar moment
#' @export
#'
#' @examples
#' # fmu.svcj = fmu.svcj
#' v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38
#'
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = 1
#'
#' par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h)
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_svcj(fmu.svcj[[i]], par)}
eval_mom_svcj <- function(fmu, par) {
  v0 = par$v0
  mu = par$mu
  k = par$k
  theta = par$theta
  sigma = par$sigma
  rho = par$rho
  lmbd = par$lmbd
  mu_v = par$mu_v
  rhoJ = par$rhoJ
  mu_s = par$mu_s
  sigma_s = par$sigma_s
  h = par$h
  #
  col1 = exp(fmu$`e^{kt}` * k * h)
  col2 = h^fmu$t
  col3 = k^(-fmu$`k^{-}`)
  col4 = mu^fmu$`mu`
  col5 = (v0-theta)^fmu$`v0-theta`
  col6 = theta^fmu$`theta`
  col7 = sigma^fmu$`sigma`
  col8 = rho^fmu$`rho`
  col9 = lmbd^fmu$`lmbd`
  col10 = mu_v^fmu$`mu_v`
  col11 = rhoJ^fmu$`rhoJ`
  col12 = mu_s^fmu$`mu_s`
  col13 = sigma_s^fmu$`sigma_s`
  #
  coef = fmu$num / fmu$den
  #
  val = sum(col1 * col2 * col3 * col4 * col5 * col6 *
              col7 * col8 * col9 * col10 * col11 * col12 * col13 * coef)
  #
  return(val)
}


#' Compute numerical unconditional moments for SVCJ
#'
#' @description
#' Compute numerical unconditional central moments from their formulas for
#' the SVCJ model.
#'
#' @param fmu moment formula of the centralized return
#' @param par list of parameters
#'
#' @return scalar moment
#' @export
#'
#' @examples
#' # fmu.svcj.unc = fmu.svcj.unc
#' v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38
#'
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = 1
#'
#' par = list(mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h)
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_svcj_unc(fmu.svcj.unc[[i]], par)}
eval_mom_svcj_unc <- function(fmu, par) {
  # v0 = par$v0
  # mu = par$mu
  k = par$k
  theta = par$theta
  sigma = par$sigma
  rho = par$rho
  lmbd = par$lmbd
  mu_v = par$mu_v
  rhoJ = par$rhoJ
  mu_s = par$mu_s
  sigma_s = par$sigma_s
  h = par$h
  #
  mean_v = theta + lmbd * mu_v / k  # E[v] = theta + lmbd * mu_v / k
  #
  col1 = exp(-fmu$`e^{-kt}` * k * h)
  col2 = h^fmu$t
  col3 = k^(-fmu$`k^{-}`)
  col4 = mean_v^fmu$`E[v]`
  col5 = theta^fmu$`theta`
  col6 = sigma^fmu$`sigma`
  col7 = rho^fmu$`rho`
  col8 = lmbd^fmu$`lmbd`
  col9 = mu_v^fmu$`mu_v`
  col10 = rhoJ^fmu$`rhoJ`
  col11 = mu_s^fmu$`mu_s`
  col12 = sigma_s^fmu$`sigma_s`
  #
  coef = fmu$num / fmu$den
  #
  val = sum(col1 * col2 * col3 * col4 * col5 * col6 *
              col7 * col8 * col9 * col10 * col11 * col12 * coef)
  #
  return(val)
}
