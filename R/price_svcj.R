#' Pricing the European Call Option Under the SVCJ Model
#'
#' @description
#' Pricing the European call option under the SVCJ model using the
#' Wu-Hu (2024) method for the simulation.
#'
#' @param N number of underling asset price samples to simulate.
#' @param S current underling asset price.
#' @param K striking price.
#' @param v0 current variance level.
#' @param tau maturity time, \eqn{T}.
#' @param r risk-less rate.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param lmbd parameter \eqn{\lambda}, arrival rate of jumps in the
#' underling asset price process and variance process.
#' @param mu_v parameter \eqn{\mu_v}, mean of the jumps
#' @param mu_b parameter \eqn{\bar{\mu}}, see Details.
#' @param rhoJ parameter \eqn{\rho_J}, correlation parameter between the
#' jumps in the price and variance.
#' @param sigma_s parameter \eqn{\sigma_s}, i.e., standard deviation of
#'  jumps in the price.
#' (exponential distribution) in the variance.
#' @param moms vector of the first six|eight|ten conditional moments of the
#' centralized return.
#' @param true_price theoretical true price of the option, calculable by
#' the Duffie et al. (2000) method.
#'
#' @details
#' The parameter \eqn{\bar{\mu}} relates to jump mean
#'  \eqn{\mu_s + \rho_J J^v} in the price through
#'  \deqn{\mu_s = \log((1+\bar{\mu})(1-\rho_J\mu_v)) - \sigma_s^2/2.}
#'  Note that \eqn{J^v} denotes jump size in the variance.
#'  Meanwhile the expected return per unit time of the SVCJ model becomes
#'  \deqn{r - \lambda\bar{\mu},} see Broadie and Kaya (2006) for details.
#'
#' @seealso [ajd.sim.wh::r_svcj()]
#' @returns vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38; true_price = 6.8619
#'
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = 1
#'
#' par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h)
#' moms = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms[i] = eval_mom_svcj(fmu.svcj[[i]], par)}
#'
#' n = 10000
#' price_svcj(n, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s,
#'            moms, true_price)
price_svcj <- function(N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ,
                       sigma_s, moms, true_price) {
  # start.time = Sys.time()
  ts = proc.time()
  Y = r_svcj(N, moms)
  #
  mu = r - lmbd * mu_b
  beta = (1 - exp(-k * tau)) / (2 * k)
  mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
  #
  Ymean = lmbd * mu_v * beta / k + (rhoJ - 1/(2*k)) * lmbd * mu_v * tau
  Ymean = Ymean + lmbd * mu_s * tau + (mu - theta/2)*tau - beta*(v0 - theta)
  #
  Y = Y + Ymean
  #
  cprice_MC = exp(-r * tau) * mean(pmax(S * exp(Y) - K, 0))
  # end.time = Sys.time()
  te = proc.time()
  tt = te - ts
  # time.taken = end.time - start.time
  # time.taken = difftime(end.time, start.time, units = "secs")
  error = cprice_MC - true_price
  #
  # return(c(error, as.numeric(time.taken)))
  return(c(error, tt[[3]]))
}
