#' Pricing the European Call Option Under the Heston Model
#'
#' @description
#' Pricing the European call option under the Heston model using the
#' Wu-Hu (2024) method for the simulation.
#'
#' @param N number of underling asset price samples to simulate.
#' @param S current underling asset price.
#' @param K striking price.
#' @param v0 the initial variance
#' @param tau maturity time, \eqn{T}.
#' @param r risk-less rate.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param moms vector of the first six|eight|ten conditional moments of the
#' centralized return.
#' @param true_price theoretical price of the option, calculable from
#' the Heston (1993) method.
#'
#' @seealso [ajd.sim.wh::r_hest()]
#' @return vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' # setting 1
#' S = 100; K = 100; v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61
#' rho = -0.7; r = 0.0319; tau = 1; true_price = 6.8061
#'
#' par_hest = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms = rep(0, 8) # centralized variable, whose first moment always = 0
#' for (i in 2:8) {moms[i] = eval_mom_hest(fmu.hest[[i]], par_hest)}
#'
#' price_hest(10000, S, K, v0, tau, r, k, theta, moms, true_price)
#'
#' # setting 2
#' S = 100; K = 100; v0 = 0.09; k = 2.00; theta = 0.09; sigma = 1.00
#' rho = -0.3; r = 0.05; tau = 5
#' true_price = 34.9998
#'
#' par_hest = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms = rep(0, 8) # centralized variable, whose first moment always = 0
#' for (i in 2:8) {moms[i] = eval_mom_hest(fmu.hest[[i]], par_hest)}
#'
#' price_hest(10000, S, K, v0, tau, r, k, theta, moms, true_price)
#'
price_hest <- function(N, S, K, v0, tau, r, k, theta, moms, true_price) {
  # start.time = Sys.time()
  ts = proc.time()
  #
  beta = (1 - exp(-k * tau)) / (2 * k)
  # adding back the non-random part
  # Y = (r - theta / 2) * tau - beta * (v0 - theta) + r_hest(N, moms)
  Y = (r - theta / 2) * tau - beta * (v0 - theta) + r_hest2(N, moms)
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
