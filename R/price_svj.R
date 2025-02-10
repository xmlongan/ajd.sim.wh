#' Pricing the European Call Option Under the SVJ Model
#'
#' @description
#' Pricing the European call option under the SVJ model using the
#' Wu-Hu (2024) method for the simulation.
#'
#' @param N number of underling asset price samples to simulate.
#' @param S current underling asset price.
#' @param K striking price.
#' @param v0 initial variance, \eqn{v_0}.
#' @param tau maturity time, \eqn{T}.
#' @param r risk-less rate.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param lmbd parameter \eqn{\lambda}, arrival rate of jumps in the
#'  underling asset price process.
#' @param mu_b parameter \eqn{\bar{\mu}}, related to the mean of the jumps,
#'   \eqn{\mu_s = \log (1+\bar{\mu}) - \sigma_s^2/2}, see Broadie-Kaya (2006).
#'   Meanwhile the expected return per unit time of the SVJ model becomes
#'   \eqn{r - \lambda\bar{\mu}}.
#' @param sigma_s parameter \eqn{\sigma_s}, standard deviation of the jumps,
#'   see Broadie-Kaya (2006).
#' @param moms vector of the first six|eight|ten conditional moments of
#'   the centralized return of the Heston SV model
#' @param true_price theoretical true price of the option, calculable by
#'  the Bates (1996) method.
#'
#' @seealso [ajd.sim.wh::r_svj()]
#' @return vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27
#' rho = -0.79; r = 0.0319; tau = 5; lmbd = 0.11; mu_b = -0.12
#' sigma_s = 0.15; true_price = 20.1642
#'
#' par_hest = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms = rep(0, 8) # centralized variable, whose first moment always = 0
#' for (i in 2:8) {moms[i] = eval_mom_hest(fmu.hest[[i]], par_hest)}
#'
#' N = 10000
#' price_svj(N, S, K, v0, tau, r, k, theta, lmbd, mu_b, sigma_s, moms,
#'           true_price)
price_svj <- function(N, S, K, v0, tau, r, k, theta, lmbd, mu_b, sigma_s,
                      moms, true_price) {
  # start.time = Sys.time()
  ts = proc.time()
  Y = r_svj(N, tau, lmbd, mu_b, sigma_s, moms)
  #
  mu = r - lmbd * mu_b
  beta = (1 - exp(-k * tau)) / (2 * k)
  # adding back the non-random part
  Y = Y + (mu - theta / 2) * tau - beta * (v0 - theta)
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
