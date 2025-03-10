#' Simulate Partially Centralized Conditional Returns of SVJ
#'
#' @description
#' Simulate partially centralized conditional return samples of the SVJ model by
#' the superposition of
#' - approximated Pearson Distribution which matches the first six|eight|ten
#'   conditional moments of the centralized return of the Heston SV model, and
#' - exact jump distribution (not centralized).
#'
#' @param n number of return samples to simulate
#' @param tau time difference or maturity term
#' @param lmbd arrival rate \eqn{\lambda} of jumps in the price process
#' @param mu_b parameter \eqn{\bar{\mu}}, related to jump mean through
#'  \eqn{\mu_s = \log(1+\bar{\mu}) - \sigma_s^2/2}
#' @param sigma_s standard deviation of the normally distributed jumps
#' @param moms vector of the first six|eight|ten conditional moments of
#'  the centralized return of the Heston SV model
#'
#' @details
#' The SVJ model is described by the following SDEs:
#' \deqn{
#'   d\log s(t) = (\mu - v(t)/2) dt + \sqrt{v(t)}dw^s(t) + dz(t),
#' }
#' \deqn{
#'   dv(t) = k(\theta - v(t))dt + \sigma_v\sqrt{v(t)}dw^v(t)
#' }
#' where \eqn{z(t)} is a Compound Poisson process with constant arrival
#' rate \eqn{\lambda}. The (conditional, with \eqn{v_0} given) return over
#' period \eqn{[0,t]} is defined as \eqn{y_t \equiv \log s(t) -\log s(0)}
#' which is decomposed as
#' \deqn{
#'   y_t = y_{hest,t} + I\!Z_t,
#' }
#' where \eqn{y_{hest,t}} denotes the return of the Heston model and
#' \eqn{I\!Z_t\equiv \int_0^tdz(u)}. And return of the Heston model is
#' decomposed as
#' \deqn{
#'   y_{hest,t} = (\mu - \theta/2)t - \beta_t (v_0 -\theta) +
#'   \overline{y}_{hest,t}|v_0,
#' }
#' See [ajd.sim.wh::rpearson()] for
#' generating \eqn{\overline{y}_{hest,t}|v_0} samples from the Heston
#' SV model.
#' The function [ajd.sim.wh::r_svj()] generates the
#' \eqn{\overline{y}_{hest,t}|v_0 + I\!Z_t} samples.
#'
#' @seealso [ajd.sim.wh::rpearson()]
#' @return vector of \eqn{\overline{y}_{hest,t}|v_0 + I\!Z_t} samples.
#' @export
#'
#' @examples
#' n = 100
#' # moms6 = c(0.0231, 0.0191, -0.0021, 0.0020, -0.0010, 0.0008)
#' v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27; rho = -0.79
#' lmbd = 0.11; mu_b = -0.12; sigma_s = 0.15; tau = 5
#'
#' par = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_hest(fmu.hest[[i]], par)}
#'
#' Y = r_svj(n, tau, lmbd, mu_b, sigma_s, moms8)
r_svj <- function(n, tau, lmbd, mu_b, sigma_s, moms) {
  # Now moms is computed with mu = r-lmbd*mu_b
  mu_s = log(1 + mu_b) - sigma_s^2 / 2
  #
  Y_c = rpearson(n, moms)
  Y_J = rep(0, n)
  #
  numJ = stats::rpois(n, lmbd * tau)
  for (i in 1:n) {
    Y_J[i] = sum(stats::rnorm(numJ[i], mean = mu_s, sd = sigma_s))
  }
  Xs = Y_c + Y_J
  return(Xs)
}
