#' Simulation of the Square-Root Jump diffusion
#'
#' @description
#' Generate samples from the Square-Root Jump Diffusion (SRJD).
#'
#' @param v0 the initial state.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma volatility of the SRJD.
#' @param lmbd jump intensity \eqn{\lambda}.
#' @param mu_v mean of the exponential jump distribution.
#' @param h terminal time.
#'
#' @return scalar
#' @export
#'
#' @examples
#' v0 = 0.007569
#' k = 3.46; theta = 0.008; sigma = 0.14; lmbd = 0.47; mu_v = 0.05; h = 1
#'
#' r_srjd(v0, k, theta, sigma, lmbd, mu_v, h)
r_srjd <- function(v0, k, theta, sigma, lmbd, mu_v, h) {
  n = rpois(1, lmbd * h)
  if (n == 0) {
    v = ajd.sim.bk::rv(v0, h, k, theta, sigma)  # only diffusion
  } else {
    jtime = sort(runif(n, 0, h))
    jsize = rexp(n, 1/mu_v)
    v = v0
    for (j in 1:n) {
      delta_t = ifelse(j==1, jtime[1], jtime[j] - jtime[j-1])
      # diffusion + jump
      v = ajd.sim.bk::rv(v, delta_t, k, theta, sigma) + jsize[j]
    }
    # the last diffusion
    v = ajd.sim.bk::rv(v, h - jtime[n], k, theta, sigma)
  }
  return(v)
}
