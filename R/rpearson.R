#' Simulate Samples from Pearson Distribution
#'
#' @description
#' Simulate samples from the approximated Pearson Distribution which matches
#' the first eight moments of the unknown true return distribution.
#'
#' @param n number of return samples to simulate
#' @param moms vector of the first eight conditional central moments of
#' the return
#'
#'
#' @returns vector of centralized conditional return samples,
#'   \eqn{\overline{y}_t|v_0}.
#' @export
#'
#' @examples
#' v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61; rho = -0.7
#' r = 0.0319; tau = 1
#'
#' par = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#'
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_hest(fmu.hest[[i]], par)}
#'
#' n = 100
#' Y = rpearson(n, moms8)
rpearson <- function(n, moms) {
  # from moments to coefficients: a, c0, c1, c2, {c3 | c4 | c4, c5}
  coefs = mom_to_coef(moms)
  # Partial Fraction Decomposition
  pfd = PFDecomp4(coefs)
  pfd_tp = 4
  dpearson = dpearson8
  #
  # Determine support of the distribution
  #
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4]) # mean, var, skew, kurt
  sd = sqrt(stdmoment[2])
  lbub = PearsonDS::qpearson(c(0.000001, 0.999999), moments = stdmoment)
  lbub = c(lbub[1] - sd, lbub[2] + sd)  # enlarge the support
  # make sure no roots lying within
  lbub = adjust_lb_ub(lbub[1], lbub[2], pfd)
  #
  N = 10000
  x = seq(lbub[1], lbub[2], length.out = N)
  # Discretize and evaluate
  #
  dx = dpearson(x, pfd) # un-normalized density
  # trapezoidal rule
  h = x[2] - x[1]
  # cumsum the middle ones except the first and the last
  cum = cumsum(dx[2:(N-1)] * h)  # may too big number inside!
  #
  px = rep(0, N)        # un-normalized cumulative probability
  px[2] = (dx[1] + dx[2]) * h / 2
  px[3:N] = (dx[1] + dx[3:N]) * h / 2 + cum
  C = px[N] # constant
  #
  Us = stats::runif(n)
  Xs = rep(0, n)
  for (i in 1:n) {
    Xs[i] = bound_refine(Us[i] * C, px, dx, x, pfd, pfd_tp)
  }
  return(Xs)
}
