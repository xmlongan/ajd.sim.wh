#' Simulate Partially Centralized Conditional Returns of SVCJ
#'
#' @description
#' Simulate the partially centralized conditional return samples of the SVCJ
#' model by the superposition of
#' - approximated Pearson Distribution which matches the first six|eight|ten
#'   conditional (partially) central moments of the unknown true SVVJ return
#'   distribution, see the Details. And
#' - exact contemporaneous jump in the price (not centralized)
#'
#' @param n number of return samples to simulate
#' @param moms vector of the first six|eight|ten conditional moments of the
#' centralized return.
#'
#' @returns vector of the partially centralized conditional return samples
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38; true_price = 6.8619
#'
#' n = 100
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = 1
#'
#' par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h)
#' moms = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms[i] = eval_mom_svcj(fmu.svcj[[i]], par)}
#'
#' Y = r_svcj(n, moms)
r_svcj <- function(n, moms) {
  # from moments to coefficients: a, c0, c1, c2, {c3 | c4 | c4, c5}
  coefs = mom_to_coef(moms)
  # Partial Fraction Decomposition
  int_p = int_p8; bisect = bisect8
  pfd = PFDecomp4(coefs)
  #
  # Determine support of the distribution
  #
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4]) # mean, var, skew, kurt
  sd = sqrt(stdmoment[2])
  skew = stdmoment[3]
  #
  if (skew < 0) {        # left-tailed
    # x = seq(mode - 7 * sd, mode + 3 * sd, length.out = N)
    lb = -7*sd; ub = 3*sd
  } else if (skew > 0) { # right-tailed
    # x = seq(mode - 3 * sd, mode + 7 * sd, length.out = N)
    lb = -7*sd; ub = 3*sd
  } else {               # symmetric
    # x = seq(mode - 5 * sd, mode + 5 * sd, length.out = N)
    lb = -5*sd; ub = 5*sd
  }
  lbub = adjust_lb_ub(lb, ub, pfd)
  lb = lbub[1]; ub = lbub[2]
  #
  logC_type = comp_logC(lb, ub, pfd, rel_err = 1e-7)
  if (logC_type$type == -2) {
    cat("resort to PearsonDS::rpearson() C either too big or small.\n")
    return(PearsonDS::rpearson(n, moments=stdmom(moms[1:4])))
  }
  logC = logC_type$logC
  #
  Us = stats::runif(n)
  # x0s = stats::qnorm(Us, mean = 0, sd = sd)
  x0s = PearsonDS::qpearson(Us, moments = stdmoment)
  #
  Xs = rep(0, n)
  for (i in 1:n) {
    U = Us[i]; x0 = x0s[i]
    j = 0
    Newton_Fail = FALSE
    repeat {
      if (x0 < lb || x0 > ub || j > 10 || Newton_Fail) {
        # if (x0 < lb || x0 > ub) { type_s = "TYPE1--out of bound" }
        # else if (j > 10)        { type_s = "TYPE2--long-iteration" }
        # else                    { type_s = "TYPE3--Newton-Condition-Failed" }
        x = bisect(lb, ub, logC, U, pfd, LB = lb, rel_err=1e-7)
        # if (x <= -1 || x >= 1) {
        #   cat(sprintf("logC = %f, type = %d\n", logC, logC_type$type))
        #   temp = "%sbisect()\n (lb, ub)= (%f, %f), U = %f, x0 = %f, x = %f\n"
        #   cat(sprintf(temp, type_s, lb, ub, U, x0, x))
        #   cat("moms = "); cat(moms, sep=","); cat("\n")
        # }
        # cat(sprintf("resort to bisect(), U = %f\n", U))
        break
      }
      prob = int_p(lb, x0, pfd, logC, rel_err = 1e-7)
      f = prob - U
      df = exp(log_dpearson8(x0, pfd) - logC) # df = dpearson(x0, pfd)/C
      f_df = f / df
      # ddf = -df * (coefs[1] + x0)/sum(coefs[2:6] * (x0^(0:4)))
      ddf_df = - (coefs[1] + x0)/sum(coefs[2:6] * (x0^(0:4)))
      delta = 1 - 2 * (f_df) * ddf_df # delta = 1 - 2 * (f/df) * (ddf/df)
      #
      # if (prob < 0 || prob > 1) {
      #   cat(sprintf("logC = %f, type = %d\n", logC, logC_type$type))
      #   temp = " f = %-e - %.7f, f/df = %-e, ddf/df = %-e, delta = %-e\n"
      #   cat(sprintf(temp, prob, U, f_df, ddf_df, delta))
      # }
      #
      if (delta >= 0 && ddf_df != 0) {
        x = x0 - (1 - sqrt(delta)) / ddf_df
      } else {
        Newton_Fail = TRUE; next
      }
      j = j + 1
      if (abs(x - x0) < 1e-5) { break }
      x0 = x
    }
    Xs[i] = x
  }
  return(Xs)
}
