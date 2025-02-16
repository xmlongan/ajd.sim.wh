#' Simulate Centralized Conditional Returns of Heston SV
#'
#' @description
#' Simulate centralized conditional return samples of the Heston SV model from
#' approximated Pearson Distribution which matches the first six|eight|ten
#' moments of the unknown true return distribution.
#'
#' @param n number of return samples to simulate
#' @param moms vector of the first six|eight|ten conditional moments of the
#'  centralized return
#'
#' @details
#' The Heston SV model is described by the following SDEs:
#' \deqn{
#'   d\log s(t) = (\mu - v(t)/2)dt + \sqrt{v(t)}dw^s(t),
#' }
#' \deqn{
#'   dv(t) = k(\theta - v(t))dt + \sigma_v\sqrt{v(t)}dw^v(t),
#' }
#' where the correlation between \eqn{w^s(t)} and \eqn{w^v(t)} is equal to
#' \eqn{\rho}, i.e., \eqn{cov(w^s(t), w^v(t)) = \rho t}.
#' The (conditional, with \eqn{v_0} given) return over period \eqn{[0,t]}
#' is defined as \eqn{y_t \equiv \log s(t) - \log s(0)}
#' which is decomposed as
#' \deqn{
#'   y_t = (\mu - \theta/2)t - \beta_t (v_0 -\theta) + \overline{y}_t|v_0,
#' }
#' where \eqn{\beta_t = (1-e^{-kt})/(2k)} and the centralized conditional
#' return
#' \deqn{
#'   \overline{y}_t|v_0 = \frac{\sigma_v}{2k}e^{-kt}I\!E_t +
#'   \left(\rho - \frac{\sigma_v}{2k}\right)I_t +
#'   \sqrt{1-\rho^2}I_t^{*}.
#' }
#' Given the initial variance \eqn{v_0}, the function generate samples of
#' the centralized conditional return \eqn{\overline{y}_t|v_0}.
#'
#' @returns vector of centralized conditional return samples,
#'   \eqn{\overline{y}_t|v_0}.
#' @export
#'
#' @examples
#' v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61; rho = -0.7
#' r = 0.0319; tau = 1
#' par = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
#' moms8 = rep(0, 8) # centralized variable whose first moment = 0
#' for (i in 2:8) {moms8[i] = eval_mom_hest(fmu.hest[[i]], par)}
#'
#' n = 100
#' Y = r_hest(n, moms8)
#'
#' # n = 2560000
#' n = 10000
#' time.start = proc.time()
#' Y = r_hest2(n, moms8)
#' time.end   = proc.time()
#' time.taken = time.end - time.start
#' time.taken
#'
#' # moms6 = c(0.0231, 0.0191, -0.0021, 0.0020, -0.0010, 0.0008)
r_hest <- function(n, moms) {
  # from moments to coefficients: a, c0, c1, c2, {c3 | c4 | c4, c5}
  coefs = mom_to_coef(moms)
  # Partial Fraction Decomposition
  if (length(moms) == 6) {
    pfd_tp = 3; dpearson = dpearson6
    pfd = PFDecomp3(coefs); x1 = pfd$x1; x2 = pfd$x2
    if (Im(x2) == 0) stop("expect 1 real root and 1 pair of conjugates!")
  } else if (length(moms) == 8) {
    pfd_tp = 4; dpearson = dpearson8
    pfd = PFDecomp4(coefs); x1 = pfd$x1; x2 = pfd$x2
    if (x1 == x2) stop("The two pairs of roots are not distinct!")
    dpearson = dpearson8
  } else if (length(moms) == 10) {
    pfd_tp = 5; dpearson = dpearson10
    pfd = PFDecomp5(coefs); x2 = pfd$x2; x3 = pfd$x3
    if (x2 == x3) stop("The two pairs of roots are not distinct!")
  } else stop("moms must be a vector of the first 6|8|10 moments")
  #
  # Determine support of the distribution
  #
  mode = -coefs[1] # -a: mode of the distribution
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4]) # mean, var, skew, kurt
  sd = sqrt(stdmoment[2])
  skew = stdmoment[3]
  #
  N = 10000
  if (skew < 0) {        # left-tailed
    x = seq(mode - 7 * sd, mode + 4.5 * sd, length.out = N)
  } else if (skew > 0) { # right-tailed
    x = seq(mode - 4.5 * sd, mode + 7 * sd, length.out = N)
  } else {               # symmetric
    x = seq(mode - 5 * sd, mode + 5 * sd, length.out = N)
  }
  # Discretize and evaluate
  #
  dx = dpearson(x, pfd) # un-normalized density
  # trapezoidal rule
  h = x[2] - x[1]
  # cumsum the middle ones except the first and the last
  cum = rep(0, N - 2)
  cum[1] = dx[2] * h
  for (i in 2:(N - 2)) {
    cum[i] = cum[i - 1] + dx[i + 1] * h # may too big number inside!
  }
  #
  px = rep(0, N)          # un-normalized cumulative probability
  px[2] = (dx[1] + dx[2]) * h / 2
  for (i in 3:N) {
    px[i] = (dx[1] + dx[i]) * h / 2 + cum[i - 2]
  }
  C = px[N] # constant
  # cat(sprintf("C = %20f\n", C))
  #
  # lb = x[1]; ub = x[length(x)]
  # logC_type = comp_logC(lb, ub, pfd, rel_err = 1e-7)
  # logC = logC_type$logC
  # cat(sprintf("C = %20f\n", exp(logC)))
  #
  Us = stats::runif(n)
  Xs = rep(0, n)
  for (i in 1:n) {
    Xs[i] = bd_refine(Us[i] * C, px, dx, x, pfd, pfd_tp)
  }
  return(Xs)
}

#' @rdname r_hest
#' @export
r_hest2 <- function(n, moms) {
  #
  # from moments to coefficients: a, c0, c1, c2, {c3 | c4 | c4, c5}
  #
  coefs = mom_to_coef(moms)
  #
  # Partial Fraction Decomposition
  #
  int_p = int_p8; bisect = bisect8
  pfd = PFDecomp4(coefs)
  # typically type 42: two distinct pairs of conjugate roots
  # but also type 44: 2 real, 2 complex: x1 != x2, (x3, x4)
  # printf("Partial Fraction Decomposition type: %i\n", pfd$type)
  #
  # Determine support of the distribution
  #
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4]) # mean, var, skew, kurt
  sd = sqrt(stdmoment[2])
  skew = stdmoment[3]
  #
  if (skew < 0) {        # left-tailed
    lb = -7*sd; ub = 4.5*sd
  } else if (skew > 0) { # right-tailed
    lb = -4.5*sd; ub = 7*sd
  } else {               # symmetric
    lb = -5*sd; ub = 5*sd
  }
  lbub = adjust_lb_ub(lb, ub, pfd)
  lb = lbub[1]; ub = lbub[2]
  #
  logC_type = comp_logC(lb, ub, pfd, rel_err = 1e-10)
  logC = logC_type$logC
  if (logC_type$type == -2) {
    printf("resort to PearsonDS::rpearson(), C = %E (too big or small)\n", exp(logC))
    return(PearsonDS::rpearson(n, moments=stdmom(moms[1:4])))
  }
  # printf("logC = %f, C = %E\n", logC, exp(logC))
  #
  Us = stats::runif(n)
  # ts = proc.time()
  mode = -coefs[1] # mode of the distribution
  x0s = stats::qnorm(Us, mean = mode, sd = sd)
  # x0s = PearsonDS::qpearson(Us, moments = stdmoment) # initial guess matters
  # setting 1: very fast
  # setting 2: consume 60 seconds more or less!
  # te = proc.time()
  # tt = te - ts
  # printf("The simulation of %i x0s consumes %.2f seconds\n", n, tt[[3]])
  #
  Xs = rep(0, n)
  for (i in 1:n) {
    U = Us[i]
    x0 = x0s[i]
    # if (U > 0.99) {
    #   xlb = max(lb, x0 - 2*sd)
    #   xub = min(ub, x0 + 2*sd)
    #   x = bisect(xlb, xub, logC, U, pfd, LB = lb, rel_err=1e-7)
    #   Xs[i] = x
    #   next
    # }
    # if (U < 0.01) {
    #   xlb = max(lb, x0 - 2*sd)
    #   xub = min(ub, x0 + 2*sd)
    #   x = bisect(xlb, xub, logC, U, pfd, LB = lb, rel_err=1e-7)
    #   Xs[i] = x
    #   next
    # }
    #
    j = 0
    Newton_Fail = FALSE
    #
    prob_next = int_p(lb, x0, pfd, logC, rel_err = 1e-10)
    repeat {
      if (Newton_Fail || x0 < lb || x0 > ub || j > 10) {
        # fmt = "bisect: U=%f, lb=% f, x0=% f, ub=%f, j=%i, Newton_Fail=%i\n"
        # printf(fmt, U, lb, x0, ub, j, Newton_Fail)
        xlb = max(lb, x0 - 2*sd)
        xub = min(ub, x0 + 2*sd)
        x = bisect(xlb, xub, logC, U, pfd, LB = lb, rel_err=1e-10)
        if (x <= lb || x >= ub) {
          x_exception(logC_type, lb, ub, U, x0, x, j, moms)
        }
        break
      }
      prob = prob_next
      f = prob - U
      df = exp(log_dpearson8(x0, pfd) - logC) # df = dpearson(x0, pfd)/C
      f_df = f / df
      # ddf = -df * (coefs[1] + x0)/sum(coefs[2:6] * (x0^(0:4)))
      ddf_df = - (coefs[1] + x0)/sum(coefs[2:6] * (x0^(0:4)))
      delta = 1 - 2 * (f_df) * ddf_df # delta = 1 - 2 * (f/df) * (ddf/df)
      #
      if (prob < 0 || prob > 1) {
        prob_exception(logC_type, prob, U, f_df, ddf_df, delta)
      }
      #
      if (delta >= 0 && ddf_df != 0) {
        x = x0 - (1 - sqrt(delta)) / ddf_df
      } else {
        Newton_Fail = TRUE; next
      }
      j = j + 1
      if (abs(x - x0) < 1e-10) { break }
      #
      if (x > x0) {
        prob_next = prob + int_p(x0, x, pfd, logC, rel_err = 1e-7)
      } else if (x < x0) {
        prob_next = prob - int_p(x, x0, pfd, logC, rel_err = 1e-7)
      }
      #
      x0 = x
    }
    Xs[i] = x
  }
  return(Xs)
}

printf <- function(fmt, ...) { cat(sprintf(fmt, ...)) }

x_exception <- function(logC_type, lb, ub, U, x0, x, j, moms) {
  if (x0 < lb || x0 > ub) { type_s = "TYPE1--out of bound" }
  else if (j > 10)        { type_s = "TYPE2--long-iteration" }
  else                    { type_s = "TYPE3--Newton-Condition-Failed" }
  #
  cat(sprintf("logC = %f, type = %d\n", logC_type$logC, logC_type$type))
  temp = "%sbisect()\n (lb, ub)= (%f, %f), U = %f, x0 = %f, x = %f\n"
  cat(sprintf(temp, type_s, lb, ub, U, x0, x))
  cat("moms = "); cat(moms, sep=","); cat("\n")
}

prob_exception <- function(logC_type, prob, U, f_df, ddf_df, delta) {
  cat(sprintf("logC = %f, type = %d\n", logC_type$logC, logC_type$type))
  temp = " f = %-e - %.7f, f/df = %-e, ddf/df = %-e, delta = %-e\n"
  cat(sprintf(temp, prob, U, f_df, ddf_df, delta))
}
