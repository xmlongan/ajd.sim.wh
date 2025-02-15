rpearson8 <- function(n, moms) {
  #
  # preliminary process
  #
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4]) # mean=0, var, skew, kurt
  sd = sqrt(stdmoment[2])
  skew = stdmoment[3]
  #
  # from moments to coefficients: a, c0, c1, c2, {c3 | c4 | c4, c5}
  #
  coefs = mom_to_coef(moms)
  mode = -coefs[1] # mode of the distribution
  #
  # Partial Fraction Decomposition (PFD)
  #
  pfd = PFDecomp4(coefs)
  # typically type 42: two distinct pairs of conjugate roots
  # but also  type 44: 2 real, 2 complex: x1 != x2, (x3, x4)
  # printf("Partial Fraction Decomposition type: %i\n", pfd$type)
  #
  # Determine support of the distribution
  #
  if (skew < 0) {        # left-tailed
    lb = mode - 7*sd; ub = mode + 4.5*sd
  } else if (skew > 0) { # right-tailed
    lb = mode - 4.5*sd; ub = mode + 7*sd
  } else {               # symmetric
    lb = mode - 5*sd; ub = mode + 5*sd
  }
  # adjust the support when necessary:
  #   real roots (from PFD) should not be included in the support
  lbub = adjust_lb_ub(lb, ub, pfd)
  if (lbub[1] != lb || lbub[2] != ub) {
    fmt = "original support (%f, %f), adjusted as (%f, %f)"
    printf(fmt, lb, ub, lbub[1], lbub[2])
  }
  lb = lbub[1]; ub = lbub[2]
  #
  # pre-computation
  #
  # q = seq(lb, ub, sd)   # 11 + 1 quantiles
  h = sd/100
  q = seq(lb, ub, h)
  p = rep(0, length(q)) # sub probabilities (un-normalized)
  for (i in 2:length(q)) {
    subprob = stats::integrate(dpearson8, q[i-1], q[i], pfd, rel.tol = 1e-7)
    p[i] = subprob$value
  }
  C = sum(p)
  logC = log(C)
  # normalize the probability
  p = cumsum(p)/C
  #
  Us = stats::runif(n)
  #
  # initial guess
  #
  # ts = proc.time()
  x0s = stats::qnorm(Us, mean = mode, sd = sd)
  x0s[x0s > ub] = ub
  x0s[x0s < lb] = lb
  # x0s = PearsonDS::qpearson(Us, moments = stdmoment) # initial guess matters
  # setting 1: very fast
  # setting 2: consume 60 seconds more or less!
  # te = proc.time()
  # tt = te - ts
  # printf("The simulation of %i x0s consumes %.2f seconds\n", n, tt[[3]])
  #
  IS = Us < 0.01     # small
  IL = Us > 0.99     # large
  IN = (!IS) & (!IL) # normal
  #
  US = Us[IS]; x0S = x0s[IS]; XS = rep(0, sum(IS))
  UL = Us[IL]; x0L = x0s[IL]; XL = rep(0, sum(IL))
  UN = Us[IN]; x0N = x0s[IN]; XN = rep(0, sum(IN))
  #
  rm(IS, IL, IN, Us, x0s)
  #
  for (i in 1:length(XS)) {
    U = US[i]
    x0 = x0S[i]
    xlb = max(lb, x0 - 2*sd)
    xub = min(ub, x0 + 2*sd)
    XS[i] = bisect8(xlb, xub, logC, U, pfd, LB = lb, rel_err=1e-7)
  }
  for (i in 1:length(XL)) {
    U = UL[i]
    x0 = x0L[i]
    xlb = max(lb, x0 - 2*sd)
    xub = min(ub, x0 + 2*sd)
    XL[i] = bisect8(xlb, xub, logC, U, pfd, LB = lb, rel_err=1e-7)
  }
  #
  for (i in 1:length(XN)) {
    U = UN[i]
    x0 = x0N[i]
    j = 0
    Newton_Fail = FALSE
    #
    i_xlb = floor((x0 - lb)/h) + 1
    prob_next = p[i_xlb] + int_p8(q[i_xlb], x0, pfd, logC, rel_err=1e-7)
    repeat {
      if (Newton_Fail || x0 < lb || x0 > ub || j > 10 ) {
        printf("Newton_Fail = %i, resort to bisect\n", Newton_Fail)
        xlb = max(lb, x0-2*sd)
        xub = min(ub, x0+2*sd)
        x = bisect8(xlb, xub, logC, U, pfd, LB = lb, rel_err = 1e-7)
        break
      }
      prob = prob_next
      f = prob - U
      df = exp(log_dpearson8(x0, pfd) - logC)
      f_df = f / df
      ddf_df = - (coefs[1] + x0)/sum(coefs[2:6] * (x0^(0:4)))
      delta = 1 - 2 * (f_df) * ddf_df
      if (prob < 0 || prob > 1) {
        fmt = "C=%E, f=%-e - %.7f, f/df=%-e, ddf/df=%-e, delta=%-e\n"
        printf(fmt, C, prob, U, f_df, ddf_df, delta)
      }
      #
      if (delta >= 0 && ddf_df != 0) {
        x = x0 - (1 - sqrt(delta)) / ddf_df
      } else {
        Newton_Fail = TRUE; next
      }
      j = j + 1
      if (abs(x - x0) < 1e-5) {break}
      #
      if (x > x0) {
        prob_next = prob + int_p8(x0, x, pfd, logC, rel_err = 1e-5)
      } else if (x < x0) {
        prob_next = prob - int_p8(x, x0, pfd, logC, rel_err = 1e-5)
      }
      #
      x0 = x
    }
    XN[i] = x
  }
  return(c(XS, XN, XL))
}
