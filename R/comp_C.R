# for matching 6 moments
comp_logC_6 <- function(lb, ub, pfd, rel_err=1e-7) {
  h = ub - lb
  logu = log_dpearson6(0, pfd)               # upper
  logl = mean(log_dpearson6(c(lb, ub), pfd)) # lower
  #
  if (logl > logu) { # not always bell shape, sometimes U shape
    temp = logu; logu = logl; logl = temp; rm(temp)
  }
  #
  threshold1 = 1000
  threshold2 = 100
  if (logl < -threshold1 || logu > threshold1) {
    # cat("\nun-normalized density is either too big or too small!\n")
    # cat("Shift to PearsonDS::rpearson()\n")
    return(list(logC = NA, type = -2))
  } else if (logl < -threshold2 || logu > threshold2) {
    # cat("\nun-normalized density is either big or small!\n")
    # cat("Resort to Rmpfr numbers\n")
    logC = comp_logC_mpfr_6(lb, ub, pfd, rel_err)
    return(list(logC = logC, type = -1))
  } else {
    prob = stats::integrate(dpearson6, lb, ub, pfd, rel.tol = rel_err)
    return(list(logC = log(prob$value), type = 0))
  }
}

comp_logC_mpfr_6 <- function(lb, ub, pfd, rel_err=1e-7, max.ord = 18) {
  h = ub - lb
  logd_lb = log_dpearson6(lb, pfd)
  logd_ub = log_dpearson6(ub, pfd)
  density_lb = exp(Rmpfr::mpfr(logd_lb, precBits = 128))
  density_ub = exp(Rmpfr::mpfr(logd_ub, precBits = 128))
  Rp = c(0.5 * h * (density_lb + density_ub))
  #
  # cat(sprintf("i = 0, 2^i, C = %e\n", Rp[1]))
  for (i in 1:max.ord) {
    Rc = Rmpfr::mpfr(rep(0, i+1), precBits = 128)
    h = h / 2
    ep = 2^(i-1) # number of equidistant points
    x = seq(lb + h, ub, by = 2*h)
    logd = log_dpearson6(x, pfd)
    Rc[1] = 0.5*Rp[1] + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
                              Rmpfr::mpfr(h, precBits = 128))
    for (j in 1:i) {
      Rc[j+1] = Rc[j] + (Rc[j] - Rp[j]) / (4^j - 1)
    }
    #
    # cat(sprintf("i = %d, 2^i, C = %e\n", i, Rc[i+1]))
    #
    if (i > 1 && abs(Rp[i] - Rc[i+1])/Rc[i+1] < rel_err) {
      return(Rmpfr::asNumeric(log(Rc[i+1])))
    }
    Rp = Rc
  }
  return(Rmpfr::asNumeric(log(Rc[max.ord + 1])))
}


log_dpearson6 <- function(x, pfd3) {
  if (pfd3$type == 31) { return(log_dpearson61(x, pfd3)) }
  if (pfd3$type == 32) { return(log_dpearson62(x, pfd3)) }
  if (pfd3$type == 33) { return(log_dpearson63(x, pfd3)) }
  if (pfd3$type == 34) { return(log_dpearson64(x, pfd3)) }
  #
  stop("The type of pfd3 must be 31-34!")
}

log_dpearson61 <- function(x, pfd3) {
  # x1: real, (x2, x3): conjugate
  x1 = pfd3$x1; A1 = pfd3$A1
  x2 = pfd3$x2; A2 = pfd3$A2; B2 = pfd3$B2
  #
  t = x - Re(x2); m = Im(x2)
  pow1 = -A1 * log(abs(x - x1))
  pow2 = -A2 * log(t^2 + m^2)/2 - ((B2 + A2*Re(x2))/m) * atan(t/m)
  # pow1 + pow2
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson62 <- function(x, pfd3) {
  # x1 = x2 = x3: real
  x1 = pfd3$x1; A2 = pfd3$A2; A3 = pfd3$A3
  #
  pow1 = A2 / (x - x1)
  pow2 = A3 / (2 * (x- x1)^2)
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson63 <- function(x, pfd3) {
  # x1 = x2 != x3: real
  x1 = pfd3$x1; A1 = pfd3$A1; A2 = pfd3$A2
  x3 = pfd3$x3; A3 = pfd3$A3
  pow1 = -A1 * log(abs(x - x1)) + A2 / (x - x1)
  pow2 = -A3 * log(abs(x - x3))
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson64 <- function(x, pfd3) {
  # x1 != x2 != x3: real
  x1 = pfd3$x1; A1 = pfd3$A1
  x2 = pfd3$x2; A2 = pfd3$A2
  x3 = pfd3$x3; A3 = pfd3$A3
  pow1 = -A1 * log(abs(x - x1))
  pow2 = -A2 * log(abs(x - x2))
  pow3 = -A3 * log(abs(x - x3))
  # return(exp(pow1 + pow2 + pow3))
  return(pow1 + pow2 + pow3)
}

# for matching 8 moments

comp_logC <- function(lb, ub, pfd, rel_err=1e-7) {
  h = ub - lb
  logu = log_dpearson8(0, pfd)               # upper
  logl = mean(log_dpearson8(c(lb, ub), pfd)) # lower
  #
  if (logl > logu) { # not always bell shape, sometimes U shape
    temp = logu; logu = logl; logl = temp; rm(temp)
  }
  #
  threshold1 = 1000
  threshold2 = 100
  if (logl < -threshold1 || logu > threshold1) {
    # cat("\nun-normalized density is either too big or too small!\n")
    # cat("Shift to PearsonDS::rpearson()\n")
    return(list(logC = NA, type = -2))
  } else if (logl < -threshold2 || logu > threshold2) {
    # cat("\nun-normalized density is either big or small!\n")
    # cat("Resort to Rmpfr numbers\n")
    logC = comp_logC_mpfr(lb, ub, pfd, rel_err)
    return(list(logC = logC, type = -1))
  } else {
    prob = stats::integrate(dpearson8, lb, ub, pfd, rel.tol = rel_err)
    return(list(logC = log(prob$value), type = 0))
  }
}

comp_logC_mpfr <- function(lb, ub, pfd, rel_err=1e-7, max.ord = 18) {
  h = ub - lb
  logd_lb = log_dpearson8(lb, pfd)
  logd_ub = log_dpearson8(ub, pfd)
  density_lb = exp(Rmpfr::mpfr(logd_lb, precBits = 128))
  density_ub = exp(Rmpfr::mpfr(logd_ub, precBits = 128))
  Rp = c(0.5 * h * (density_lb + density_ub))
  #
  # cat(sprintf("i = 0, 2^i, C = %e\n", Rp[1]))
  for (i in 1:max.ord) {
    Rc = Rmpfr::mpfr(rep(0, i+1), precBits = 128)
    h = h / 2
    ep = 2^(i-1) # number of equidistant points
    x = seq(lb + h, ub, by = 2*h)
    logd = log_dpearson8(x, pfd)
    Rc[1] = 0.5*Rp[1] + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
                            Rmpfr::mpfr(h, precBits = 128))
    for (j in 1:i) {
      Rc[j+1] = Rc[j] + (Rc[j] - Rp[j]) / (4^j - 1)
    }
    #
    # cat(sprintf("i = %d, 2^i, C = %e\n", i, Rc[i+1]))
    #
    if (i > 1 && abs(Rp[i] - Rc[i+1])/Rc[i+1] < rel_err) {
      return(Rmpfr::asNumeric(log(Rc[i+1])))
    }
    Rp = Rc
  }
  return(Rmpfr::asNumeric(log(Rc[max.ord + 1])))
}

log_dpearson8 <- function(x, pfd) {
  if (pfd$type == 41) { return(log_dpearson81(x, pfd)) }
  if (pfd$type == 42) { return(log_dpearson82(x, pfd)) }
  if (pfd$type == 43) { return(log_dpearson83(x, pfd)) }
  if (pfd$type == 44) { return(log_dpearson84(x, pfd)) }
  if (pfd$type == 45) { return(log_dpearson85(x, pfd)) }
  if (pfd$type == 46) { return(log_dpearson86(x, pfd)) }
  if (pfd$type == 47) { return(log_dpearson87(x, pfd)) }
  if (pfd$type == 48) { return(log_dpearson88(x, pfd)) }
  if (pfd$type == 49) { return(log_dpearson89(x, pfd)) }
  #
  stop("The type of pfd4 must be 41-49!")
}

log_dpearson81 <- function(x, pfd) {
  # pfd$type = 41
  A = pfd$A1; B = pfd$B1
  #
  t = x - Re(pfd$x1); m = Im(pfd$x1); B = B + A * Re(pfd$x1)
  #
  pow = A/(2*(t^2 + m^2)) - (B/(2*m^2)) * (t/(t^2 + m^2))
  pow = pow - (B/(2*m^3)) * atan(t/m)
  # return(exp(pow))
  return(pow)
}

log_dpearson82 <- function(x, pfd) {
  # pfd$type = 42
  x1 = pfd$x1; A1 = pfd$A1; B1 = pfd$B1
  x2 = pfd$x2; A2 = pfd$A2; B2 = pfd$B2
  #
  t1 = x - Re(x1); m1 = Im(x1)
  t2 = x - Re(x2); m2 = Im(x2)
  pow1 = -A1 * log(t1^2 + m1^2)/2 - ((B1 + A1*Re(x1))/m1) * atan(t1/m1)
  pow2 = -A2 * log(t2^2 + m2^2)/2 - ((B2 + A2*Re(x2))/m2) * atan(t2/m2)
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson83 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; B3 = pfd$B3
  #
  t3 = x - Re(x3); m3 = Im(x3)
  #
  pow1 = -A1 * log(abs(x - x1)) + A2/(x-x1)
  pow2 = -A3 * log(t3^2 + m3^2)/2 - ((B3 + A3*Re(x3))/m3) * atan(t3/m3)
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson84 <- function(x, pfd) {
  # pfd$type = 44
  x1 = pfd$x1; A1 = pfd$A1
  x2 = pfd$x2; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; B3 = pfd$B3
  #
  t3 = x - Re(x3); m3 = Im(x3);
  #
  pow1 = -A1 * log(abs(x - x1)) - A2 * log(abs(x - x2))
  pow2 = -A3 * log(t3^2 + m3^2)/2 - ((B3 + A3*Re(x3))/m3) * atan(t3/m3)
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson85 <- function(x, pfd) {
  x1 = pfd$x1; A3 = pfd$A3; A4 = pfd$A4
  #
  pow = A3 /(2*(x - x1)^2) + A4 /(3*(x - x1)^3)
  # return(exp(pow))
  return(pow)
}

log_dpearson86 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x - x1)) + A2/(x - x1) + A3/(2*(x- x1)^2)
  pow2 = -A4 * log(abs(x - x4))
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson87 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) + A2/(x-x1)
  pow2 = -A3 * log(abs(x-x3)) + A4/(x-x3)
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson88 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) + A2/(x-x1)
  pow2 = -A3 * log(abs(x-x3)) - A4 * log(abs(x-x4))
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

log_dpearson89 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1
  x2 = pfd$x2; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) - A2 * log(abs(x-x2))
  pow2 = -A3 * log(abs(x-x3)) - A4 * log(abs(x-x4))
  # return(exp(pow1 + pow2))
  return(pow1 + pow2)
}

# #
# logd_lb = log_dpearson8(lb, pfd)
# logd_ub = log_dpearson8(ub, pfd)
# denlb = exp(Rmpfr::mpfr(logd_lb, precBits = 128))
# denub = exp(Rmpfr::mpfr(logd_ub, precBits = 128))
# h = ub - lb
# prob0 = h * (denlb + denub) / 2
# # first fine approximation
# N = 2^5                             # number of segments
# h = (ub - lb) / N; prob = prob0 / N # update h, prob
# #
# x = seq(lb, ub, length.out = N+1)
# logd = log_dpearson8(x[2:N], pfd)
# prob = prob + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
#                     Rmpfr::mpfr(h, precBits = 128))
# # count = length(x)
# #
# # temp = "C = %e, C = %e, h = %f, (lb, ub) = (%.5f, %.5f)\n"
# # cat(sprintf(temp, prob0, prob, h, lb, ub))
# err = abs(prob - prob0)
# while (TRUE) {
#   if (err/prob < rel_err || h < 1e-5) {break}
#   prob0 = prob
#   h = h / 2
#   x = seq(lb + h, ub, by = 2*h)
#   logd = log_dpearson8(x, pfd)
#   prob = prob / 2
#   prob = prob + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
#                       Rmpfr::mpfr(h, precBits = 128))
#   err = abs(prob - prob0)
#   # count = count + length(x)
#   # cat(sprintf("C = %e\n", prob))
# }
# # cat(sprintf("%d evaluations of dpearson(), C = %e\n", count, prob))
# logC = log(prob)
# return(Rmpfr::asNumeric(logC))

# #
# Id <- Rmpfr::integrateR(dpearson8,
#                         lb, #Rmpfr::mpfr(lb, precBits = 128),
#                         ub, #Rmpfr::mpfr(ub, precBits = 128),
#                         pfd, rel.tol = rel_err, verbose = F)
# return(list(logC = log(Id$value), type = 0))
# #
# # normal precision: double number
# h = ub - lb
# prob0 = h*(exp(log_dpearson8(lb, pfd)) + exp(log_dpearson8(ub, pfd)))/2
# # first fine approximation
# N = 2^5                             # number of segments
# h = (ub - lb) / N; prob = prob0 / N # update h, prob
# #
# x = seq(lb, ub, length.out = N+1)
# prob = prob + sum(exp(log_dpearson8(x[2:N], pfd)) * h)
# # count = length(x)
# #
# # temp = "\nC = %e, C = %e, h = %f, (lb, ub) = (%.5f, %.5f)\n"
# # cat(sprintf(temp, prob0, prob, h, lb, ub))
# err = abs(prob - prob0)
# while (TRUE) {
#   if (err/prob < rel_err || h < 1e-5) {break}
#   prob0 = prob
#   h = h/2
#   x = seq(lb + h, ub, by = 2*h)
#   prob = prob/2  + sum(exp(log_dpearson8(x, pfd)) * h)
#   err = abs(prob - prob0)
#   # count = count + length(x)
#   # cat(sprintf("C = %e, h = %f, abs_err = %f, rel_err = %f\n", prob, h, err, err/prob))
# }
# if (err/prob > rel_err) {
#   x = seq(lb, ub, 0.001)
#   par(mfrow=c(2,1))
#   par(cex=0.8, mai=c(0.6,0.8,0.6,0.5))
#   plot(x, dpearson8(x, pfd), type = 'l', col='blue', ylab = "density",
#        main = "un-normalized p(x)")
#   plot(x, log_dpearson8(x, pfd), type = 'l', col='red', ylab = "log-density",
#        main = "un-normalized log(p(x))")
#   print(pfd)
#   #
#   h = ub - lb
#   prob0 = h*(exp(log_dpearson8(lb, pfd)) + exp(log_dpearson8(ub, pfd)))/2
#   # first fine approximation
#   N = 2^5                             # number of segments
#   h = (ub - lb) / N; prob = prob0 / N # update h, prob
#   #
#   x = seq(lb, ub, length.out = N+1)
#   prob = prob + sum(exp(log_dpearson8(x[2:N], pfd)) * h)
#   count = length(x)
#   #
#   temp = "\nC = %e, C = %e, h = %f, (lb, ub) = (%.5f, %.5f)\n"
#   cat(sprintf(temp, prob0, prob, h, lb, ub))
#   err = abs(prob - prob0)
#   while (TRUE) {
#     if (err/prob < rel_err || h < 1e-5) {break}
#     prob0 = prob
#     h = h/2
#     x = seq(lb + h, ub, by = 2*h)
#     prob = prob/2  + sum(exp(log_dpearson8(x, pfd)) * h)
#     err = abs(prob - prob0)
#     count = count + length(x)
#     temp = "C = %e, h = %f, abs_err = %f, \trel_err = %.7f\n"
#     cat(sprintf(temp, prob, h, err, err/prob))
#   }
#   cat(sprintf("%d evaluations of dpearson(), \tC = %e\n", count, prob))
# }
# return(list(logC = log(prob), type = 0))
