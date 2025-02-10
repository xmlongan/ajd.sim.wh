#' Dynamically Compute the Cumulative Probability
#'
#' @name int_p
#' @description
#' Dynamically compute the cumulative probability over an interval
#'
#' @param lb lower bound
#' @param ub upper bound
#' @param pfd3,pfd4,pfd5 solution of Partial Fraction Decomposition
#' @param eps accuracy to achieve
#' @param logC log of the constant \eqn{C} (density normalization).
#' @param rel_err relative error to accept
#' @param max.ord maximum order of the Romberg Algorithm
#'
#' @return scalar
#' @export
#'
#' @examples
#' lb = -0.5952014; ub = 0.3030635
int_p8 <- function(lb, ub, pfd4, logC, rel_err=1e-7, max.ord = 18) {
  if (abs(logC) < 100) {
    prob = stats::integrate(dpearson8, lb, ub, pfd4, rel.tol = rel_err)
    prob = prob$value
    return(exp(log(prob) - logC)) # return(prob/exp(logC))
  } else { # special case large number or small number: Romberg Algorithm
    h = ub - lb
    logd_lb = log_dpearson8(lb, pfd4)
    logd_ub = log_dpearson8(ub, pfd4)
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
      logd = log_dpearson8(x, pfd4)
      Rc[1] = 0.5*Rp[1] + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
                              Rmpfr::mpfr(h, precBits = 128))
      for (j in 1:i) {
        Rc[j+1] = Rc[j] + (Rc[j] - Rp[j]) / (4^j - 1)
      }
      #
      # cat(sprintf("i = %d, 2^i, C = %e\n", i, Rc[i+1]))
      #
      if (i > 1 && abs(Rp[i] - Rc[i+1])/Rc[i+1] < rel_err) {
        return(Rmpfr::asNumeric(exp(log(Rc[i+1]) - logC)))
      }
      Rp = Rc
    }
    return(Rmpfr::asNumeric(exp(log(Rc[max.ord + 1]) - logC)))
  }
}

# # special case
# logd_lb = log_dpearson8(lb, pfd4)
# logd_ub = log_dpearson8(ub, pfd4)
# denlb = exp(Rmpfr::mpfr(logd_lb, precBits = 128))
# denub = exp(Rmpfr::mpfr(logd_ub, precBits = 128))
# h = ub - lb
# prob0 = h * (denlb + denub) / 2
# # first fine approximation
# N = 2^5                             # number of segments
# h = (ub - lb) / N; prob = prob0 / N # update h, prob
# #
# x = seq(lb, ub, length.out = N+1)
# logd = log_dpearson8(x[2:N], pfd4)
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
#   logd = log_dpearson8(x, pfd4)
#   prob = prob / 2
#   prob = prob + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
#                       Rmpfr::mpfr(h, precBits = 128))
#   err = abs(prob - prob0)
#   # count = count + length(x)
#   # cat(sprintf("C = %e\n", prob))
# }
# # cat(sprintf("%d evaluations of dpearson(), C = %e\n", count, prob))
# prob_normalized = prob/exp(Rmpfr::mpfr(logC, precBits = 128))
# return(Rmpfr::asNumeric(prob_normalized))

# # Id <- Rmpfr::integrateR(dpearson8,
# #                         lb,#Rmpfr::mpfr(lb, precBits = 128),
# #                         ub,#Rmpfr::mpfr(ub, precBits = 128),
# #                         pfd4, rel.tol = 1e-7, abs.tol = 1e-5, verbose = T)
# # return(Id$value/exp(logC))
# h = ub - lb
# #
# prob0 = h*(exp(log_dpearson8(lb, pfd4)) + exp(log_dpearson8(ub, pfd4)))/2
# # first fine approximation
# N = ifelse(h > 0.001, 4, 2)         # number of segments
# h = (ub - lb) / N; prob = prob0 / N # update h, prob
# #
# x = seq(lb, ub, length.out = N+1)
# prob = prob + sum(exp(log_dpearson8(x[2:N], pfd4)) * h)
# count = length(x)
# #
# # cat(sprintf("(lb, ub) = (%.5f, %.5f)\n", lb, ub))
# # temp = "First estimation:\n2 points, \tprob = %f\n%d points, \tprob = %f\n"
# # cat(sprintf(temp, prob0, N+1, prob))
# #
# err = abs(prob - prob0)
# while (TRUE) {
#   if (err < abs_err || err/prob < rel_err || h < 1e-5) {break}
#   prob0 = prob
#   h = h / 2
#   x = seq(lb + h, ub, by = 2*h)
#   prob = prob/2 + sum(exp(log_dpearson8(x, pfd4)) * h)
#   # prob = prob/2 + mean(exp(log_dpearson8(x, pfd4))) * h * length(x)
#   #
#   # cat(sprintf("prob =%f, prob0=%f\n", prob, prob0))
#   err = abs(prob - prob0)
#   count = count + length(x)
# }
# # cat(sprintf("total times of dpearson() evaluation: %d\n", count))
# # if (prob < 0 || prob > 1) {
# #   cat(sprintf("\t%d evaluations of dpearson(),\tprob = %.5f\n", count, prob))
# # }
# return(prob/exp(logC))

#' @rdname int_p
#' @export
int_p6 <- function(lb, ub, pfd3, logC, rel_err=1e-7, max.ord = 18) {
  if (abs(logC) < 100) {
    prob = stats::integrate(dpearson6, lb, ub, pfd3, rel.tol = rel_err)
    prob = prob$value
    return(exp(log(prob) - logC)) # return(prob/exp(logC))
  } else { # special case large number or small number: Romberg Algorithm
    h = ub - lb
    logd_lb = log_dpearson6(lb, pfd3)
    logd_ub = log_dpearson6(ub, pfd3)
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
      logd = log_dpearson6(x, pfd3)
      Rc[1] = 0.5*Rp[1] + sum(exp(Rmpfr::mpfr(logd, precBits = 128)) *
                                Rmpfr::mpfr(h, precBits = 128))
      for (j in 1:i) {
        Rc[j+1] = Rc[j] + (Rc[j] - Rp[j]) / (4^j - 1)
      }
      #
      # cat(sprintf("i = %d, 2^i, C = %e\n", i, Rc[i+1]))
      #
      if (i > 1 && abs(Rp[i] - Rc[i+1])/Rc[i+1] < rel_err) {
        return(Rmpfr::asNumeric(exp(log(Rc[i+1]) - logC)))
      }
      Rp = Rc
    }
    return(Rmpfr::asNumeric(exp(log(Rc[max.ord + 1]) - logC)))
  }
}

#' @rdname int_p
#' @export
int_p10 <- function(lb, ub, pfd5, eps) {
  h = ub - lb
  prob0 = h * (dpearson10(lb, pfd5) + dpearson10(ub, pfd5)) / 2
  # first fine approximation
  N = ifelse(h > 0.001, 4, 2)     # number of segments
  h = (ub - lb) / N; prob = prob0 / N # update h, prob
  x = lb + h
  while (x < ub) {
    prob = prob + dpearson10(x, pfd5) * h
    x = x + 2 * h
  }
  #
  while (abs(prob - prob0) > eps && h > eps * 0.01) {
    prob0 = prob
    h = h / 2; prob = prob0 / 2
    x = lb + h
    while (x < ub) {
      prob = prob + dpearson10(x, pfd5) * h
      x = x + 2 * h
    }
  }
  return(prob)
}
