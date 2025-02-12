#' Bisection to Find the Target Quantile
#'
#' @name bisect
#' @description
#' Use bisection to find the target quantile
#'
#' @param lb lower bound of the quantile
#' @param ub upper bound of the quantile
#' @param U target probability
#' @param pfd3,pfd4,pfd5 solution of Partial Fraction Decomposition
#' @param LB lower bound of the support of the distribution
#' @param logC log of the constant \eqn{C} (density normalization).
#' @param rel_err relative error to accept
#' @param eps accuracy to achieve
#'
#' @return scalar
#' @export
#'
#' @examples
#' lb = -0.01; ub = 0.01
bisect6 <- function(lb, ub, logC, U, pfd3, LB, rel_err=1e-5) {
  x = (lb + ub) / 2
  while (ub - lb > 1e-5) {
    x = (lb + ub) / 2
    prob = int_p6(LB, x, pfd3, logC, rel_err)
    if (abs(prob - U) < 1e-5) { break }
    else if (prob < U) { lb = x }
    else { ub = x }
  }
  return(x)
}

#' @rdname bisect
#' @export
bisect8 <- function(lb, ub, logC, U, pfd4, LB, rel_err=1e-5) {
  x = (lb + ub) / 2
  prob_next = int_p8(LB, x, pfd4, logC, rel_err)
  while (ub - lb > 1e-7) {
    x = (lb + ub) / 2
    # prob = int_p8(LB, x, pfd4, logC, rel_err)
    prob = prob_next
    if (abs(prob - U) < 1e-7) {
      break
    } else if (prob < U) {
      lb = x
      prob_next = prob + int_p8(lb, (lb+ub)/2, pfd4, logC, rel_err)
    } else {
      ub = x
      prob_next = prob - int_p8((lb+ub)/2, ub, pfd4, logC, rel_err)
    }
    # x = (lb + ub) / 2
    # cat(sprintf("U = %f, prob = %f, x = %f\n", U, prob, x))
  }
  return(x)
}

#' @rdname bisect
#' @export
bisect10 <- function(lb, ub, U, pfd5, LB, eps) {
  x = (lb + ub) / 2
  while (ub - lb > eps) {
    x = (lb + ub) / 2
    prob = int_p10(LB, x, pfd5, eps)
    if (prob < U) { lb = x }
    else if (prob > U) { ub = x }
    else { break }
  }
  return(x)
}
