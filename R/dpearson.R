#' Un-normalized Density of Pearson Distribution
#'
#' @name dpearson
#' @description
#' Un-normalized density of Pearson Distribution
#'
#' @param pfd3 list of x1(real), A1;        x2(complex), A2, B2
#' @param pfd5 list of x1(real), A1;        x2(complex), A2, B2;
#'  x3(complex), A3, B3
#' @param x vector of quantiles
#' @param pfd Partial Fraction Decomposition encoded in a list
#'
#' @return scalar or vector of the density|densities
#' @export
#'
#' @examples
#' x = 0
dpearson6 <- function(x, pfd3) {
  if (pfd3$type == 31) { return(dpearson61(x, pfd3)) }
  if (pfd3$type == 32) { return(dpearson62(x, pfd3)) }
  if (pfd3$type == 33) { return(dpearson63(x, pfd3)) }
  if (pfd3$type == 34) { return(dpearson64(x, pfd3)) }
  #
  stop("The type of pfd3 must be 31-34!")
}

dpearson61 <- function(x, pfd3) {
  # x1: real, (x2, x3): conjugate
  x1 = pfd3$x1; A1 = pfd3$A1
  x2 = pfd3$x2; A2 = pfd3$A2; B2 = pfd3$B2
  #
  t = x - Re(x2); m = Im(x2)
  pow1 = -A1 * log(abs(x - x1))
  pow2 = -A2 * log(t^2 + m^2)/2 - ((B2 + A2*Re(x2))/m) * atan(t/m)
  # pow1 + pow2
  return(exp(pow1 + pow2))
}

dpearson62 <- function(x, pfd3) {
  # x1 = x2 = x3: real
  x1 = pfd3$x1; A2 = pfd3$A2; A3 = pfd3$A3
  #
  pow1 = A2 / (x - x1)
  pow2 = A3 / (2 * (x- x1)^2)
  return(exp(pow1 + pow2))
}

dpearson63 <- function(x, pfd3) {
  # x1 = x2 != x3: real
  x1 = pfd3$x1; A1 = pfd3$A1; A2 = pfd3$A2
  x3 = pfd3$x3; A3 = pfd3$A3
  pow1 = -A1 * log(abs(x - x1)) + A2 / (x - x1)
  pow2 = -A3 * log(abs(x - x3))
  return(exp(pow1 + pow2))
}

dpearson64 <- function(x, pfd3) {
  # x1 != x2 != x3: real
  x1 = pfd3$x1; A1 = pfd3$A1
  x2 = pfd3$x2; A2 = pfd3$A2
  x3 = pfd3$x3; A3 = pfd3$A3
  pow1 = -A1 * log(abs(x - x1))
  pow2 = -A2 * log(abs(x - x2))
  pow3 = -A3 * log(abs(x - x3))
  return(exp(pow1 + pow2 + pow3))
}


#' @rdname dpearson
#' @export
dpearson8 <- function(x, pfd) {
  if (pfd$type == 41) { return(dpearson81(x, pfd)) }
  if (pfd$type == 42) { return(dpearson82(x, pfd)) }
  if (pfd$type == 43) { return(dpearson83(x, pfd)) }
  if (pfd$type == 44) { return(dpearson84(x, pfd)) }
  if (pfd$type == 45) { return(dpearson85(x, pfd)) }
  if (pfd$type == 46) { return(dpearson86(x, pfd)) }
  if (pfd$type == 47) { return(dpearson87(x, pfd)) }
  if (pfd$type == 48) { return(dpearson88(x, pfd)) }
  if (pfd$type == 49) { return(dpearson89(x, pfd)) }
  #
  stop("The type of pfd4 must be 41-49!")
}

dpearson81 <- function(x, pfd) {
  # pfd$type = 41
  A = pfd$A1; B = pfd$B1
  #
  t = x - Re(pfd$x1); m = Im(pfd$x1); B = B + A * Re(pfd$x1)
  #
  pow = A/(2*(t^2 + m^2)) - (B/(2*m^2)) * (t/(t^2 + m^2))
  pow = pow - (B/(2*m^3)) * atan(t/m)
  return(exp(pow))
}

dpearson82 <- function(x, pfd) {
  # pfd$type = 42
  x1 = pfd$x1; A1 = pfd$A1; B1 = pfd$B1
  x2 = pfd$x2; A2 = pfd$A2; B2 = pfd$B2
  #
  t1 = x - Re(x1); m1 = Im(x1)
  t2 = x - Re(x2); m2 = Im(x2)
  pow1 = -A1 * log(t1^2 + m1^2)/2 - ((B1 + A1*Re(x1))/m1) * atan(t1/m1)
  pow2 = -A2 * log(t2^2 + m2^2)/2 - ((B2 + A2*Re(x2))/m2) * atan(t2/m2)
  return(exp(pow1 + pow2))
}

dpearson83 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; B3 = pfd$B3
  #
  t3 = x - Re(x3); m3 = Im(x3)
  #
  pow1 = -A1 * log(abs(x - x1)) + A2/(x-x1)
  pow2 = -A3 * log(t3^2 + m3^2)/2 - ((B3 + A3*Re(x3))/m3) * atan(t3/m3)
  return(exp(pow1 + pow2))
}

dpearson84 <- function(x, pfd) {
  # pfd$type = 44
  x1 = pfd$x1; A1 = pfd$A1
  x2 = pfd$x2; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; B3 = pfd$B3
  #
  t3 = x - Re(x3); m3 = Im(x3);
  #
  pow1 = -A1 * log(abs(x - x1)) - A2 * log(abs(x - x2))
  pow2 = -A3 * log(t3^2 + m3^2)/2 - ((B3 + A3*Re(x3))/m3) * atan(t3/m3)
  return(exp(pow1 + pow2))
}

dpearson85 <- function(x, pfd) {
  x1 = pfd$x1; A3 = pfd$A3; A4 = pfd$A4
  #
  pow = A3 /(2*(x - x1)^2) + A4 /(3*(x - x1)^3)
  return(exp(pow))
}

dpearson86 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x - x1)) + A2/(x - x1) + A3/(2*(x- x1)^2)
  pow2 = -A4 * log(abs(x - x4))
  return(exp(pow1 + pow2))
}

dpearson87 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) + A2/(x-x1)
  pow2 = -A3 * log(abs(x-x3)) + A4/(x-x3)
  return(exp(pow1 + pow2))
}

dpearson88 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) + A2/(x-x1)
  pow2 = -A3 * log(abs(x-x3)) - A4 * log(abs(x-x4))
  return(exp(pow1 + pow2))
}

dpearson89 <- function(x, pfd) {
  x1 = pfd$x1; A1 = pfd$A1
  x2 = pfd$x2; A2 = pfd$A2
  x3 = pfd$x3; A3 = pfd$A3
  x4 = pfd$x4; A4 = pfd$A4
  #
  pow1 = -A1 * log(abs(x-x1)) - A2 * log(abs(x-x2))
  pow2 = -A3 * log(abs(x-x3)) - A4 * log(abs(x-x4))
  return(exp(pow1 + pow2))
}


#' @rdname dpearson
#' @export
dpearson10 <- function(x, pfd5) {
  # roots: 1 real, 2 pairs of conjugate complex numbers
  x1 = pfd5$x1; x2 = pfd5$x2; x3 = pfd5$x3
  A1 = pfd5$A1; A2 = pfd5$A2; B2 = pfd5$B2; A3 = pfd5$A3; B3 = pfd5$B3
  #
  t2 = x - Re(x2); m2 = Im(x2)
  t3 = x - Re(x3); m3 = Im(x3)
  pow1 = -A1 * log(abs(x - x1))
  pow2 = -A2 * log(t2^2 + m2^2)/2 - ((B2 + A2*Re(x2))/m2) * atan(t2/m2)
  pow3 = -A3 * log(t3^2 + m3^2)/2 - ((B3 + A3*Re(x3))/m3) * atan(t3/m3)
  return(exp(pow1 + pow2 + pow3))
}
