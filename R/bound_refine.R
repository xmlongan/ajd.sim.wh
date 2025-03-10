#' Bound and Refine to Find the Target Quantile
#'
#' @name bound_refine
#' @description
#' Bound and refine the target quantile, `bd_refine`.
#' Get the lower and upper bound of the target quantile, `bound`.
#' Refine within the bounds to find the target quantile, `refine`.
#'
#' @param UC random number U times the normal constant C
#' @param px vector of un-normalized cumulative probability
#' @param dx vector of un-normalized density
#' @param x vector of quantiles
#' @param lb lower bound of the target quantile
#' @param ub upper bound of the target quantile
#' @param px_lb un-normalized cumulative probability of the lower bound `lb`
#' @param dx_lb un-normalized density of the lower bound `lb`
#' @param pfd list of Partial Fraction Decomposition results
#' @param pfd_tp type of `pfd`:
#'  - 3: list of x1(real), A1;        x2(complex), A2, B2
#'  - 4: list of x1(complex), A1, B1; x2(complex), A2, B2
#'  - 5: list of x1(real), A1;        x2(complex), A2, B2; x3(complex), A3, B3
#'
#' @return
#'  - `bd_refine`: scalar, quantile point \eqn{x} such that \eqn{P(X\le x)=U}.
#'  - `bound`: vector (lower_bound_index, upper_bound_index, target_index).
#'  - `refine`: scalar as that in `bd_refine`.
#' @export
#'
#' @examples
#' U = runif(1)
bound_refine <- function(UC, px, dx, x, pfd, pfd_tp) {
  lui = bound(UC, px)
  # lower bound index, upper bound index, target index
  l = lui[1]; u = lui[2]; i = lui[3]
  if (i != 0) return(x[i])
  #
  lb = x[l]; ub = x[u]; px_lb = px[l]; dx_lb = dx[l]
  #
  if (pfd_tp != 3 && pfd_tp != 4 && pfd_tp != 5) {
    stop("pfd_tp must be one of 3, 4, 5")
  }
  return(refine(UC, lb, ub, px_lb, dx_lb, pfd, pfd_tp))
}

#' @rdname bound_refine
#' @export
bound <- function(UC, px) {
  l = 1; u = length(px); # i = 0: not found
  #
  while (u - l > 1) {
    i = floor((l + u) / 2)
    if (px[i] < UC) { l = i }
    else if (px[i] > UC) { u = i }
    else { return(c(l, u, i)) }
  }
  return(c(l, u, 0))
}

#' @rdname bound_refine
#' @export
refine <- function(UC, lb, ub, px_lb, dx_lb, pfd, pfd_tp) {
  # Newton-Raphson Method
  x = lb - (px_lb - UC) / dx_lb
  x = ifelse(x < ub, x, (lb + ub) / 2) # make sure it lives within the bounds
  #
  if (pfd_tp == 3) dpearson = dpearson6
  else if (pfd_tp == 4) dpearson = dpearson8
  else if (pfd_tp == 5) dpearson = dpearson10
  else stop("pfd_tp must be one of 3, 4, 5")
  #
  while (ub - lb > 1e-10) {
    x0 = x
    dx = dpearson(x, pfd)
    px = px_lb + (x - lb) * (dx_lb + dx) / 2
    delta = px - UC
    # if (delta > 1e-5) {
    if (delta > 0) {
      ub = x
      x = x - delta / dx # when delta/P_X is very small, X stuck to ub
      x = ifelse(lb < x, x, (lb + ub) / 2)
    # } else if (delta < -1e-5) {
    } else if (delta < 0) {
      lb = x
      px_lb = px; dx_lb = dx
      x = x - delta / dx
      x = ifelse(x < ub, x, (lb + ub) / 2)
    } else { return(x) }
    #
    if (abs(x - x0) < 1e-10) { break }
  }
  return(x)
}
