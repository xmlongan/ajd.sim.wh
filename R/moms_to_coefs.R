#' Pearson Coefficients from Moments
#'
#' @description
#' Compute the coefficients in the following by matching the first 6|8|10
#' moments:
#' \deqn{
#'   \frac{d p(x)} {d x} = - \frac{a+x}{c_0x^0 + \cdots + c_nx^n} p(x)
#' }
#' where \eqn{n=3, 4, 5}, \eqn{p(x)} is the (un-normalized) density of the
#' Pearson Distribution.
#'
#' @param moms vector of the first six|eight|ten moments
#'
#' @return vector of the coefficients:
#' - \eqn{\{a, c_0,\cdots, c_3\}} corresponding to the first 6 moments.
#' - \eqn{\{a, c_0,\cdots, c_4\}} corresponding to the first 8 moments.
#' - \eqn{\{a, c_0,\cdots, c_5\}} corresponding to the first 10 moments.
#' @export
#'
#' @examples
#' moms6 = c(0.02310, 0.01915, -0.00206, 0.00199, -0.00097, 0.00077)
#' moms8 = c(0.0231070307, 0.0191518038, -0.0020677784, 0.0019956050,
#'          -0.0009780987, 0.0007711325, -0.0006488957, 0.0006454927)
#' moms10 = c(0.0231070307, 0.0191518038, -0.0020677784, 0.0019956050,
#'           -0.0009780987, 0.0007711325, -0.0006488957, 0.0006454927,
#'           -0.0007204021, 0.0008985861)
#' coef5 = moms_to_coefs(moms6)
#' coef6 = moms_to_coefs(moms8)
#' coef7 = moms_to_coefs(moms10)
moms_to_coefs <- function(moms) {
  N = length(moms)
  if (N < 6) stop("at least match the first 6 moments")
  m1 = moms[1]; m2 = moms[2]; m3 = moms[3]
  m4 = moms[4]; m5 = moms[5]; m6 = moms[6]
  if (N >= 8)  { m7 = moms[7]; m8 = moms[8] }
  if (N >= 10) { m9 = moms[9]; m10 = moms[10] }
  #
  if (N == 6) {
    A = c(
          1,    m1,    m2,    m3,    m4,
          0,    -1, -2*m1, -3*m2, -4*m3,
         -1, -2*m1, -3*m2, -4*m3, -5*m4,
      -2*m1, -3*m2, -4*m3, -5*m4, -6*m5,
      -3*m2, -4*m3, -5*m4, -6*m5, -7*m6
    )
    A = matrix(A, nrow = 5, ncol = 5, byrow = FALSE)
    b = moms[1:5]
  } else if (N == 8) {
    A = c(
          1,    m1,    m2,    m3,    m4,    m5,
          0,    -1, -2*m1, -3*m2, -4*m3, -5*m4,
         -1, -2*m1, -3*m2, -4*m3, -5*m4, -6*m5,
      -2*m1, -3*m2, -4*m3, -5*m4, -6*m5, -7*m6,
      -3*m2, -4*m3, -5*m4, -6*m5, -7*m6, -8*m7,
      -4*m3, -5*m4, -6*m5, -7*m6, -8*m7, -9*m8
    )
    A = matrix(A, nrow = 6, ncol = 6, byrow = FALSE)
    b = moms[1:6]
  } else if (N == 10) {
    A = c(
      1,      0,    -1, -2*m1, -3*m2, -4*m3, -5*m4,
      m1,    -1, -2*m1, -3*m2, -4*m3, -5*m4, -6*m5,
      m2, -2*m1, -3*m2, -4*m3, -5*m4, -6*m5, -7*m6,
      m3, -3*m2, -4*m3, -5*m4, -6*m5, -7*m6, -8*m7,
      m4, -4*m3, -5*m4, -6*m5, -7*m6, -8*m7, -9*m8,
      m5, -5*m4, -6*m5, -7*m6, -8*m7, -9*m8, -10*m9,
      m6, -6*m5, -7*m6, -8*m7, -9*m8, 10*m9, -11*m10
    )
    A = matrix(A, nrow = 7, ncol = 7, byrow = FALSE)
    b = moms[1:7]
  } else { stop("either match the first 6|8|10 moments, no other options") }
  x = solve(A, -b)
  return(x)
}
