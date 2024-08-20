#' Partial Fraction Decomposition
#'
#' @name PFDecomp
#' @description
#' Partial Fraction Decomposition of the following polynomial fraction:
#' \deqn{
#'   \frac{a+x}{c_0x^0 + \cdots + c_nx^n},
#' }
#' where \eqn{n=3, 4, 5}. The decomposition depends on the roots of the
#' polynomial in the denominator \eqn{c_0x^0 + \cdots + c_nx^n}.
#'
#' @param coef5 vector of the coefficients, \eqn{\{a, c_0, \cdots, c_3\}}
#' @param coef6 vector of the coefficients, \eqn{\{a, c_0, \cdots, c_4\}}
#' @param coef7 vector of the coefficients, \eqn{\{a, c_0, \cdots, c_5\}}
#'
#' @return solution of the partial fraction decomposition:
#' - `PFDecomp3()`: list of x1 (real),    x2 (complex), A1,     A2, B2
#' - `PFDecomp4()`: list of x1 (complex), x2 (complex), A1, B1, A2, B2
#' - `PFDecomp5()`: list of x1 (real),    x2 (complex), x3 (complex),
#'   A1, A2, B2, A3, B3
#' @export
#'
#' @examples
#' coef5 = c(-0.08538000, 0.01466260, -0.07073738, 0.10918559, 0.05940828)
#' PFDecomp3(coef5)
PFDecomp3 <- function(coef5) {
  a = coef5[1]; c0 = coef5[2]; c1 = coef5[3]; c2 = coef5[4]; c3 = coef5[5]
  z = polyroot(coef5[2:5])
  # the first two roots are the conjugate complex roots
  # - real root z[3]
  # - a pair of conjugate complex roots: z[1], z[2]
  x1 = Re(z[3]) #: label the last as x1
  x2 = z[1]     #: label the first as x2
  #
  A1 = (1/c3) * (x1+a)/(x1^2 + (Mod(x2))^2 - 2*x1*Re(x2))
  A2 = -A1
  B2 = (A1 * (Mod(x2)^2) - a/c3)/x1
  #
  pfd = list(x1=x1, x2=x2, A1=A1, A2=A2, B2=B2)
  return(pfd)
}

#' @rdname PFDecomp
PFDecomp4 <- function(coef6) {
  a  = coef6[1]; c0 = coef6[2]; c1 = coef6[3];
  c2 = coef6[4]; c3 = coef6[5]; c4 = coef6[6]
  #
  z = polyroot(coef6[2:6])
  #
  # two distinct pairs of conjugate complex roots
  # - one pair    : z[1], z[2]
  # - another pair: z[3], z[4]
  #
  p1 = -2 * Re(z[1]); q1 = (Mod(z[1]))^2
  p2 = -2 * Re(z[3]); q2 = (Mod(z[3]))^2
  #
  A = c(
    1, p2, q2, 0,
    0,  1, p2, q2,
    1, p1, q1, 0,
    0,  1, p1, q1
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = FALSE)
  b = c(0, 0, 1/c4, a/c4)
  #
  A1B1_A2B2 = solve(A, b)
  pfd = list(x1 = z[1], x2 = z[3],
             A1 = A1B1_A2B2[1], B1 = A1B1_A2B2[2],
             A2 = A1B1_A2B2[3], B2 = A1B1_A2B2[4])
  return(pfd)
}

#' @rdname PFDecomp
PFDecomp5 <- function(coef7) {
  a  = coef7[1]; c0 = coef7[2]; c1 = coef7[3]; c2 = coef7[4]
  c3 = coef7[5]; c4 = coef7[6]; c5 = coef7[7]
  #
  # one real root and two distinct pairs of conjugate complex roots:
  #    - one real root:                       x[3]
  #    - one pair of conjugate complex roots: x[0], x[1]
  #    - another pair:                        x[2], x[4]
  #
  z = polyroot(coef7[2:7])
  #
  x1 = Re(z[4])
  p2 = -2 * Re(z[1]); q2 = (Mod(z[1]))^2
  p3 = -2 * Re(z[3]); q3 = (Mod(z[3]))^2
  #
  A = c(
    1,                  1,       0,        1,        0,
    p2+p3,          p3-x1,       1,    p2-x1,        1,
    q3+p2*p3+q2, q3-x1*p3,    p3-x1, q2-x1*p2,    p2-x1,
    p2*q3+q2*p3,   -x1*q3, q3-x1*p3,   -x1*q2, q2-x1*p2,
    q2*q3,              0,   -x1*q3,        0,   -x1*q2
  )
  A = matrix(A, nrow = 5, ncol = 5, byrow = FALSE)
  b = c(0, 0, 0, 1/c5, a/c5)
  A1_A2B2_A3B3 = solve(A, b)
  #
  pfd = list(x1 = x1, x2 = z[1], x3 = z[3],
             A1 = A1_A2B2_A3B3[1],
             A2 = A1_A2B2_A3B3[2], B2 = A1_A2B2_A3B3[3],
             A3 = A1_A2B2_A3B3[4], B3 = A1_A2B2_A3B3[5])
  return(pfd)
}
