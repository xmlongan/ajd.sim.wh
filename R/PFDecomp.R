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
#' - `PFDecomp3()`: list of the roots and decomposition results, four types,
#'   i.e., `31`, `32`, `33`, `34`
#' - `PFDecomp4()`: list of the roots and decomposition results, nine types,
#'   i.e., `41`, `42`, ..., `49`
#' - `PFDecomp5()`: list of the roots and decomposition results, twelve types,
#'   i.e., `51`, `52`, ..., `512`, not fully implemented.
#' @export
#'
#' @examples
#' coef5 = c(-0.08538000, 0.01466260, -0.07073738, 0.10918559, 0.05940828)
#' PFDecomp3(coef5)
PFDecomp3 <- function(coef5) {
  a = coef5[1]; c3 = coef5[5]; z = polyroot(coef5[2:5])
  zlist = classify3(z)
  #
  if (zlist$type == 31) { return(pfd31(zlist$ordered_z, a, c3)) }
  if (zlist$type == 32) { return(pfd32(zlist$ordered_z, a, c3)) }
  if (zlist$type == 33) { return(pfd33(zlist$ordered_z, a, c3)) }
  if (zlist$type == 34) { return(pfd34(zlist$ordered_z, a, c3)) }
  #
  cat("coef5 = "); cat(coef5, sep = ','); cat('\n')
  cat("polyroot(z) = \n")
  cat('\t'); cat(z[1]); cat('\n'); cat('\t'); cat(z[2]); cat('\n')
  cat('\t'); cat(z[3]); cat('\n')
  warning("PFDecomp3() not implemented for this case")
}

classify3 <- function(z, eps = 1e-10) {
  I = abs(Im(z)) < eps; rroot = z[I]; croot = z[!I]
  #
  if (length(croot) == 2) {
    # one pair of conjugates
    if (is.conj(z[1], z[2])) {
      zlist = list(ordered_z = z[c(3,1,2)], type = 31)
      # real first, then complex
    } else if (is.conj(z[1], z[3])) {
      zlist = list(ordered_z = z[c(2,1,3)], type = 31)
    } else {
      zlist = list(ordered_z = z[c(1,2,3)], type = 31)
    }
  } else {
    # all real roots
    rt = sort(Re(rroot)) # rt[1] <= rt[2] <= rt[3]
    if (abs(rt[1] - rt[2]) < eps) {
      if (abs(rt[2] - rt[3]) < eps) {
        # x1 = x2 = x3
        ordered_z = rt; type = 32
      } else {
        # x1 = x2 != x3
        ordered_z = rt; type = 33
      }
    } else {
      if (abs(rt[2] - rt[3]) < eps) {
        # x1 != x2 = x3
        ordered_z = rt[c(2,3,1)]; type = 33
      } else {
        # x1 != x2 != x3
        ordered_z = rt; type = 34
      }
    }
    zlist = list(ordered_z = ordered_z + 0i, type = type) # keep consistent
  }
  return(zlist)
}

pfd31 <- function(z, a, c3) {
  # x1: real, (x2, x3): conjugate
  x1 = Re(z[1])
  p = -2 * Re(z[2]); q = (Mod(z[2]))^2
  A1 = (x1 + a) / (c3 * (x1^2 + p*x1 + q))
  A2 = -A1
  B2 = (A1 * q - a / c3) / x1
  pfd = list(x1 = x1, A1 = A1,
             x2 = z[2], A2 = A2, B2 = B2,
             type = 31)
  return(pfd)
}

pfd32 <- function(z, a, c3) {
  # x1 = x2 = x3: real
  x1 = Re(z[1])
  A1 = 0
  A2 = 1/c3
  A3 = (a + x1) / c3
  pfd = list(x1 = x1, A1 = 0, A2 = A2, A3 = A3,
             type = 32)
  return(pfd)
}

pfd33 <- function(z, a, c3) {
  # x1 = x2 != x3: real
  x1 = Re(z[1]); x3 = Re(z[3])
  A1 = - (a + x3) / (c3 * (x3 - x1)^2)
  A2 = 1 / c3 + (a + x3) / (c3 * (x1- x3))
  A3 = -A1 # (a + x3) / (c3 * (x3 - x1)^2)
  pfd = list(x1 = x1, A1 = A1, A2 = A2,
             x3 = x3, A3 = A3,
             type = 33)
  return(pfd)
}

pfd34 <- function(z, a, c3) {
  # x1 != x2 != x3: real
  x1 = Re(z[1]); x2 = Re(z[2]); x3 = Re(z[3])
  term = (x1^2 + x1*x2 + a*(x1+x2)) / (x1^2 - x2^2)
  A1 = term /  (c3 * (x1 - x3))
  A2 = (1 - term) / (c3 * (x2 - x3))
  A3 = -(A1 + A2)
  pfd = list(x1 = x1, A1 = A1,
             x2 = x2, A2 = A2,
             x3 = x3, A3 = A3,
             type = 34)
  return(pfd)
}


#' @rdname PFDecomp
#' @export
PFDecomp4 <- function(coef6) {
  a = coef6[1]; c4 = coef6[6]; z = polyroot(coef6[2:6])
  zlist = classify(z)
  #
  if (zlist$type == 41) { return(pfd41(zlist$ordered_z, a, c4)) }
  if (zlist$type == 42) { return(pfd42(zlist$ordered_z, a, c4)) }
  if (zlist$type == 43) { return(pfd43(zlist$ordered_z, a, c4)) }
  if (zlist$type == 44) { return(pfd44(zlist$ordered_z, a, c4)) }
  if (zlist$type == 45) { return(pfd45(zlist$ordered_z, a, c4)) }
  if (zlist$type == 46) { return(pfd46(zlist$ordered_z, a, c4)) }
  if (zlist$type == 47) { return(pfd47(zlist$ordered_z, a, c4)) }
  if (zlist$type == 48) { return(pfd48(zlist$ordered_z, a, c4)) }
  if (zlist$type == 49) { return(pfd49(zlist$ordered_z, a, c4)) }
  #
  cat("coef6 = "); cat(coef6, sep = ','); cat('\n')
  cat("polyroot(z) = \n")
  cat('\t'); cat(z[1]); cat('\n'); cat('\t'); cat(z[2]); cat('\n')
  cat('\t'); cat(z[3]); cat('\n'); cat('\t'); cat(z[4]); cat('\n')
  warning("PFDecomp4() not implemented for this case")
}

classify <- function(z, eps=1e-10) {
  I = abs(Im(z)) < eps; rroot = z[I]; croot = z[!I]
  #
  if (length(croot) == 4) {
    # two pairs of conjugates
    if (is.conj(z[1], z[2])) {
      if (is.equal(z[1], z[3]) || is.equal(z[1], z[4])) {
        zlist = list(ordered_z = z, type = 41)
      } else {
        zlist = list(ordered_z = z, type = 42)
      }
    } else if (is.conj(z[1], z[3])) {
      if (is.equal(z[1], z[2]) || is.equal(z[1], z[4])) {
        zlist = list(ordered_z = z[c(1,3,2,4)], type = 41)
      } else {
        zlist = list(ordered_z = z[c(1,3,2,4)], type = 42)
      }
    } else {
      if (is.equal(z[1], z[2]) || is.equal(z[1], z[3])) {
        zlist = list(ordered_z = z[c(1,4,2,3)], type = 41)
      } else {
        zlist = list(ordered_z = z[c(1,4,2,3)], type = 42)
      }
    }
  } else if (length(croot) == 2) {
    # one pair of conjugates
    if (abs(Re(rroot[1]) - Re(rroot[2])) < eps) {
      zlist = list(ordered_z = c(rroot, croot), type = 43)
    } else {
      if (Re(rroot[1]) < Re(rroot[2])) {
        zlist = list(ordered_z = c(rroot, croot), type = 44)
      } else {
        zlist = list(ordered_z = c(rroot[2], rroot[1], croot), type = 44)
      }
    }
  } else {
    # all real roots
    rt = sort(Re(rroot)) # rt[1] <= rt[2] <= rt[3] <= rt[4]
    if (abs(rt[1] - rt[2]) < eps) {
      if (abs(rt[2] - rt[3]) < eps) {
        if (abs(rt[3] - rt[4]) < eps) {
          # x1 = x2 = x3 = x4
          ordered_z = rt; type = 45
        } else {
          # x1 = x2 = x3 != x4
          ordered_z = rt; type = 46
        }
      } else {
        if (abs(rt[3] - rt[4]) < eps) {
          # x1 = x2 != x3 = x4
          ordered_z = rt; type = 47
        } else {
          # x1 = x2 != x3 != x4
          ordered_z = rt; type = 48
        }
      }
    } else {
      if (abs(rt[2] - rt[3]) < eps) {
        if (abs(rt[3] - rt[4]) < eps) {
          # x1 != x2 = x3 = x4
          ordered_z = rt[c(2,3,4,1)]; type = 46
        } else {
          # x1 != x2 = x3 != x4
          ordered_z = rt[c(2,3,1,4)]; type = 48
        }
      } else {
        if (abs(rt[3] - rt[4]) < eps) {
          # x1 != x2 != x3 = x4
          ordered_z = rt[c(3,4,1,2)]; type = 48
        } else {
          # x1 != x2 != x3 != x4
          ordered_z = rt; type = 49
        }
      }
    }
    zlist = list(ordered_z = ordered_z + 0i, type = type) # keep consistent
  }
  return(zlist)
}

pfd41 <- function(z, a, c4) {
  # all complex: (x1, x2) = (x3, x4)
  pfd = list(x1 = z[1], A1 = 1/c4, B1 = a/c4, type = 41)
  return(pfd)
}

pfd42 <- function(z, a, c4) {
  # all complex: (x1, x2) != (x3, x4)
  p1 = -2 * Re(z[1]); q1 = (Mod(z[1]))^2
  p2 = -2 * Re(z[3]); q2 = (Mod(z[3]))^2
  #
  # cat(sprintf("p1 = %f, q1 = %f\n", p1, q1))
  # cat(sprintf("p2 = %f, q2 = %f\n", p2, q2))
  #
  A = c(
    1, 0, 1, 0,
    p2,1, p1,1,
    q2,p2,q1,p1,
    0, q2,0,q1
  )
  A = matrix(A, nrow=4, ncol=4, byrow=TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1B1_A2B2 = solve(A, b)
  pfd = list(x1 = z[1], A1 = A1B1_A2B2[1], B1 = A1B1_A2B2[2],
             x2 = z[3], A2 = A1B1_A2B2[3], B2 = A1B1_A2B2[4],
             type = 42)
  return(pfd)
}

pfd43 <- function(z, a, c4) {
  # 2 real, 2 complex: x1 = x2, (x3, x4)
  x1 = Re(z[1])
  p = -2 * Re(z[3]); q = (Mod(z[3]))^2
  #
  A = c(
    1, 0, 1, 0,
    p-x1, 1, -2*x1, 1,
    q-x1*p, p, x1^2, -2*x1,
    -x1*q, q, 0, x1^2
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2_A3B3 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2_A3B3[1], A2 = A1A2_A3B3[2],
             x3 = z[3], A3 = A1A2_A3B3[3], B3 = A1A2_A3B3[4],
             type = 43)
  return(pfd)
}

pfd44 <- function(z, a, c4) {
  # 2 real, 2 complex: x1 != x2, (x3, x4)
  x1 = Re(z[1]); x2 = Re(z[2])
  #
  x3 = z[3]
  p = -2 * Re(x3); q = (Mod(x3))^2
  #
  A = c(
    1, 1, 1, 0,
    p-x2, p-x1, -(x1+x2), 1,
    q-x2*p, q-x1*p, x1*x2, -(x1+x2),
    -x2*q, -x1*q, 0, x1*x2
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2_A3B3 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2_A3B3[1],
             x2 = x2, A2 = A1A2_A3B3[2],
             x3 = x3, A3 = A1A2_A3B3[3], B3 = A1A2_A3B3[4],
             type = 44)
  return(pfd)
}

pfd45 <- function(z, a, c4) {
  # 4 real: x1 = x2 = x3 = x4
  x1 = Re(z[1])
  pfd = list(x1 = x1, A3 = 1/c4, A4 = a/c4 + x1/c4, type = 45)
  return(pfd)
}

pfd46 <- function(z, a, c4) {
  # 4 real: x1 = x2 = x3 != x4
  x1 = Re(z[1]); x4 = Re(z[4])
  A = c(
    1, 0, 0, 1,
    -(2*x1+x4), 1, 0, -3*x1,
    x1^2+2*x1*x4, -2*(x1+x4), 1, 3*x1^2,
    -x1^2*x4, x1*x4, -x4, -x1^3
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2A3A4 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2A3A4[1], A2 = A1A2A3A4[2], A3 = A1A2A3A4[3],
             x4 = x4, A4 = A1A2A3A4[4],
             type = 46)
  return(pfd)
}

pfd47 <- function(z, a, c4) {
  # 4 real: x1 = x2 != x3 = x4
  x1 = Re(z[1]); x3 = Re(z[3])
  A = c(
    1, 0, 1, 0,
    -(x1+2*x3), 1, -(x3+2*x1), 1,
    2*x1*x3+x3^2, -2*x3, 2*x3*x1+x1^2, -2*x1,
    -x1*x3^2, x3^2, -x3*x1^2, x1^2
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2A3A4 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2A3A4[1], A2 = A1A2A3A4[2],
             x3 = x3, A3 = A1A2A3A4[3], A4 = A1A2A3A4[4],
             type = 47)
  return(pfd)
}

pfd48 <- function(z, a, c4) {
  # 4 real: x1 = x2 != x3 != x4
  x1 = Re(z[1]); x3 = Re(z[3]); x4 = Re(z[4])
  A = c(
    1, 0, 1, 1,
    -(x1+x3+x4), 1, -(2*x1+x4), -(2*x1+x3),
    x3*x4+x1*(x3+x4), -(x3+x4), x1^2+2*x1*x4, x1^2+2*x1*x3,
    -x1*x3*x4, x3*x4, -x1^2*x4, -x1^2*x3
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2A3A4 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2A3A4[1], A2 = A1A2A3A4[2],
             x3 = x3, A3 = A1A2A3A4[3],
             x4 = x4, A4 = A1A2A3A4[4],
             type = 48)
  return(pfd)
}

pfd49 <- function(z, a, c4) {
  # 4 real: x1 != x2 != x3 != x4
  x1 = Re(z[1]); x2 = Re(z[2]); x3 = Re(z[3]); x4 = Re(z[4])
  A = c(
    1, 1, 1, 1,
    -(x2+x3+x4), -(x1+x3+x4), -(x1+x2+x4), -(x1+x2+x3),
    x3*x4 + x2*(x3+x4), x3*x4 + x1*(x3+x4), x1*x2 + x4*(x1+x2), x1*x2 + x3*(x1+x2),
    -x2*x3*x4, -x1*x3*x4, -x1*x2*x4, -x1*x2*x3
  )
  A = matrix(A, nrow = 4, ncol = 4, byrow = TRUE)
  b = c(0, 0, 1/c4, a/c4)
  A1A2A3A4 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1A2A3A4[1],
             x2 = x2, A2 = A1A2A3A4[2],
             x3 = x3, A3 = A1A2A3A4[3],
             x4 = x4, A4 = A1A2A3A4[4],
             type = 49)
  return(pfd)
}


#' @rdname PFDecomp
PFDecomp5 <- function(coef7) {
  a = coef7[1]; c0 = coef7[2]; c1 = coef7[3]; c2 = coef7[4]
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
    1, 1, 0, 1, 0,
    p2 + p3, p3 - x1, 1, p2 - x1, 1,
    q3 + p2 * p3 + q2, q3 - x1 * p3, p3 - x1, q2 - x1 * p2, p2 - x1,
    p2 * q3 + q2 * p3, -x1 * q3, q3 - x1 * p3, -x1 * q2, q2 - x1 * p2,
    q2 * q3, 0, -x1 * q3, 0, -x1 * q2
  )
  A = matrix(A, nrow = 5, ncol = 5, byrow = FALSE)
  b = c(0, 0, 0, 1 / c5, a / c5)
  A1_A2B2_A3B3 = solve(A, b)
  #
  pfd = list(x1 = x1, A1 = A1_A2B2_A3B3[1],
             x2 = z[1], A2 = A1_A2B2_A3B3[2], B2 = A1_A2B2_A3B3[3],
             x3 = z[3], A3 = A1_A2B2_A3B3[4], B3 = A1_A2B2_A3B3[5])
  return(pfd)
}
