#' Checking Conjugate or Equal of Two Complex Numbers
#'
#' @name conjugate_or_equal
#' @param z1 complex number
#' @param z2 complex number
#' @param eps relative error
#'
#' @returns bool `TRUE` or `FALSE`
#'
#' @export
#'
#' @examples
#' z1 = 1+1.000001i
#' z2 = 1-1.000000i
#' is.conj(z1, z2)
is.conj <- function(z1, z2, eps = 1e-10) {
  # if (!is.complex(z1) || !is.complex(z2)) {
  #   stop("The two numbers must be complex!")
  # }
  #
  if (abs(Re(z1) - Re(z2)) < eps && abs(Im(z1) + Im(z2)) < eps) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' @rdname conjugate_or_equal
#' @export
is.equal <- function(z1, z2, eps = 1e-10) {
  if (abs(Re(z1) - Re(z2)) < eps && abs(Im(z1) - Im(z2)) < eps) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


adjust_lb_ub <- function(lb, ub, pfd) {
  if (pfd$type == 43 || pfd$type == 45) {
    # two real: x1 = x2 or 4 real: x1 = x2 = x3 = x4
    if (0 <= pfd$x1) {
      ub = ifelse(ub < pfd$x1, ub, pfd$x1 - 1e-5)
    } else {
      lb = ifelse(pfd$x1 < lb, lb, pfd$x1 + 1e-5)
    }
    lbub = c(lb, ub)
  } else if (pfd$type == 44) {
    # two real: x1 < x2
    lbub = adjust_lb_ub_2(lb, ub, pfd$x1, pfd$x2)
  } else if (pfd$type == 46) {
    # x1 = x2 = x3 != x4 | x2 = x3 = x4 != x1
    if (pfd$x1 < pfd$x4) {
      lbub = adjust_lb_ub_2(lb, ub, pfd$x1, pfd$x4)
    } else {
      lbub = adjust_lb_ub_2(lb, ub, pfd$x4, pfd$x1)
    }
  } else if (pfd$type == 47) {
    # x1 = x2 < x3 = x4
    lbub = adjust_lb_ub_2(lb, ub, pfd$x1, pfd$x3)
  } else if (pfd$type == 48) {
    # x1 = x2 != x3 != x4 | x2 = x3 != x1 != x4 | x3 = x4 != x1 != x2
    x1x2x3 = sort(pfd$x1, pfd$x3, pfd$x4)
    lbub = adjust_lb_ub_3(lb, ub, x1x2x3[1], x1x2x3[2], x1x2x3[3])
  } else if (pfd$type == 49) {
    # x1 != x2 != x3 != x4
    lbub = adjust_lb_ub_4(lb, ub, pfd$x1, pfd$x2, pfd$x3, pfd$x4)
  } else {
    lbub = c(lb, ub)
  }
  return(lbub)
}

adjust_lb_ub_2 <- function(lb, ub, x1, x2) {
  # x1 < x2
  if (0 <= x1) {
    ub = ifelse(ub < x1, ub, x1 - 1e-5)
  } else if (x2 <= 0) {
    lb = ifelse(x2 < lb, lb, x2 + 1e-5)
  } else {
    # x1 < 0 < x2
    lb = ifelse(x1 < lb, lb, x1 + 1e-5)
    ub = ifelse(ub < x2, ub, x2 - 1e-5)
  }
  return(c(lb, ub))
}

adjust_lb_ub_3 <- function(lb, ub, x1, x2, x3) {
  # x1 < x2 < x3
  if (0 <= x1) {
    ub = ifelse(ub < x1, ub, x1 - 1e-5)
  } else if (x3 <= 0) {
    lb = ifelse(x3 < lb, lb, x3 + 1e-5)
  } else {
    # x1 < 0 < x3
    if (x2 <= 0) {
      # x1 < x2 <= 0 < x3
      lb = ifelse(x2 < lb, lb, x2 + 1e-5)
      ub = ifelse(ub < x3, ub, x3 - 1e-5)
    } else {
      # x1 < 0 < x2 < x3
      lb = ifelse(x1 < lb, lb, x1 + 1e-5)
      ub = ifelse(ub < x2, ub, x2 - 1e-5)
    }
  }
  return(c(lb, ub))
}

adjust_lb_ub_4 <- function(lb, ub, x1, x2, x3, x4) {
  # x1 < x2 < x3 < x4
  if (0 <= x1) {
    ub = ifelse(ub < x1, ub, x1 - 1e-5)
  } else if (x4 <= 0) {
    lb = ifelse(x4 < lb, lb, x4 + 1e-5)
  } else {
    # x1 < 0 < x4
    if (0 <= x2) {
      # x1 < 0 <= x2 < x3 < x4
      lb = ifelse(x1 < lb, lb, x1 + 1e-5)
      ub = ifelse(ub < x2, ub, x2 - 1e-5)
    } else if (x3 <= 0) {
      # x1 < x2 < x3 <= 0 < x4
      lb = ifelse(x3 < lb, lb, x3 + 1e-5)
      ub = ifelse(ub < x4, ub, x4 - 1e-5)
    } else {
      # x1 < x2 < 0 < x3 < x4
      lb = ifelse(x2 < lb, lb, x2 + 1e-5)
      ub = ifelse(ub < x3, ub, x3 - 1e-5)
    }
  }
  return(c(lb, ub))
}

# adjust lb and ub for PFD3
adjust3_lb_ub <- function(lb, ub, pfd3) {
  if (pfd3$type == 31 || pfd3$type == 32) {
    # x1: real, (x2, x3): conjugate | x1 = x2 = x3: real
    if ( 0 <= pfd3$x1) {
      ub = ifelse(ub < pfd3$x1, ub, pfd3$x1 - 1e-5)
    } else {
      lb = ifelse(pfd3$x1 < lb, lb, pfd3$x1 + 1e-5)
    }
    lbub = c(lb, ub)
  } else if (pfd3$type == 33) {
    # x1 = x2 != x3: real | x2 = x3 != x1
    if (pfd3$x1 < pfd3$x3) {
      lbub = adjust_lb_ub_2(lb, ub, pfd3$x1, pfd3$x3)
    } else {
      lbub = adjust_lb_ub_2(lb, ub, pfd3$x3, pfd3$x1)
    }
  } else {
    # x1 < x2 < x3
    lbub = adjust_lb_ub_3(lb, ub, pfd3$x1, pfd3$x2, pfd3$x3)
  }
  return(lbub)
}
