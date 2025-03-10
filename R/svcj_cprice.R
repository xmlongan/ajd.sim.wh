#' Evaluate Mean of the SVCJ Model
#'
#' @description
#' Evaluate mean of the SVCJ model
#'
#' @param par list of the parameters
#'
#' @return scalar
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38; true_price = 6.8619
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = tau
#' par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h,
#'            S=S, K=K, r=r)
#' ymean = mean_y(par)
mean_y <- function(par) {
  v0 = par$v0; mu = par$mu; k = par$k; theta = par$theta; sigma = par$sigma
  rho = par$rho; lmbd = par$lmbd; mu_v = par$mu_v; rhoJ = par$rhoJ
  mu_s = par$mu_s; sigma_s = par$sigma_s; tau = par$h

  beta = (1 - exp(-k * tau)) / (2 * k)
  #
  Ymean = lmbd * mu_v * beta / k + (rhoJ - 1/(2*k)) * lmbd * mu_v * tau
  Ymean = Ymean + lmbd * mu_s * tau + (mu - theta/2)*tau - beta*(v0 - theta)
  return(Ymean)
}

#' Call Option Pricing under the SVCJ Model
#'
#' @description
#' Pricing the European call option under the SVCJ model, using the numerical
#' integral with the Pearson density approximation of the unknown true density
#'
#' @param par list of the parameters
#' @param N number of samples to employ to estimate the unconditional call
#'  option price
#'
#' @return scalar
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38; true_price = 6.8619
#' mu = r - lmbd * mu_b
#' mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
#' h = tau
#' par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
#'            mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h,
#'            S=S, K=K, r=r)
#'
#' svcj_cprice(par) # normal call option pricing
#'
#' # For measuring errors of rpearson(N, moms) when unconditional moments
#' # supplied
#' N = 100
#' svcj_cprice_unc(par, N)
svcj_cprice <- function(par) {
  # call option price
  moms = rep(0, 8)
  for (i in 2:8) {moms[i] = ajd.sim.wh::eval_mom_svcj(ajd.sim.wh::fmu.svcj[[i]], par)}

  coefs = ajd.sim.wh::mom_to_coef(moms)
  pfd_tp = 4; dpearson = ajd.sim.wh::dpearson8
  pfd = ajd.sim.wh::PFDecomp4(coefs)
  #
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4])
  sd = sqrt(stdmoment[2]); skew = stdmoment[3]
  N = 10000
  if (skew < -1) {
    x = seq(-8*sd, 4*sd, length.out = N)
  } else {
    x = seq(-7*sd, 5*sd, length.out = N)
  }
  #
  dx = dpearson(x, pfd)
  h = x[2] - x[1]
  #
  cum = rep(0, N - 2)
  cum[1] = dx[2] * h
  for (i in 2:(N - 2)) {
    cum[i] = cum[i - 1] + dx[i + 1] * h # may too big number inside!
  }
  #
  px = rep(0, N)          # un-normalized cumulative probability
  px[2] = (dx[1] + dx[2]) * h / 2
  for (i in 3:N) {
    px[i] = (dx[1] + dx[i]) * h / 2 + cum[i - 2]
  }
  C = px[N] # constant
  dx = dx/C # normalize
  #
  ymean = mean_y(par)
  lb = log(par$K/par$S) - ymean
  lui = ajd.sim.wh::bound(lb, x)
  i = lui[2]
  if (lb < x[i-1] || lb > x[i]) {
    cat(sprintf("x[i] = %f, lb = %f, x[i-1] = %f\n", x[i], lb, x[i-1]))
  }
  if (i > N) {stop("i > N")}

  # lb, x[i]: \int_lb^x[i] dy
  lf = (dpearson(lb, pfd)/C) * (par$S * exp(lb   + ymean) - par$K)
  uf = dx[i] *                 (par$S * exp(x[i] + ymean) - par$K)
  val = ((lf + uf)/2) * (x[i] - lb)
  #
  if (i == N) {
    return(val)
  }
  for (n in i:(N-1)) {
    lf = dx[n] *   (par$S * exp(x[n]   + ymean) - par$K)
    uf = dx[n+1] * (par$S * exp(x[n+1] + ymean) - par$K)
    val = val + ((lf+uf)/2) * h
  }
  val = exp(-par$r * par$h) * val
  return(val)
}

#' @rdname svcj_cprice
#' @export
svcj_cprice_unc <- function(par, N) {
  ts = proc.time()
  # vs = rep(0, N)
  prices = rep(0, N)
  #
  for (n in 1:N) {
    # long period: 4 times long to reach steady state
    v = r_srjd(par$v0, par$k, par$theta, par$sigma, par$lmbd, par$mu_v, 4*par$h)
    par$v0 = v
    #
    # vs[n] = v
    prices[n] = svcj_cprice(par)
  }
  # return(list(Ys=Ys, vs=vs))
  te = proc.time()
  tt = te - ts
  return(c(mean(prices), tt[3]))
}
