library(ajd.sim.kbf)
# library(ajd.sim.wh)
library(devtools)
load_all()
library(PearsonDS)
source("script/analysis-Rmd/comp_C.R")
# Heston case 1
moms_1 = c(0.0000000000, 0.0186178689, -0.0033707270, 0.0022472255, -0.0012220384, 0.0009233124, -0.0007854778, 0.0007776901)
# Heston case 2
moms_2 = c(0.0000000, 0.5346568, -0.3868392, 1.6363802, -3.9943968, 16.5271025, -70.1374506, 365.1618607)

get_lb_ub <- function(moms) {
  # mean, var, skew, kurt
  stdmoment = ajd.sim.kbf::stdmom(moms[1:4])
  sd = sqrt(stdmoment[2]); skew = stdmoment[3]; kurt = stdmoment[4]
  #
  printf("sd = %f, skewness = %f, kurtosis = %f\n", sd, skew, kurt)
  #
  if (skew < 0) {         # left-tailed
    if (skew < -1) { lb = -8*sd; ub = 4*sd } else { lb = -7*sd; ub = 5*sd }
  } else if (skew > 0) {  # right-tailed
    if (skew > 1)  { lb = -4*sd; ub = 8*sd } else { lb = -5*sd; ub = 7*sd }
  } else {                # symmetric
    lb = -6*sd; ub = 6*sd
  }
  return(c(lb, ub))
}

shift <- function(moms) {
  lbub = get_lb_ub(moms)
  lb = lbub[1]; ub = lbub[2]

  coef6 = mom_to_coef(moms)
  pfd4 = PFDecomp4(coef6)
  lbub = adjust_lb_ub(lb, ub, pfd4)

  fmt = "before adjust: (lb,ub) = (%5f,%5f)\nafter adjust: (lb,ub) = (%5f,%5f)\n"
  cat(sprintf(fmt, lb, ub, lbub[1], lbub[2]))
  logC_type = comp_logC(lbub[1], lbub[2], pfd4)
  logC = logC_type$logC
  #
  # par(mfrow = c(2,1))
  #
  # x = seq(lbub[1], lbub[2], 0.001)
  x = seq(lb, ub, 0.001)
  d8 = dpearson8(x, pfd4) / exp(logC)
  #
  coef5 = mom_to_coef(moms[1:6])
  pfd3 = PFDecomp3(coef5)
  lbub = adjust3_lb_ub(lb, ub, pfd3)
  logC_type = comp_logC_6(lbub[1], lbub[2], pfd3)
  logC = logC_type$logC
  #
  d6 = dpearson6(x, pfd3) / exp(logC)
  d4 = PearsonDS::dpearson(x, moments = ajd.sim.kbf::stdmom(moms[1:4]))
  return(list(x=x, d4=d4, d6=d6, d8=d8))
}


par(mfrow = c(1, 2))
x_d4_d6_d8 = shift(moms_1)
x = x_d4_d6_d8$x
d4 = x_d4_d6_d8$d4
d6 = x_d4_d6_d8$d6
d8 = x_d4_d6_d8$d8
plot(x, d8,  type = 'l', lty = "solid", col='blue', ylab="density")
lines(x, d6, type = 'l', lty = "dashed", col='green4')
lines(x, d4, type = 'l', lty = "dotted", col='red')

cols = c("blue", "green4", "red")
ltys = c("solid", "dashed", "dotted")
lgd_txt = c("1st - 8th", "1st - 6th", "1st - 4th")
#
cols = rev(cols); ltys = rev(ltys); lgd_txt = rev(lgd_txt)
#
legend("topleft", inset = 0.02, title = "moments matched",
       legend = lgd_txt, col = cols, lty = ltys, cex=0.8)

x_d4_d6_d8 = shift(moms_2)
x = x_d4_d6_d8$x
d4 = x_d4_d6_d8$d4
d6 = x_d4_d6_d8$d6
d8 = x_d4_d6_d8$d8
plot(x, d8,  type = 'l', lty = "solid", col='blue', ylab="density")
lines(x, d6, type = 'l', lty = "dashed", col='green4')
lines(x, d4, type = 'l', lty = "dotted", col='red')

cols = c("blue", "green4", "red")
ltys = c("solid", "dashed", "dotted")
lgd_txt = c("1st - 8th", "1st - 6th", "1st - 4th")
#
cols = rev(cols); ltys = rev(ltys); lgd_txt = rev(lgd_txt)
#
legend("topleft", inset = 0.02, title = "moments matched",
       legend = lgd_txt, col = cols, lty = ltys, cex=0.8)
