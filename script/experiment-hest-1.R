library(ajd.sim.wh)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------Heston----------------------------------------------------
# parameter setting 1
S = 100; K = 100; v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61
rho = -0.7; r = 0.0319; tau = 1
true_price = 6.8061 # true option price

par_hest = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
moms = rep(0, 8) # centralized variable, whose first moment always = 0
for (i in 2:8) {moms[i] = eval_mom_hest(fmu.hest[[i]], par_hest)}

cl = makeCluster(10)
registerDoParallel(cl)
#
N = 10000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_hest(
  N, S, K, v0, tau, r, k, theta, moms, true_price)
write_to_csv(err_dur, "./script/heston-s1-moms8-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_hest(
  N, S, K, v0, tau, r, k, theta, moms, true_price)
write_to_csv(err_dur, "./script/heston-s1-moms8-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_hest(
  N, S, K, v0, tau, r, k, theta, moms, true_price)
write_to_csv(err_dur, "./script/heston-s1-moms8-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_hest(
  N, S, K, v0, tau, r, k, theta, moms, true_price)
write_to_csv(err_dur, "./script/heston-s1-moms8-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_hest(
  N, S, K, v0, tau, r, k, theta, moms, true_price)
write_to_csv(err_dur, "./script/heston-s1-moms8-2560K.csv")
#
stopCluster(cl)
