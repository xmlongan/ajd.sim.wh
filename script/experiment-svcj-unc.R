library(ajd.sim.wh)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------SVCJ------------------------------------------------------
S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1; sigma_s = 0.0001
mu_v = 0.05; rhoJ = -0.38
# true_price = 6.8619 # true option price
#
cl = makeCluster(12)
registerDoParallel(cl)
#
N = 10000
# err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
#   N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v,
#   mu_b, rhoJ, sigma_s, true_price)
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v, mu_b, rhoJ, sigma_s)
write_to_csv(err_dur, "./script/unc-svcj-moms8-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v, mu_b, rhoJ, sigma_s)
write_to_csv(err_dur, "./script/unc-svcj-moms8-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v, mu_b, rhoJ, sigma_s)
write_to_csv(err_dur, "./script/unc-svcj-moms8-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v, mu_b, rhoJ, sigma_s)
write_to_csv(err_dur, "./script/unc-svcj-moms8-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v, mu_b, rhoJ, sigma_s)
write_to_csv(err_dur, "./script/unc-svcj-moms8-2560K.csv")
#
stopCluster(cl)
