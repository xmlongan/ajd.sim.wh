library(ajd.sim.wh)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")
#---------------------SVCJ------------------------------------------------------
S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1; sigma_s = 0.0001
mu_v = 0.05; rhoJ = -0.38
true_price = 6.8619 # true option price
#
mu = r - lmbd * mu_b
mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
h = 1

par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
           mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h)
moms = rep(0, 8) # centralized variable whose first moment = 0
for (i in 2:8) {moms[i] = eval_mom_svcj(fmu.svcj[[i]], par)}
#
# err = rep(0, 200); dur = rep(0, 200)
# for (g in 1:200) {
#   err.dur = ajd.sim.wh::price_svcj(
#     N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v,
#     mu_b, rhoJ, sigma_s, true_price)
#   err[g] = err.dur[1]
#   dur[g] = err.dur[2]
# }
# df = data.frame(err = err, dur = dur)
# write.csv(df, file = "./script/svcj-moms8-0010K.csv")
#
cl = makeCluster(12)
registerDoParallel(cl)
#
N = 10000
# err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
#   N, S, K, v0, tau, r, k, theta, sigma, rho, lmbd, mu_v,
#   mu_b, rhoJ, sigma_s, true_price)
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s, moms,
  true_price)
write_to_csv(err_dur, "./script/svcj-moms8-0010K.csv")
#
N = 40000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s, moms,
  true_price)
write_to_csv(err_dur, "./script/svcj-moms8-0040K.csv")
#
N = 160000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s, moms,
  true_price)
write_to_csv(err_dur, "./script/svcj-moms8-0160K.csv")
#
N = 640000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s, moms,
  true_price)
write_to_csv(err_dur, "./script/svcj-moms8-0640K.csv")
#
N = 2560000
err_dur = foreach(g = 1:200) %dopar% ajd.sim.wh::price_svcj(
  N, S, K, v0, tau, r, k, theta, lmbd, mu_v, mu_b, rhoJ, sigma_s, moms,
  true_price)
write_to_csv(err_dur, "./script/svcj-moms8-2560K.csv")
#
stopCluster(cl)
