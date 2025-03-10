library(ajd.sim.wh)
library(foreach)
library(doParallel)
source("./script/write_to_csv.R")

S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
rho = -0.82; r = 0.0319; tau = 1; lmbd = 0.47; mu_b = -0.1
sigma_s = 0.0001; mu_v = 0.05; rhoJ = -0.38; true_price = 6.8619

mu = r - lmbd * mu_b
mu_s = log((1 + mu_b) * (1 - rhoJ*mu_v)) - sigma_s^2 / 2
h = 1

par = list(v0=v0, mu=mu, k=k, theta=theta, sigma=sigma, rho=rho, lmbd=lmbd,
           mu_v=mu_v, rhoJ=rhoJ, mu_s=mu_s, sigma_s=sigma_s, h=h,
           S=S, K=K, r=r)

cl = makeCluster(12)
registerDoParallel(cl)
#
G = 4^4
N = 10000
# calculate the unconditional call option price
# for measuring the accuracy of the rpearson(N, moms) when unconditional
# SVCJ moments supplied
err_dur = foreach(g = 1:G) %dopar% ajd.sim.wh::svcj_cprice_unc(par, N)
write_to_csv(err_dur, "./script/price_time.csv")
