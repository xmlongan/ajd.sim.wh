require(devtools)
load_all()

S = 100; K = 100; v0 = 0.09; k = 2.00; theta = 0.09; sigma = 1.00
rho = -0.3; r = 0.05; tau = 5
true_price = 34.9998

par_hest = list(v0=v0, k=k, theta=theta, sigma=sigma, rho=rho, h=tau)
moms = rep(0, 8) # centralized variable, whose first moment always = 0
for (i in 2:8) {moms[i] = eval_mom_hest(fmu.hest[[i]], par_hest)}

N = 160000
G = 20
err = rep(0, G)
dur = rep(0, G)

for (g in 1:G) {
  err_dur = price_hest(N, S, K, v0, tau, r, k, theta, moms, true_price)
  err[g] = err_dur[1]
  dur[g] = err_dur[2]
  # printf("error = % f, cptm = %.2f\n", err[g], dur[g])
}
printf("RMSE = %f, cptm = %f\n", sqrt(mean(err[err != 0]^2)), mean(dur[dur != 0]))
hist(err[err!=0])
hist(dur[dur!=0])

