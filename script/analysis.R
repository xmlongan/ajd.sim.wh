library(readr)

# heston-s1-WH
No_simulation = c(1, 4, 16, 64, 256) * 10000
txts = c("0010","0040","0160","0640","2560")
fnames = paste0("script/results/heston-s1-moms8-", txts, "K.csv")
RMSE = rep(0, length(txts))
cptm = rep(0, length(txts))
for (i in 1:length(fnames)) {
  error_dur <- read_csv(fnames[i])
  error = error_dur$error
  ctime = error_dur$computing_time
  RMSE[i] = sqrt(mean(error^2))
  cptm[i] = mean(ctime)
}
df_h1 = data.frame(No_simulation=No_simulation, RMSE=RMSE, Computing_time=cptm)
write.csv(df_h1, file="script/df_h1_wh.csv")

# heston-s2-WH
No_simulation = c(1, 4, 16, 64, 256) * 10000
txts = c("0010","0040","0160","0640","2560")
fnames = paste0("script/results/G-200/heston-s2-2/heston-s2-moms8-", txts, "K.csv")
# fnames = paste0("script/results/heston-s2-moms8-", txts, "K.csv")
RMSE = rep(0, length(txts))
cptm = rep(0, length(txts))
for (i in 1:length(fnames)) {
  error_dur <- read_csv(fnames[i])
  error = error_dur$error
  ctime = error_dur$computing_time
  RMSE[i] = sqrt(mean(error^2))
  cptm[i] = mean(ctime)
}
df_h2 = data.frame(No_simulation=No_simulation, RMSE=RMSE, Computing_time=cptm)
write.csv(df_h2, file="script/df_h2_wh.csv")

# SVJ-WH
# No_simulation = c(4, 16, 64, 256) * 10000
# txts = c("0040","0160","0640","2560")
# fnames = paste0("script/results/SVJ-BK-", txts, "K.csv")
No_simulation = c(1, 4, 16, 64, 256) * 10000
txts = c("0010","0040","0160","0640", "2560")
fnames = paste0("script/svj-moms8-", txts, "K.csv")
RMSE = rep(0, length(txts))
cptm = rep(0, length(txts))
for (i in 1:length(fnames)) {
  error_dur <- read_csv(fnames[i])
  error = error_dur$error
  ctime = error_dur$computing_time
  RMSE[i] = sqrt(mean(error^2))
  cptm[i] = mean(ctime)
}
df_svj = data.frame(No_simulation=No_simulation, RMSE=RMSE, Computing_time=cptm)
write.csv(df_svj, file="script/df_svj_wh.csv")

# SVCJ-WH
No_simulation = c(1, 4, 16, 64, 256) * 10000
txts = c("0010","0040","0160","0640","2560")
fnames = paste0("script/results/svcj-moms8-", txts, "K.csv")
RMSE = rep(0, length(txts))
cptm = rep(0, length(txts))
for (i in 1:length(fnames)) {
  error_dur <- read_csv(fnames[i])
  error = error_dur$error
  ctime = error_dur$computing_time
  RMSE[i] = sqrt(mean(error^2))
  cptm[i] = mean(ctime)
}
df_svcj = data.frame(No_simulation=No_simulation, RMSE=RMSE, Computing_time=cptm)
write.csv(df_svcj, file="script/df_svcj_wh.csv")

# SVCJ-WH steady state distribution
No_simulation = c(1, 4, 16, 64, 256) * 10000
txts = c("0010","0040","0160","0640","2560")
fnames = paste0("script/unc-svcj-moms8-", txts, "K.csv")
RMSE = rep(0, length(txts))
cptm = rep(0, length(txts))
for (i in 1:length(fnames)) {
  error_dur <- read_csv(fnames[i])
  error = error_dur$error - 7.109439 #7.103647
  ctime = error_dur$computing_time
  RMSE[i] = sqrt(mean(error^2))
  cptm[i] = mean(ctime)
}
df_svcj_unc = data.frame(No_simulation=No_simulation, RMSE=RMSE, Computing_time=cptm)
write.csv(df_svcj_unc, file="script/df_svcj_wh_unc.csv")
