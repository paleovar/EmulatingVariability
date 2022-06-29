dir.runs = "output-current/detrended_both_850_1850_emics/"
files = list.files(dir.runs, pattern = "ME.*_v")

for (file in files) {
  res = loadRData(paste0(dir.runs, file))
  means = get_post_means(res)
  print(paste0("lambda: ", round(means$lambda, 3), ", weights: ", round(means$weights, 3)))
  print(res$diagnostics$gelman_diag_mpsrf)
}

res1 = loadRData("~/maybritt-msc/EBMSpectrum/output-current/detrended_both_850_1850_emics/emic_DC_2box_v7.rda")
res2 = loadRData("~/maybritt-msc/EBMSpectrum/output-current/detrended_both_850_1850_emics/emic_DC_2box_v3.rda")
plot_two_fits(res1, res2)
