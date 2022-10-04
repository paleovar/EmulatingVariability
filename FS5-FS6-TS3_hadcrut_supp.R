source("init.R")

# ------ 1 vs 2 vs 3 box ------------------

fit_1box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_1box.rda")
fit_2box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_2box.rda")
fit_3box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_3box.rda")

RMSE_1box <- sqrt(sum((fit_1box$input_params$y_obs - fit_1box$posteriors$model_fit$mean)^2) / length(fit_1box$input_params$y_obs))
RMSE_2box <- sqrt(sum((fit_2box$input_params$y_obs - fit_2box$posteriors$model_fit$mean)^2) / length(fit_2box$input_params$y_obs))
RMSE_3box <- sqrt(sum((fit_3box$input_params$y_obs - fit_3box$posteriors$model_fit$mean)^2) / length(fit_3box$input_params$y_obs))

data_3box = data.frame(year = 1850:2000,
                       temp = fit_3box$posteriors$model_fit$mean,
                       temp_lower = fit_3box$posteriors$model_fit$lower_quant,
                       temp_upper = fit_3box$posteriors$model_fit$upper_quant,
                       type = "3-box")

save <- F

lsize <- 0.7 #line size
tsize <- 12.5 #textsize

col1 = "deepskyblue2"
col2 = COL[["fit"]]
col3 = "black"
p1 <- ClimBayes::plot_two_fits(fit_1box, fit_2box, line_size=lsize, text_size=tsize, 
                               cols = c("darkgrey", "skyblue", COL[["fit"]]),
                               print_RMSE = F) +
  geom_line(aes(x = year, y = temp, col = type), data = data_3box, size = lsize) +
  geom_ribbon(data = data_3box,
              mapping = aes(x = year, ymin = temp_upper, ymax = temp_lower, fill= type),
              alpha = 0.2, linetype = 0) +
  scale_color_manual(values=c("Observations"=COL[["simulation"]],
                              "1-box"=col1,
                              "2-box"=col2,
                              "3-box"=col3),
                     breaks = c("Observations", "1-box", "2-box", "3-box"),
                     labels=c("HadCRUT5", "forced response (1-box)", "forced response (2-box)", "forced response (3-box)"), name="") +
  scale_fill_manual(values=c("Observations"=COL[["simulation"]],
                             "1-box"=col1,
                             "2-box"=col2,
                             "3-box"=col3),
                    breaks = c("Observations", "1-box", "2-box", "3-box"),
                    guide = "none") +
  scale_x_continuous(limits=c(1850,2000), expand=c(0.,0.), name = latex2exp::TeX("Time (yr CE)"))+
  scale_y_continuous(name = latex2exp::TeX("Temperature anomaly $(K)$")) +
  theme(legend.position=c(0.25,0.85),
        legend.box.background =  element_rect(colour = "black"),
        legend.title= element_blank(),
        legend.margin = ggplot2::margin(t=-0.02,l=0.05,b=0.05,r=0.1, unit='cm')) +
  theme_td(tsize) + theme(plot.margin = margin(0,0.5,0,0, "cm")) 
print(p1)

p1_ann <- p1 +
  ggplot2::annotate("text",
                    x=1970,
                    y=-0.11,
                    label=paste0("RMSE (1-box): ", round(RMSE_1box, 3), "0 K"), col=col1, size=tsize/3) +
  ggplot2::annotate("text",
                    x=1970,
                    y=-0.18,
                    label=paste0("RMSE (2-box): ", round(RMSE_2box, 3), " K"), col=col2, size=tsize/3) +
  ggplot2::annotate("text",
                    x=1970,
                    y=-0.25,
                    label=paste0("RMSE (3-box): ", round(RMSE_3box, 3), " K"), col=col3, size=tsize/3)
p1_ann

if(save) ggsave("plots/FS6_hadcrut_123box.pdf", device = cairo_pdf, width=6, height=4, dpi=900)


# --------- get parameter table -----------

tbb <- tibble()
fits <- list(fit_1box, fit_2box, fit_3box)

for(i in 1:3) {
  fit = fits[[i]]
  tbb <- rbind(tbb, list_to_tibble(fit$posteriors$parameters) %>%
                 unnest(data) %>%
    add_column(box=paste0(i, "-box")))
}

table_tbb <- list()
for(i in unique(tbb$name)){
  table_tbb[[i]] <- tbb %>% filter(name==i) %>% mutate(
    !!i:=paste0(as.character(round(mean,2)), " (", as.character(round(lower_quant,2)), ",", as.character(round(upper_quant,2)), ")")
  ) %>% select(-median, -mean, -variance, -name, -lower_quant, -upper_quant)
}

table_tbb <- purrr::reduce(table_tbb, left_join)
table_tbb %>% select(box, lambda1, lambda2, lambda3, weights1, weights2, T0, F0, Cap) -> table_tbb

if(save) write.csv2(table_tbb %>% select(-T0, -F0), "tableS3_hadcrut_params.csv", row.names = FALSE, quote=F,  na="")

# ----------- noise -------------------

noise2 <- gen_noise_from_ebm_fit(3, fit_2box)
noise_df2 <- as.data.frame(t(noise2))
noise_df2 <- as.data.frame(apply(noise_df2, 2, function(x) x + fit_2box$posteriors$model_fit$mean))
noise_df2$year <- 1850:2000
pivot_longer(noise_df2, cols = starts_with("V"), names_to = "realisation", values_to = "y") ->
  noise_df2

fit_df = data.frame(year = 1850:2000, y = fit_2box$posteriors$model_fit$mean) %>% add_column(realisation="2-box")
sol_df = data.frame(year = 1850:2000, y = fit_2box$input_params$y_obs) %>% add_column(realisation="Observation")
my_colors <- RColorBrewer::brewer.pal(5, "OrRd")[3:5]
gg_noisy <- ggplot(rbind(noise_df2, sol_df, fit_df)) +
  geom_line(aes(x = year, y = y, col=realisation, size = realisation, alpha=realisation)) +
  labs(x = "Time (yr CE)", y = "Temperature anomaly (K)") +
  scale_x_continuous(expand=c(0,0)) +
  theme_bw() +
  theme_td(tsize) +
  theme(legend.position = "none") +
  scale_color_manual(values=c("2-box"=col2,
                              "Observation"=COL[["simulation"]], 
                              "V1"=my_colors[1],
                              "V2"=my_colors[2],
                              "V3"=my_colors[3]),
                     breaks=c("Observation","2-box","V1","V2","V3"),
                     labels = c("HadCRUT5", "forced", "forced+internal (E1)", "forced+internal (E2)", "forced+internal (E3)")) +
  scale_size_manual(values=c("Observation"=1.4,
                              "2-box"=1.4,
                              "V1"=lsize,
                              "V2"=lsize,
                              "V3"=lsize)) +
  scale_alpha_manual(values=c("Observation"=1,
                              "2-box"=1,
                              "V1"=0.7,
                              "V2"=0.7,
                              "V3"=0.7)) +
  theme(legend.position=c(0.25,0.8),
        legend.box.background =  element_rect(colour = "black"),
        legend.title= element_blank(),
        legend.margin = ggplot2::margin(t=-0.02,l=0.05,b=0.05,r=0.1, unit='cm')) +
  theme_td(tsize) + theme(plot.margin = margin(0,0.5,0,0, "cm")) +
  guides(size="none", alpha="none")
gg_noisy

if(save) 
ggsave("plots/FS5_hadcrut_fit+noise_v2.pdf", device = cairo_pdf, width=6, height=4, dpi=900)
