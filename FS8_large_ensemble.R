source("init.R")
load("output/lme-fits/LME_ensemble_mean_2000.rda")
load("output/lme-fits/cesm_le_2000_member_1.rda")

saveToPDF <- F

lsize <- 0.7 #line size
tsize <- 12.5 #textsize

temp_vals = fit_2box$input_params$y_obs
start_year = fit_2box$meta$start_year

offset = ens_mean$ens_mean[1] - temp_vals[1]
ens_mean_norm = ens_mean$ens_mean - ens_mean$ens_mean[1] + offset

data = rbind(
        tibble(year =start_year:(start_year+length(temp_vals)-1),
               temp = temp_vals - temp_vals[1],
               type = "CESM LME (E1)",
                temp_lower = NA, 
                temp_upper = NA),
        tibble(year = ens_mean$year,
                temp = ens_mean_norm,
                type = "CESM LME (mean)",
                temp_lower = NA, 
                temp_upper = NA),
        tibble(year = start_year:(start_year+length(temp_vals)-1),
                                 temp = fit_2box$posteriors$model_fit$mean,
                                 temp_lower = fit_2box$posteriors$model_fit$lower_quant,
                                 temp_upper = fit_2box$posteriors$model_fit$upper_quant,
                                 type = "forced response")
)

col_forc = COL[["forcing"]]
COL[["mean"]] = "deepskyblue2"

scl=1/20
y_shift= -1

ggplot() + 
    geom_line(aes(x=start_year:(start_year+length(temp_vals)-1),y=fit_2box$input_params$forc_vals*scl +y_shift), size=lsize, colour=col_forc) +
    geom_line(data=data, aes(x = year, y = temp, col=type), lwd=lsize) +
    scale_color_manual(values= c("forced response"=COL[["fit"]], 
                                 "CESM LME (mean)"=COL[["mean"]],
                                 "CESM LME (E1)"=COL[["simulation"]])) +
    scale_x_continuous(limits=c(850,2000), expand=c(0,0), name = latex2exp::TeX("Time (yr CE)")) +
    geom_ribbon(data = data %>% filter(type == "forced response"),
                  ggplot2::aes(x = year, ymax = temp_upper, ymin = temp_lower,
                               fill = type),
                  linetype = 0, alpha=.3) +
      ggplot2::scale_fill_manual(values= c("forced response" = COL[["fit"]])) +
      ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(fill=NA)), fill="none") +
    scale_y_continuous(name = latex2exp::TeX("Temperature anomaly $(K)$"), 
    sec.axis = sec_axis(~ (. - y_shift)/scl, 
                breaks=c(-20,-10,0),
                name = latex2exp::TeX("Radiative forcing anomaly $(W/m^2)$"))) +
    theme_bw() +
    theme(legend.position=c(0.77,0.12),
        legend.box.background =  element_rect(colour = "black"),
        legend.title= element_blank(),
        axis.text.y.right=element_text(margin = unit(c(t = 0, r = 0, b = 0, l = 2.), "mm"), colour=col_forc),
        axis.ticks.y.right=element_line(colour=col_forc),
        axis.title.y.right=element_text(colour=col_forc, hjust=1)) + 
        theme_td(tsize) +
        theme(plot.margin = margin(0.1,0.1,0,0, "cm")) 

if(saveToPDF) ggsave("plots/FS8_lme_comparison.pdf", device = cairo_pdf, width=6, height=4.5, dpi=900)

