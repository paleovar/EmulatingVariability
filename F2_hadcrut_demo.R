source("init.R") #load required packages, fix parameters and functions

saveToPdf <- F

#load data
fit_1box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_1box.rda")
fit_2box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_2box.rda")
fit_3box = loadRData("output/hadcrut-fits-cea-sbf/hadcrut_3box.rda")

#useful functions
create_labs_dict <- function(X_names) {
  d = length(X_names)
  labs_dict = list("T0" = unname(latex2exp::TeX("$T_0 \\, (K)$")),
                   "F0" = unname(latex2exp::TeX("$F_0 \\, (W m^{-2})$")),
                   "Cap" = unname(latex2exp::TeX("$C \\, (W yr m^{-2} K^{-1})$")))
  for(j in 1:d) {
    var_name = X_names[j]
    if(startsWith(var_name, "lambda")) {
      ind = sub("lambda", "", var_name)
      labs_dict[[var_name]] = unname(latex2exp::TeX(paste0("$\\lambda_", ind, " \\, (yr^{-1})$")))
    } else if(startsWith(var_name, "weights")) {
      ind = sub("weights", "", var_name)
      labs_dict[[var_name]] = unname(latex2exp::TeX(paste0("$\\w_", ind, "$")))
    }
  }
  return(labs_dict)
}

#print parameters
print("Estimated parameters for beta prior")

for (par in names(fit_2box$posteriors$parameters)) {
  posteriors = fit_2box$posteriors$parameters[[par]]
  print(paste0(par, ": ", round(posteriors$mean,2),
               " (", round(posteriors$lower_quant, 2), ", ",
               round(posteriors$upper_quant, 2), ")"))
}

#define plot parameters
lsize <- 0.7 #line size
tsize <- 12.5 #textsize

plots2 = plot_fit(fit_2box)
gg_marginals2  = plots2$gg_marginals

test_2box <- noisy_spectra(fit_2box, df.log=0.08, median=F, debug=F, n_samples = 1000)

#CIs
test_2box <- rbind(test_2box %>% filter(!name=="fit+noise") %>% group_by(name, process),
                   test_2box %>% unnest(data) %>% filter(name=="fit+noise")  %>%
               select(-lim.1, -lim.2) %>% rename(lim.1="5%", lim.2 = "95%")%>% group_by(name, process) %>% nest())


test_2box <- test_2box %>% filter(process=="spec")

#upper left pannel
scl=1/4
y_shift= -0.5
col_forc = COL[["forcing"]]
p1 <- ClimBayes::plot_fit(fit_2box, line_size=lsize, textsize=tsize) + scale_fill_manual(values=c("Observations"=COL[["simulation"]], "Model est."=COL[["fit"]]), guide=F) +
  geom_line(aes(y=rep(fit_2box$input_params$forc_vals*scl +y_shift,2)), size=lsize, colour=col_forc) +
  scale_x_continuous(limits=c(1850,2000), expand=c(0,0), name = latex2exp::TeX("Time (yr CE)"))+
  scale_y_continuous(name = latex2exp::TeX("Temperature anomaly $(K)$"), sec.axis = sec_axis(~ (. - y_shift)/scl, name = latex2exp::TeX("Radiative forcing anomaly $(W/m^2)$"))) +
  scale_color_manual(values=c("Observations"=COL[["simulation"]], "Model est."=COL[["fit"]]), labels=c("HadCRUT5", "forced response (2-box)"), name="") +
  theme(legend.position=c(0.3,0.88),
        legend.box.background =  element_rect(colour = "black"),
        legend.title= element_blank(),
        axis.text.y.right=element_text(margin = unit(c(t = 0, r = 0, b = 0, l = 2.), "mm"), colour=col_forc),
        axis.ticks.y.right=element_line(colour=col_forc),
        axis.title.y.right=element_text(colour=col_forc)) + theme_td(tsize)
print(p1)

#upper right pannel
gg <- ClimBayes::plot_marginals(fit_2box, return_all=T, line_size=lsize, textsize=tsize, labs_dict=create_labs_dict(fit_2box$input_params$X_names))
gg$gg <- NULL
bounds <- list(
  c(fit_2box$input_params$lambda1_lb, fit_2box$input_params$lambda1_ub),
  c(fit_2box$input_params$lambda2_lb, fit_2box$input_params$lambda2_ub),
  c(fit_2box$input_params$weights1_lb,fit_2box$input_params$weights1_ub),
  c(fit_2box$input_params$T0_lb,fit_2box$input_params$T0_ub),
  c(fit_2box$input_params$F0_lb,fit_2box$input_params$F0_ub)
)

plots <- list()
for(i in seq_along(gg$plot_list)){
  plots[[i]] <- gg$plot_list[[i]] + theme_td(tsize) +
      scale_color_manual(values=c("Posterior"=COL[["Posterior"]], "Prior"=COL[["Prior"]])) +
      scale_x_continuous(limits=bounds[[i]], breaks=c(bounds[[i]][1], sum(bounds[[i]])/2, bounds[[i]][2])) +
      theme(axis.title.y=element_blank(),
            legend.position="none")
}
leg <- get_legend(gg$plot_list[[i]] + theme_td(tsize) +
                    scale_color_manual(values=c("Prior"=COL[["Prior"]], "Posterior"=COL[["Posterior"]]), labels=c("prior", "posterior")) +
                    scale_x_continuous(limits=bounds[[i]],  breaks=c(bounds[[i]][1], sum(bounds[[i]])/2, bounds[[i]][2])) +
                    theme(legend.box.background = element_rect(colour = "black"),
                          legend.direction = "vertical")
)
plots[[i+1]] <- ggpubr::as_ggplot(leg)

p2 <- cowplot::plot_grid(plotlist=plots, align="hv")
print(p2)

#lower left pannel
test_2box$name <- factor(test_2box$name, c("fit+noise", "fit", "simulation"))
p3 <- plot_spec(test_2box %>% arrange(desc(name)), c(0.00004,0.7), c(2.1,153), name.y=latex2exp::TeX('Power spectral density (PSD) $(K^2 yr)$'), name.x=latex2exp::TeX('Period $(yr)$'), name.fill="name", name.data="data") +
  theme_td(tsize) +
  scale_fill_manual(values=c( "fit+noise"=COL[["fit+noise"]],"fit"=COL[["fit"]],  "simulation"=COL[["simulation"]]), labels=c( "forced + internal variability", "forced variability", "HadCRUT5")) +
  scale_color_manual(values=c( "fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]],  "simulation"=COL[["simulation"]]), labels=c("forced + internal variability", "forced variability",  "HadCRUT5")) +
  theme(legend.position = c(0.32, 0.18),
        legend.box.background = element_rect(colour = "black")
        )
print(p3)

#lower right pannel
var_tibble <- tscale_var(test_2box %>% rename(Spec=data), tscales=list("interannual"=list("t1"=2, "t2"=5),
                                                 "decadal"=list("t1"=5,"t2"=20),
                                                 "multidecadal"=list("t1"=20,"t2"=100)))#,
                                                 #"centennial"=list("t1"=100, "t2"=150)))
var_tibble$var <- as.numeric(var_tibble$var)
var_tibble$dof <- as.numeric(var_tibble$dof)

get_ratios <- function(tibble, tar="simulation"){
  idx <- which(tibble$name == tar)
  vartar <- tibble$var[[idx]]
  doftar <- tibble$dof[[idx]]
  tibble <- tibble %>% mutate(varratio=var/vartar)
  tibble <- tibble %>% rowwise %>% mutate(q.low=ConfRatio(varratio=varratio, df.1=doftar, df.2=dof)[1], q.up=ConfRatio(varratio=varratio, df.1=doftar, df.2=dof)[2])

  return(tibble)
}

var_split <- var_tibble %>% group_by(process, tscale) %>% group_split()
var_joint <- lapply(var_split, function(x){get_ratios(x, tar="simulation")})
var_tibble <- do.call(rbind, var_joint)  %>% ungroup()

var_tibble$name <- factor(var_tibble$name, c("simulation", "fit+noise",  "fit"))
var_tibble$tscale <- factor(var_tibble$tscale, rev(c("interannual", "decadal", "multidecadal")))
#var_tibble$tscale <- factor(var_tibble$tscale, rev(c("interannual", "decadal", "multidecadal", "centennial")))

p4 <- ggplot(var_tibble %>% filter(name!="simulation"), aes(x=tscale, y=varratio, color=name, fill=name)) +
  geom_errorbar(aes(ymin=q.low, ymax=q.up), width=0., alpha=0.8) +
  geom_point(size = 4, alpha = .5) +
  ylab("Varriance ratio (emulated / target)") +
  theme_bw() + theme_td(tsize) +
  scale_y_log10(limits=c(0.03, 1.3), breaks=c(0.03, 0.1, 0.3, 1), labels=c(0.03, 0.1, 0.3, 1)) +
  xlab("Timescale")+
  geom_hline(yintercept = 1) +
  scale_fill_manual(values=c( "fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]]), labels=c( "forced + internal variability",  "forced variability")) +
  scale_color_manual(values=c("fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]]), labels=c("forced + internal variability",  "forced variability")) +
  theme(legend.position=c(0.29, 0.15),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(.1, "cm"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", hjust = 0, vjust=-1),
        legend.box.background = element_rect(colour = "black")
        )
print(p4)

#plot all together
y.grob <- grid::textGrob("Density", gp=grid::gpar(col="black", fontsize=tsize), rot=90)
p2 <- gridExtra::grid.arrange(gridExtra::arrangeGrob(p2, left = y.grob))

t1 <- cowplot::plot_grid(p1,p3, nrow=2, ncol=1, align="v", labels=c("a", "c"))
t2 <- cowplot::plot_grid(p2,p4, nrow=2, ncol=1, align="v", labels=c("b", "d"))

cowplot::plot_grid(t1, t2, ncol=2, align="hv")

#save
if(saveToPdf) ggsave("plots/F2_hadcrut_framework.pdf", width=12, height=8, dpi=900, device = cairo_pdf)
