source("init.R")
#load data
dir.runs <- "output/cmip-emic-fits/"  #source data directory
files_2box = list.files(dir.runs, pattern="_2box.rda")  #files in source data directory

txtsize <- 12.5

saveToPDF <- F #TRUE for saving the plot as .pdf  in the `./plots/` folder

signal_tbb <- signal_tbb %>% mutate(type=replace(type, name=="LO0", "EMIC"))

cnt=0
for(modeltype in c("EMIC", "EMIC_LO", "CMIP")){
  runs <- sort(c(signal_tbb %>% filter(type==modeltype))$name)
  p <- list()
  for (i in seq_along(runs)){
    cnt <- cnt+1
    print(i)
    if(!is.null(grep(paste0(runs[[i]]), files_2box))){
      res_2box <-  loadRData(paste0(dir.runs, files_2box)[[grep(paste0(runs[[i]]), files_2box)]])
    }

    p[[i]] <- plot_fit(res_2box) +
      scale_x_continuous(limits= c(850,1850), name="Time (yr CE)", expand=c(0.01,0.01)) +
      scale_y_continuous(limits=c(-1.9,0.6), breaks=c(0.5, 0,-.5, -1,-1.5, -2), labels=c("", "0", "", "-1", "", "")) +
      theme_td(txtsize) +
      annotate("text", 900, 0.5, label=letters[[cnt]], size=3.0, fontface =2, hjust = 0) +
      theme(legend.position=c(0.78,0.15), plot.margin = margin(l=0.1, b=0.0, unit="cm")) + #plot.margin
      scale_color_manual(
        values = c("Observations" = COL[["simulation"]], "Model est." = COL[["fit"]]),
        breaks=c("Observations", "Model est."),
        labels=c("target simulation", "forced response (2-box)")) +
      scale_fill_manual(
        values = c("Model est." = COL[["fit"]])) +
      guides(fill="none")

      p[[i]] <- p[[i]] + annotate("text", 900, -1.6, label=c(signal_tbb %>% filter(name==runs[[i]]))$alt_name, size=3.0, hjust = 0)

      if(modeltype!="EMIC_LO"){
        if(!i %in% c(7,8,9,10)){
          p[[i]] <- p[[i]]  + theme(
            axis.title.x=element_blank(),
            axis.text.x = element_blank())
        }
        if(!i %in% c(1,5,9)){
          p[[i]] <- p[[i]]  + theme(
            axis.text.y = element_blank())
        }
      }

      if(modeltype=="EMIC_LO"){
        if(i == 1){
          p[[i]] <- p[[i]]  + theme(
            axis.title.x=element_blank(),
            axis.text.x = element_blank())
        }
        if(!i %in% c(1,5)){
          p[[i]] <- p[[i]]  + theme(
            axis.text.y = element_blank())
        }
      }
  }
  pall <- lapply(p, function(x) x +
                 theme(
                   legend.position = "none",
                   strip.text.x = element_blank(),
                   axis.title.y=element_blank()
                 )
              )

  #plot AR5 EMICs
  if(modeltype=="EMIC"){
    leg <- get_only_legend(p[[1]] + theme(legend.position=c(1,0.4),
                                          legend.key.size = unit(1.25, 'cm'),
                                          legend.box.background = element_rect(colour = "black"),
                                          plot.margin = unit(c(1, 1, 1, 1), "cm")))

    p1 <- cowplot::plot_grid(
      plot_grid(plotlist=pall[c(1,2,3,4)], ncol=4, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      plot_grid(plotlist=pall[c(5,6,7,8)], ncol=4, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      NULL,
      plot_grid(pall[[9]], pall[[10]], leg, NULL, ncol=4, axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      ncol=1, align="v", rel_heights = c(1,1.3,-0.35,1.3)
    )

    p1 <- ggpubr::annotate_figure(p1, top=ggpubr::text_grob("AR5 EMICs", x=0.05, face="bold"),
                                  left=ggpubr::text_grob(latex2exp::TeX('Temperature anomaly (K)'), rot=90, size=txtsize))
  }

  #plot LOVECLIM ensemble members
  if(modeltype=="EMIC_LO"){
    p2 <- cowplot::plot_grid(
      plot_grid(plotlist=pall[c(1,2,3,4)], ncol=4, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      NULL,
      plot_grid(pall[[5]], NULL, NULL, NULL, ncol=4, axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      ncol=1, align="v", rel_heights = c(1.3,-0.35,1.3)
    )

    p2 <- ggpubr::annotate_figure(p2, top=ggpubr::text_grob("LOVECLIM ensemble members", x=0.12, face="bold"),
                                  left=ggpubr::text_grob(latex2exp::TeX('Temperature anomaly (K)'), rot=90, size=txtsize))
  }

  #plot CMIP5 simulations
  if(modeltype=="CMIP"){
    p3 <- cowplot::plot_grid(
      plot_grid(plotlist=pall[c(1,2,3,4)], ncol=4, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      plot_grid(plotlist=pall[c(5,6,7,8)], ncol=4, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      NULL,
      plot_grid(pall[[9]], pall[[10]], NULL, NULL, ncol=4, axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
      ncol=1, align="v", rel_heights = c(1,1.3,-0.35,1.3)
        )

    p3 <- ggpubr::annotate_figure(p3, top=ggpubr::text_grob("CMIP5 models", x=0.06, face="bold"),
                                left=ggpubr::text_grob(latex2exp::TeX('Temperature anomaly (K)'), rot=90, size=txtsize))
    }

  if(modeltype=="CMIP") cowplot::plot_grid(p1, p2, p3, nrow=3, rel_heights = c(1,0.7,1))
  if(saveToPDF) ggsave(paste0("plots/FS7_fits_supp.pdf"), width=12, height=12, dpi=800, device=cairo_pdf)
}

