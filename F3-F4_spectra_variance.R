###THIS SCRIPT IS ORGANIZED AS FOLLOWS###
# 1. CMIP and EMIC spectra deterministic + internal part (Fig. 3)
# 2. Timescale dependent variance ratio of simulated to emulated spectra (Fig. 4)

source("init.R")

saveToPDF <- T

#GET DATA
median=F #consider the mean spectrum
result <- list()
dir.runs <- "output/cmip-emic-fits/"
files_2box <- list.files(dir.runs, pattern="_2box.rda")

for(modeltype in c("CMIP", "EMIC", "EMIC_LO")){
  print(modeltype)

  file <- paste0("output-spectra/spectra_fit_and_noise_", modeltype, ".Rds")

  if(!file.exists(file)){

    #----------------get spectra cmip-------------#
    runs <- c(signal_tbb %>% filter(type==modeltype))$name

    res <- tibble()
    for(i in seq_along(runs)){ #
        print(i)
        if(!is.null(grep(paste0(runs[[i]]), files_2box))){
          res_cmip <-  loadRData(paste0(dir.runs, files_2box)[[grep(paste0(runs[[i]]), files_2box)]])
          res <- rbind(res, noisy_spectra(res_cmip, median=median, debug=T) %>% add_column(run=runs[[i]], box="2box"))
        }
    }
    rm(i, res_cmip)
    res <- res %>% unnest(data) %>% select(run,box,name,process, freq,spec,dof,lim.1, lim.2, "5%", "95%") %>% rename(temperature=spec) %>% group_by(run, box, name, process) %>% nest()
    saveRDS(res, file)
    } else {
    res <- readRDS(file)
    }
  result[[modeltype]] <- res
  rm(res, file)
}

txtsize <- 12.5

#---------------1. CMIP and EMIC spectra deterministic + internal part (Fig. 3) ---------------#
res <- bind_rows(result) %>% filter(process=="spec", box=="2box")
res <- rbind(res %>% filter(!name=="fit+noise"),
             res %>% unnest(data) %>% filter(name=="fit+noise")  %>%
               select(-lim.1, -lim.2) %>% rename(lim.1="5%", lim.2 = "95%")%>% nest()) %>%
  rename(Spec=data)
res <- inner_join(res, signal_tbb %>% filter(!name %in% c("LO2", "LO3", "LO4", "LO5")) %>%
                    rename(run=name)) %>%
                    mutate(type= case_when(type=="EMIC_LO" ~ "EMIC", TRUE ~ type)) %>% ungroup()

#order
runs <- list()
runs[["CMIP"]] <- unique(res %>% filter(type=="CMIP") %>% select(run))
runs[["EMIC"]] <- unique(res %>% filter(type=="EMIC") %>% select(run) %>% arrange(run))

p <- list()

txtsize <- 12.5

res$name <- factor(res$name, c("fit+noise", "fit", "simulation"))

cnt <- 1
for(type in c("EMIC", "CMIP")){
  print(type)
  for(i in runs[[type]]$run){
    print(i)
    p[[type]][[i]] <- plot_spec(res %>% arrange(name)  %>% filter(run==i), c(0.00008,4), c(2,501), name.fill="name",  name.x=latex2exp::TeX('Period $(yr)$')) + theme_td(txtsize) +
      scale_color_manual(values=c("simulation"=COL[["simulation"]], "fit"=COL[["fit"]], "fit+noise"=COL[["fit+noise"]])) +
      scale_fill_manual(values=COL) + annotate("text", 400, 0.0004, label=c(signal_tbb %>% filter(name==i))$alt_name, size=3.0, hjust = 0) +
      annotate("text", 400, 1.8, label=letters[[cnt]], size=3.0, fontface =2, hjust = 0) +
      scale_y_log10(breaks=c(0.001, 0.1), name=NULL, labels = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0,0), limits=c(0.00008,4)) +
      theme(legend.position="bottom", legend.direction='vertical', legend.title = element_blank(), legend.text = element_text(size=14))
    print(p[[type]][[i]])
    cnt <- cnt + 1
  }
}

pall <- lapply(p, function(x) lapply(x, function(y) y +
                 theme(
                   legend.position = "none",
                   strip.text.x = element_blank(),
                   plot.margin = unit(c(0,0,0,0), "cm")
                 )
  )
)

leg <- get_only_legend(p[[1]][[1]] +
                         scale_color_manual(values=c( "fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]],  "simulation"=COL[["simulation"]]), labels=c("forced + internal variability", "forced variability","target simulation")) +
                         theme( legend.box.background = element_rect(colour = "black"))
                       )

ncol=4
for(type in c("CMIP", "EMIC")){
  N = length(pall[[type]])

  for(i in 1:N){
    pall[[type]][[i]] <- pall[[type]][[i]] +
      annotation_logticks(sides = "rl",
                          short = unit(0,"mm"),
                          mid = unit(0,"mm"),
                          long = unit(0.5,"mm"))
    if(!i %in% seq(1, N, ncol)){
      pall[[type]][[i]] <- pall[[type]][[i]]  + theme(
        axis.ticks.y=element_line(colour="black"),
        axis.title.y=element_blank(),
        axis.text.y = element_blank())
    }
    if(!i %in% c(N-3, N-2, N-1, N)){
      pall[[type]][[i]] <- pall[[type]][[i]]  +
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank())
    }
  }
}

type="CMIP"
pend <- pall[[type]][[length(pall[[type]])]]
p1 <- cowplot::plot_grid(
  plot_grid(plotlist=pall[[type]][c(1,2,3,4)], ncol=ncol, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
  plot_grid(plotlist=pall[[type]][c(5,6,7,8)], ncol=ncol, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
  NULL,
  plot_grid(pall[[type]][[9]], pend, NULL, NULL, ncol=ncol, axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,0.98,1.02,1)),
  ncol=1, align="v", rel_heights = c(1,1.3,-0.47,1.3))

type="EMIC"
pend <- pall[[type]][[length(pall[[type]])]]
p2 <- cowplot::plot_grid(
  plot_grid(plotlist=pall[[type]][c(1,2,3,4)], ncol=ncol, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
  plot_grid(plotlist=pall[[type]][c(5,6,7,8)], ncol=ncol, align="h", axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,1,1,1)),
  NULL,
  plot_grid(pall[[type]][[9]], pall[[type]][[10]], pend, leg, ncol=ncol, axis="tblr", scale=1.0, greedy=F, rel_widths = c(1.2,.99,.98,1.01)),
  ncol=1, align="v", rel_heights = c(1,1.3,-0.47,1.3))

p1 <- ggpubr::annotate_figure(p1, top=ggpubr::text_grob("CMIP5 models", x=0.05, face="bold"),
                              left=ggpubr::text_grob(latex2exp::TeX('Power spectral density (PSD) $(K^2 yr)$'), rot=90, size=txtsize))
print(p1)
p2 <- ggpubr::annotate_figure(p2, top=ggpubr::text_grob("AR5 EMICs", x=0.05, face="bold"),
                              left=ggpubr::text_grob(latex2exp::TeX('Power spectral density (PSD) $(K^2 yr)$'), rot=90, size=txtsize))
print(p2)

cowplot::plot_grid(p2,p1, nrow=2)
if(saveToPDF) ggsave(paste0("plots/F3_CMIP_EMIC_spectra.pdf"), width=12, height=8, dpi=900, device = cairo_pdf)


#-------------2. Timescale dependent variance ratio of simulated to emulated spectra (Fig. 4)---------------#
res <- bind_rows(result) %>% filter(process=="spec", box=="2box")
res <- rbind(res %>% filter(!name=="fit+noise") %>% unnest(data) %>% rename(spec="temperature") %>% nest(),
             res %>% unnest(data) %>% filter(name=="fit+noise")  %>%
               select(-lim.1, -lim.2) %>% rename(lim.1="5%", lim.2 = "95%", spec="temperature") %>% nest()) %>%
               rename(Spec=data)
res <- inner_join(res, signal_tbb %>% rename(run=name))


var_tibble <- tscale_var(res)
var_tibble$var <- as.numeric(var_tibble$var)

get_ratios <- function(tibble, tar="simulation"){
  idx <- which(tibble$name == tar)
  vartar <- tibble$var[[idx]]
  tibble <- tibble %>% mutate(varratio=var/vartar)
  return(tibble)
}

var_split <- var_tibble %>% group_by(run, box, tscale, process) %>% group_split()
var_joint <- lapply(var_split, function(x){get_ratios(x, tar="simulation")})
var_tibble <- do.call(rbind, var_joint)

var_tibble$name <- factor(var_tibble$name, c("simulation", "fit+noise",  "fit"))
var_tibble$tscale <- factor(var_tibble$tscale, rev(c("interannual", "decadal", "multidecadal", "centennial")))
var_tibble$type <- factor(var_tibble$type, c("EMIC", "EMIC_LO", "CMIP"))

facet.labs <- c("a AR5 EMICs", "b LOVECLIM ensemble members", "c CMIP5 models")
names(facet.labs) <-  c("EMIC", "EMIC_LO", "CMIP")

var_tibble <- var_tibble %>% mutate(type=replace(type, run=="LO0", "EMIC"))

ggplot(var_tibble %>% filter(name!="simulation"), aes(x=tscale, y=varratio, color=name, fill=name)) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.3, dodge.width=0.), size = 4, alpha = .5) +
  stat_summary(
    fun = base::mean, geom = "point",
    shape = 95, size = 30, stroke=0.6, alpha=1
  ) +
  ylab("Varriance ratio (emulated / target)") +
  theme_bw() + theme_td(txtsize) +
  scale_y_log10(limits=c(0.09, 3), breaks=c(0.1, 0.3,1,3), labels=c( 0.1, 0.3, 1, 3)) +
  scale_fill_manual(values=c( "fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]]), labels=c("forced + internal variability",  "forced variability")) +
  scale_color_manual(values=c("fit+noise"=COL[["fit+noise"]], "fit"=COL[["fit"]]), labels=c("forced + internal variability",  "forced variability")) +
  facet_wrap(~type, scale="fix",  labeller=as_labeller(facet.labs)) +
  theme(legend.position="none",
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(.5, "cm"),
        axis.title.x=element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(size=txtsize-1),
        strip.text.x = element_text(face = "bold", hjust = 0, vjust=-1, size=txtsize))

if(saveToPDF) ggsave(paste0("plots/F4_varriance_ratios.pdf"), width=12, height=4., dpi=900, device = cairo_pdf)

