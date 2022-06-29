#required packages
library(tibble)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(zoo)
library(PaleoSpec)
library(ClimBayes)
library(purrr)
library(cowplot)
library(stringr)
library(latex2exp)
library(ggpubr)

#dictionary
cmip_tbb <- tibble(
  name=c("bcc_csm1", "csiro", "giss", "hadgem", "miroc", "ccsm4", "fgoals", "hadcm3", "ipsl", "mpi_esm"),
  alt_name=c("BCC-CSM1-1", "CSIRO-Mk2L-1-2", "GISS-E2-R", "HadGEM2-ES", "MIROC-ESM", "CCSM4", "FGOALS-s2", "HadCM3", "IPSL-CM5A-LR", "MPI-ESM-P")
) %>% add_column(type="CMIP")

emic_tbb1 <- tibble(
  name=c("B3", "C2", "C3", "DC", "I2", "ME", "MI", "UM", "UV"),
  alt_name=c("Bern 3D", "CLIMBER2 (mean)", "CLIMBER-3alpha", "DCESS ESM v1", "IGSM 2.2",
             "MESMO 1.0", "MIROC-lite", "UMD", "UVic v2.9")
) %>% add_column(type="EMIC")

emic_tbb2 <- tibble(
  name=c("LO0", "LO1", "LO2","LO3", "LO4", "LO5"),
  alt_name=c("LOVECLIM V.1.2 (mean)","LOVECLIM V.1.2 (E1)",  "LOVECLIM V.1.2 (E2)",  "LOVECLIM V.1.2 (E3)",  "LOVECLIM V.1.2 (E4)",  "LOVECLIM V.1.2 (E5)")
) %>% add_column(type="EMIC_LO")

signal_tbb <- full_join(cmip_tbb, rbind(emic_tbb1, emic_tbb2))

#define theme
#colors
COL <- c(
  "simulation" = "grey50",
  "fit" = "#0F2080",
  "fit+noise" = "#F5793A",
  "Prior" = "grey50",
  "Posterior" = "#0F2080",
  "forcing" = "#2AAF74"
)

#plotting style
theme_td <- function(textsize=12){
theme(axis.title = element_text(size = textsize),
      axis.text = element_text(size = textsize),
      axis.text.x = element_text(margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm"), color="black"),
      axis.text.y = element_text(margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"), color="black"),
      axis.ticks.length = unit(-1.4, "mm")) +
  theme(legend.key = element_rect(color = "transparent"),
        legend.background=element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=textsize),
        legend.key.height = unit(0.2, "cm")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.minor = element_blank())
}

#--------------------------------------------------------------#
#' @title PSD plot
#' @description Plots the power spectral density in log-log space.
#' taken from Ellerhoff & Rehfeld, 2021 (https://zenodo.org/record/5923109, DOI: 10.1103/PhysRevE.104.064136)
#' @param tibble nested tibble (name, data) with column Spec(freq, spec, dof, lim.1, lim.2) is the nested spectrum
#' @param main title
#' @param name.y y-axis label
#' @param name.x x-axis label
#' @param name.data string, name of the nested spectrum column, i.e. either "Spec" or "Smooth_spec"
#' @param name.col string, grouping feature for colored lines
#' @param name.fill string, grouping feature for filled confidence bands
#' @param name.line string, grouping feature for linetype
#' @param name.alpha string, grouping feature for transparence of filled ribbons
#' @param conf.band TRUE / FALSE denotes whether confidence bands of spectra are plotted or not
#' @return ggplot object
#' @export
plot_spec <- function(tibble, ylims, xlims, name.y=latex2exp::TeX('PSD $S(\\tau)\\, (K^2 yr)$ '), name.x=latex2exp::TeX('period $\\tau\\,(yr)$'), name.data = "Spec", name.col="name", name.fill=NULL, name.line=NULL, name.alpha=NULL, conf.bands=T){

  #--------------------------------------------------------------#
  #' @title Reverse log transform
  #' @description reverse log transformation for plotting in logarithmic space. Useful when plotting power spectral densities.
  #' adapted from the metR package https://rdrr.io/cran/metR/src/R/reverselog_trans.R by Elio Campitelli
  #' @param base Base of the logarithm, default value 10
  #' @return axis object
  #' @export
  reverselog_trans <- function(base = 10) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
  }
  #--------------------------------------------------------------#

  if(is.null(name.fill)){name.fill<-name.col}
  xlims <- rev(sort(xlims))
  ylims <- sort(ylims)
  yrs.period <- rev(c(0.0001, 0.001, 0.01,  0.1, 1, 10, 100, 1000, 10000, 100000, 1000000))
  yrs.labels <- rev(c(latex2exp::TeX('$10^{-4}$'),latex2exp::TeX('$10^{-3}$'), latex2exp::TeX('$10^{-2}$'), latex2exp::TeX('$10^{-1}$'), latex2exp::TeX('$10^{0}$'), latex2exp::TeX('$10^{1}$'), latex2exp::TeX('$10^{2}$'), latex2exp::TeX('$10^{3}$'), latex2exp::TeX('$10^{4}$'), latex2exp::TeX('$10^{5}$'), latex2exp::TeX('$10^{6}$')))
  if("temperature" %in% names(tibble %>% tidyr::unnest(name.data))){
    tmp <- tibble %>% tidyr::unnest(name.data) %>%
      ggplot2::ggplot(ggplot2::aes(x = 1/freq, y=temperature))
  } else if ("spec" %in% names(tibble %>% tidyr::unnest(name.data))){
    tmp <- tibble %>% tidyr::unnest(name.data) %>% rename(temperature=spec) %>%
      ggplot2::ggplot(ggplot2::aes(x = 1/freq, y=temperature))
  } else {
    print("no column with spectral data found")
  }
  if(conf.bands==T){
    if(!is.null(name.alpha)){
      tmp <- tmp + ggplot2::geom_ribbon(ggplot2::aes_string(ymin="lim.1", ymax="lim.2", fill=name.fill, alpha=name.alpha), linetype = 0) +
        ggplot2::guides(fill=FALSE)
    } else {
      tmp <- tmp + ggplot2::geom_ribbon(alpha=0.2, ggplot2::aes_string(ymin="lim.1", ymax="lim.2", fill=name.fill), linetype = 0) +
        ggplot2::guides(fill=FALSE)
    }
  }
  tmp + ggplot2::geom_line(size=0.7, ggplot2::aes_string(color=name.col, linetype=name.line)) + ggplot2::theme_bw() +
    ggplot2::scale_y_log10(name=name.y, labels = scales::trans_format("log10", scales::math_format(10^.x)), expand=c(0.0, 0.0), limits=ylims)  +
    ggplot2::scale_x_continuous(trans=reverselog_trans(10), breaks = yrs.period, labels = yrs.labels, name=name.x,expand=c(0.0, 0.0), limits=xlims)
}

#--------------------------------------------------------------#
#' @title Equidistant timeseries from tibble
#' @description  Wrapper for applying PaleoSpec::MakeEquidistant() to a tibble with data column that contains the timeseries
#' taken from Ellerhoff & Rehfeld, 2021 (https://zenodo.org/record/5923109, DOI: 10.1103/PhysRevE.104.064136)
#' @param y tibble object with nested data column
#' @return tibble object with two nested tibble. The "data" column the raw data and the "EquiTS" contains the equidistant timeseries.
#' @export
equidistant <- function(y){
  cnt <- 0
  if(!"interp.res" %in% names(y)){
    y <- y %>% tibble::add_column(interp.res=1)
    print("timestep set to 1")
  }
  y_equi <- y %>% tibble::add_column(EquiTS =
                               .$data %>% lapply(., function(x){
                                 cnt <<- cnt+1
                                 tmp <- x %>% tibble::as_tibble()
                                 tmpdf <- list()
                                 for (i in names(tmp)){
                                   if(i=="year"){next}
                                    tmpdf[[i]] <- PaleoSpec::MakeEquidistant(tmp$year, tmp[[i]], dt = y$interp.res[cnt])
                                 }
                                 tmpdf %>% tibble::as_tibble()
                               }
                               )
  )
  return(y_equi)
}

#--------------------------------------------------------------#
#' @title Computing the spectrum over a tibble
#' @description wrapper for cmputing the spectrum from a tibble using the PaleoSpec pacakge
#' taken from Ellerhoff & Rehfeld, 2021 (https://zenodo.org/record/5923109, DOI: 10.1103/PhysRevE.104.064136)
#' @param y tibble with data EquiTS column that contains interpolated time series
#' @param k a positive integer, the number of tapers, often 2*nw.
#' @param nw a positive double precision number, the time-bandwidth parameter.
#' @return tibble with nested spectrum
#' @export
tibble_spec <- function(y, k=3, nw=2){
  cnt <- 0
  y_spec <- y %>% tibble::add_column(Spec =
                                       .$EquiTS %>% lapply(., function(x){
                                         cnt <<- cnt+1
                                         tmp <- x %>% tibble::as_tibble()
                                         tmpdf <- list()
                                         for (i in names(tmp)){
                                           target <- tmp[[i]]
                                           if(any(is.na(tmp[[i]]))){
                                             if("model" %in% names(y)){
                                               print(paste0(y$model[[cnt]], " ", i, " contains NA values, which were removed"))
                                               target = tseries::na.remove(tmp[[i]])
                                             }
                                             if("Name" %in% names(y)){
                                               print(paste0(y$Name[[cnt]], " ", i, " contains NA values, which were removed"))
                                               target = tseries::na.remove(tmp[[i]])
                                             }
                                           }
                                           target <- target - mean(target)
                                           restmp <- PaleoSpec::SpecMTM(target, k, nw, detrend=TRUE)
                                           restmp <- PaleoSpec::LogSmooth(restmp, df.log=0.04)
                                           tmpdf[[i]] <- restmp$spec
                                         }
                                         tmpdf[["freq"]] <- restmp$freq
                                         tmpdf[["dof"]] <- restmp$dof
                                         tmpdf[["lim.1"]] <- restmp$lim.1
                                         tmpdf[["lim.2"]] <- restmp$lim.2
                                         tmpdf %>% tibble::as_tibble()
                                       }
                                       )
  )
  return(y_spec)
}

#----------------------------#
#' @title Mean Spectrum
#' @description Wrapper for the PaleoSpec::MeanSpectrum function. We introduced a correction to the number of records,
#' in case spectra where only partially overlapping.
#' taken from Ellerhoff et al. (2022)
#' @param speclist list of spectra
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @param weights vector of weights (same length as elements in speclist)
#' @return list(spec,nRecords) spec=average spectrum, nRecords = number of records contributing to each spectral estimate
#' @export
MeanSpec <- function(specList, iRemoveLowest = 1, weights = rep(1, length(specList))){
  meanspec <- PaleoSpec::MeanSpectrum(specList, iRemoveLowest, weights)
  meanspec$spec$spec <- meanspec$spec$spec * length(specList)/meanspec$nRecord
  meanspec$spec <- AddConfInterval(meanspec$spec)
  return(meanspec)
}

#--------------------------------------------------------------#
#' @title Speclist to tibble spec converter
#' @description Converts a list with objects of class "spec" to tibble with nested tibbles that contain the spectrum.
#' Inverts the tibble_to_list() functions.
#' taken from Ellerhoff and Rehfeld (2021) (https://zenodo.org/record/5923109, DOI: 10.1103/PhysRevE.104.064136)
#' @param list_spec_obj  list of objects of class "spec"
#' @return tibble object that contains the spectrum
#' @export
list_to_tibble <- function(list_spec_obj){

  df <- tibble(name=character(), data=list())

  for (i in 1:length(list_spec_obj)){

    if("mtm" %in% class(list_spec_obj[[i]])){
      class(list_spec_obj[[i]]) <- "spec"
      list_spec_obj[[i]]$mtm <- NULL
    }

    if(class(list_spec_obj[[i]])=="zoo"){
      time <- index(list_spec_obj[[i]])
      temp <- coredata(list_spec_obj[[i]])
      list_spec_obj[[i]] <- list(time=time, temp=temp)
    }

    class(list_spec_obj[[i]]) <- "list"

    df <- df %>% add_row(
      name = names(list_spec_obj)[[i]],
      data = list_spec_obj[[i]] %>% as.data.frame() %>% as_tibble() %>% list()
    )
  }
  return(df)
}

#--------------------------------------------------------------#
#' @title Spectral Gain
#' @description Computes the spectral gain by dividing the output by the input spectrum, based on functions from the PaleoSpec package. Taken from
#' Ellerhoff & Rehfeld, 2021, https://zenodo.org/record/5923109, DOI: 10.1103/PhysRevE.104.064136
#' @param specList list that contains at least two objects of class "spec"
#' @param input index of the input spectrum
#' @param output index of the output spectrum
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @return object of class "spec"
#' @export
transferSpec <- function (specList, input=input.spec, output=output.spec, iRemoveLowest = 1){
  remove.lowestFreq <- function (spec, iRemove){
    {
      if (iRemove == 0)
        index = seq(spec$spec)
      else index <- (-(1:iRemove))
      spec$spec <- spec$spec[index]
      spec$freq <- spec$freq[index]
      spec$dof <- spec$dof[index]
      return(spec)
    }
  }
  get.fend.existing <- function (x){
    return(max(x$freq[!is.na(x$spec)]))
  }
  get.fstart.existing <- function (x) {
    return(min(x$freq[!is.na(x$spec)]))
  }
  get.df <- function (x){
    return(mean(diff(x$freq)))
  }
  AddConfInterval_Fdist <- function(transferspec, var.dof1, var.dof2, pval = 0.05){
    {
      if (!(length(transferspec$spec) == length(var.dof1)) && (!length(transferspec$spec) ==
                                                               length(var.dof2))) {
        stop("same lengths must be provided")
      }
      if (!(is.numeric(transferspec$spec)) || !(is.numeric(var.dof1)) ||
          !is.numeric(var.dof2)) {
        stop("non-numeric arguments")
      }
      res <- matrix(NA, nrow = length(transferspec$spec), ncol = 2)
      for (i in 1:length(transferspec$spec)) {
        QF <- qf(p = c(pval/2, (1 - pval/2)), df1 = var.dof1[i],
                 df2 = var.dof2[i])
        tmp <- QF * transferspec$spec[i]
        res[i, ] <- tmp
      }
      transferspec$lim.1 <- res[,2]
      transferspec$lim.2 <- res[,1]
      class(transferspec) <- "spec"
      return(transferspec)
    }
  }
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)
  freqRef <- seq(from = min(unlist(lapply(specList, get.fstart.existing))),
                 to = max(unlist(lapply(specList, get.fend.existing))),
                 by = min(unlist(lapply(specList, get.df))))
  specList.interpolated <- list()
  for (i in 1:length(specList)) specList.interpolated[[i]] <- SpecInterpolate(freqRef,
                                                                              specList[[i]])
  NSpectra <- length(specList.interpolated)
  result <- list(freq = specList.interpolated[[1]]$freq, spec = rep(0,
                                                                    length(specList.interpolated[[1]]$spec)))
  specMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  dofMatrix <- matrix(NA, NSpectra, length(specList.interpolated[[1]]$spec))
  for (i in 1:length(specList.interpolated)) {
    if (sum((result$freq - specList.interpolated[[i]]$freq)^2) >
        0.1)
      stop("Different spectra length or resolutions")
    specMatrix[i, ] <- specList.interpolated[[i]]$spec
    dofMatrix[i, ] <- specList.interpolated[[i]]$dof
  }

  var.dof1 <- dofMatrix[output, ]
  var.dof2 <- dofMatrix[input,]
  result$spec <- mapply('/', specMatrix[output, ], specMatrix[input, ]) #na.rm=TRUE
  result <- AddConfInterval_Fdist(result, dofMatrix[output, ], dofMatrix[input,])
  class(result) <- "spec"
  return(list(spec = result))
}

#--------------------------------------------------------------#
#' @title emulation of forced and internal variance
#' @description Wrapper for sampling internal variance based on ClimBayes::gen_noise_from_ebm_fit()
#' @param ebm_fit object of class ebm_fit from ClimBayes package
#' @param df.log decibels for logsmoothing in PaleoSpec::LogSmooth()
#' @param median defines whether mean or median is computet
#' @param probs.vec vector defining the quantile output
#' @param n_samples number of bootstrapping samples
#' @param debug plots spectra if TRUE
#' @return tibble with nested data that contains spectrum
#' @export
noisy_spectra <- function(ebm_fit, df.log=0.04, median=F, probs.vec=c(0.05, 0.5, 0.95), n_samples=1000, debug=F){
  transfer <- list()
  l <- list()
  tbb <- tibble()

  noise <- ClimBayes::gen_noise_from_ebm_fit(n_samples, ebm_fit)

  spec_list <- lapply(1:n_samples, function(j){
    noise_real <- as.vector(noise[j, ])
    full <- ebm_fit$posteriors$model_fit$median + noise_real
    if(debug){
      plot(full, type = "l")
      lines(ebm_fit$posteriors$model_fit$median, col = "red")
    }
    colnames(full) <- NULL
    sp <- PaleoSpec::SpecMTM(
      ts(full, deltat = 1))
    sp <- PaleoSpec::LogSmooth(sp, df.log=df.log)
    return(sp)
  }
  )

  #----quantiles----#
  spec.full <- simplify2array(lapply(spec_list, function(x) x$spec))
  q <- function(x){quantile(x, probs=probs.vec, na.rm=T)}
  spec.range.q <- apply(spec.full, 1, q)
  spec.range.q <- as_tibble(t(spec.range.q), name_repair=T) %>% add_column(freq=spec_list[[1]]$freq)
  #-----------------#
  if(median){
    l[["fit+noise"]] <- list(spec=spec.range.q$`50%`, freq=spec.range.q$freq)
  } else{
    l[["fit+noise"]] <- MeanSpec(spec_list, iRemoveLowest = 0)$spec
  }
  tbb <- rbind(tbb, as_tibble(cbind(list_to_tibble(l) %>% unnest(data), spec.range.q[!names(spec.range.q) %in% c("50%", "freq")]) %>% group_by(name) %>% nest())) %>% add_column(process="spec")

  l[["simulation"]] <- PaleoSpec::SpecMTM(ts(ebm_fit$input_params$y_obs, deltat = 1))
  l[["simulation"]] <- PaleoSpec::LogSmooth(l[["simulation"]], df.log=df.log)

  l[["fit"]] <- PaleoSpec::SpecMTM(ts(ebm_fit$posterior$model_fit$median, deltat = 1))
  l[["fit"]] <- PaleoSpec::LogSmooth(l[["fit"]], df.log=df.log)

  transfer[["fit+noise"]] <- transferSpec(l, input=2, output=1, iRemoveLowest = 1)$spec
  transfer[["fit"]] <- transferSpec(l, input=2, output=3, iRemoveLowest = 1)$spec

  l[["fit+noise"]] <- NULL

  tbb <- rbind(tbb, as_tibble(list_to_tibble(l)) %>% add_column(process="spec"), as_tibble(list_to_tibble(transfer)) %>% add_column(process="transfer"))

  if(debug){
    show(LPlot(tbb %>% filter(process=="spec", name=="fit+noise") %>% unnest(data), ylim = c(1e-4, 10), col = "blue"))
    show(LLines(tbb %>% filter(process=="spec", name=="simulation") %>% unnest(data), col="black"))
    show(LLines(tbb %>% filter(process=="spec", name=="fit") %>% unnest(data), col = "red"))
  }
  return(tbb)
}

#----------------------------#
#' @title Variance computation by integration of the spectrum
#' @description Simple computation of variance on frequency intervals by integration of the spectrum, akin to PaleoSpec::GetVarFromSpectra(),
#' with consideration of degrees of freedom, taken from Ellerhoff et al. (2022) (https://doi.org/10.5281/zenodo.6074747)
#' @param specs list of objects of class "spec"
#' @param f target timescales
#' @param dfreq default NULL
#' @param output
#' @return  returns list with variance and degrees of freedom
#' @export
getvaranddof_from_spec <- function(spec, f, dfreq = NULL, output="var"){
  out <- PaleoSpec::GetVarFromSpectra(spec, f, dfreq = dfreq)

  var.out <- out$var
  dof.out <- out$dof

  if(output=="var"){
    res <- var.out
  }
  if(output=="dof"){
    res <- dof.out
  }
  if(output=="all"){
    res =list(var=as.numeric(var.out), dof=as.numeric(dof.out))
  }
  return(res)
}

#--------------------------------------------------------------#
#' @title Confidence interval for variance ratios from F-distribution
#' @description Confidence interval for variance ratios computed from a F-distribution (see *PaleoSpec*-package) for specific p values
#' @param varratio variance ratio estimate
#' @param df.1 degrees of freedom of denominator
#' @param df.2 degrees of freedom of numerator
#' @param pval p-value, e.g. 0.05 for 95% confidence interval
#' @return vector of quantiles
#' @export
ConfRatio <- function (varratio, df.1, df.2, pval = 0.05)
{
  return(varratio * qf(c(pval/2, (1 - pval/2)), df1 = df.1,
                       df2 = df.2))
}

#--------------------------------------------------------------#
#' @title Variance computation over tibble
#' @description Wrapper for variance computation over tibble using GetVar()
#' @param tibble input data object of class tibble
#' @param tscales list with vectors c(f1,f2) that contains timescales for variance computation
#' @param output
#' @return object of class tibble()
#' @export
tscale_var <- function(tibble,
                       tscales=list("interannual"=list("t1"=2, "t2"=5),
                                            "decadal"=list("t1"=5,"t2"=20),
                                            "multidecadal"=list("t1"=20,"t2"=100),
                                            "centennial"=list("t1"=100,"t2"=200))){

  tscales <- list_to_tibble(tscales) %>% unnest(data) %>% rename(tscale=name)

  tmp <- tibble()
    for(i in seq_along(tscales$tscale)){
      freq.end <- tscales$t1[[i]]
      freq.start <- tscales$t2[[i]]
      tmp <- rbind(tmp, tibble %>% mutate(var = purrr::map(Spec, getvaranddof_from_spec, f=c(1/tscales$t2[[i]], 1/tscales$t1[[i]]), output="var"),
                                          tscale = tscales$tscale[[i]]) %>% select(-Spec) %>% ungroup())
    }

  tmp2 <- tibble()
    for(i in seq_along(tscales$tscale)){
      freq.end <- tscales$t1[[i]]
      freq.start <- tscales$t2[[i]]
      tmp2 <- rbind(tmp2, tibble %>% mutate(dof = purrr::map(Spec, getvaranddof_from_spec, f=c(1/tscales$t2[[i]], 1/tscales$t1[[i]]), output="dof"),
                                          tscale = tscales$tscale[[i]]) %>% select(-Spec) %>% ungroup())
    }

  tmp <- full_join(tmp, tscales)
  tmp2 <- full_join(tmp2, tscales)

  return(left_join(tmp, tmp2))
}

#-----------external helpers------------#
#taken from https://www.geeksforgeeks.org/add-common-legend-to-combined-ggplot2-plots-in-r/
get_only_legend <- function(plot) {

  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot))

  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")

  # extract legend
  legend <- plot_table$grobs[[legend_plot]]

  # return legend
  return(legend)
}

#taken from: https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
