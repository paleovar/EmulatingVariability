source("init.R")

dir.runs <- "output/cmip-emic-fits/" #source data directory
save <- F #TRUE for saving output as .csv file in the project's main directory

for(modeltype in c("CMIP", "EMIC", "EMIC_LO")){

  filename = paste0("T2_params_", modeltype, "_2box.csv") #name of output file

  #load data
  files_1box = list.files(dir.runs, pattern="_1box.rda")
  files_2box = list.files(dir.runs, pattern="_2box.rda")

  runs <- c(signal_tbb %>% filter(type==modeltype))$name
  tbb <- tibble()
  sds <- tibble()
  for (i in seq_along(runs)){
      print(runs[[i]])
      if(!is.null(grep(paste0(runs[[i]]), files_2box))){
        res_cmip <-  loadRData(paste0(dir.runs, files_2box)[[grep(paste0(runs[[i]]), files_2box)]])
        tbb <- rbind(tbb, list_to_tibble(res_cmip$posteriors$parameters) %>% unnest(data) %>%
                       add_column(run=runs[[i]], box="2box"))
        sds <- rbind(sds, tibble(name=runs[[i]],
                                 sd=round(sd(res_cmip$posteriors$model_fit$mean - res_cmip$input_params$y_obs),2)))
    }
  }

  #create table
  table_tbb <- list()
  for(i in unique(tbb$name)){
    table_tbb[[i]] <- tbb %>% filter(name==i) %>% mutate(
      !!i:=paste0(as.character(format(round(mean,2), nsmall = 2)), " (",
                  as.character(format(round(lower_quant,2), nsmall = 2)), ",",
                  as.character(format(round(upper_quant,2), nsmall = 2)), ")")
    ) %>% select(-median, -mean, -variance, -name, -lower_quant, -upper_quant)
  }

  table_tbb <- purrr::reduce(table_tbb, left_join) %>% rename(name=run) %>% left_join(., signal_tbb) %>%
    relocate(alt_name, type, box, lambda1, lambda2, F0, weights1) %>%
    left_join(sds)

  table_tbb_nice <- table_tbb %>% filter(box=="2box") %>%
    select(-name, -box, -type) %>%
    relocate(alt_name, lambda1, lambda2, weights1, T0, F0, sd) %>%
    arrange(alt_name) %>%
    mutate(sd = format(round(sd, 2), nsmall = 2))

  #save tabble
  if(save){
    write.table(table_tbb_nice,
               filename, row.names = FALSE, quote=F,  na="", sep = ";", dec = ".")
  }
}

