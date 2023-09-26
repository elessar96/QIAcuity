
#' @title Load necessary packages
#' @description This function loads all necessary packages for the script to run.
#' @return NULL

setup <- function(){

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  if(!require("dplyr", quietly=TRUE))
    install.packages("dplyr")

  if(!require("openxlsx", quietly=TRUE))
    install.packages("openxlsx")

  if(!require("tidyr", quietly=TRUE))
    install.packages("tidyr")

  if(!require("tidyverse", quietly=TRUE))
    install.packages("tidyverse")

  if(!require("changepoint", quietly = TRUE))
    install.packages("changepoint")

  if(!require("pastecs", quietly=TRUE))
    install.packages("pastecs")

  library(pastecs)

  library(tidyr)
  library(tidyverse)
  library(ggplot2)
  library(dplyr)
  library(openxlsx)

  return(TRUE)
}



#' @title Analyze output of a QIAcuity experiment
#' @description This function takes the files exported from the QIAcuity software for one or more experiments and analyzes them in an automated way.
#' @param input_path path to where the .zip or .csv files exported from the QIAcuity software are saved. No files are written here. The function can handle multiple experiments saved in the same directory. It will always analyze all of them and create separate results for them.
#' @param output_path path to where any results and plots should be saved. Directory for each experiment will be created there.
#' @return A list containing a list of dataframes for each experiment.
#' @details This function returns a list with an element for each experiment that was in the directory set by path. In each element there are the following data frames: The data for all channels in raw data, after baseline correction, cross talk correction and competition correction; cross talk and competition analyses that give an estimate of the magnitude of the effects the different reactions have on each other; thresholds gives the thresholds for each channel at each step; Volumes contains information about the total volume contained in the partitions of each well.


QIAcuityAnalysis <- function(input_path, output_path, smooth=FALSE){
  if(setup_conducted ==FALSE){
    setup_conducted <- setup()
  }
  setwd(input_path)

  files <- dir()
  # find all experiments in directory
  plates <- gregexec("^.*_RFU_img.*", dir()) %>% regmatches(dir(), .) %>% unlist() %>% gsub("\\.\\w+", "",. ) %>% gsub("RFU_img\\d_[[:alpha:]]+_.+", "", .)  %>% unique()# %>% gsub("_RFU_img", "", .) %>% unique()
  plates <- plates[!plates=="character(0)"]

  rem <- c()
  for(plate in 1:length(plates)){
    files <- dir()[grep(plates[[plate]], dir())]
    files <- files[grep(".zip", files)]
    if(length(files) > 1){
      rem <- c(rem, plate)
      for(file in 1:length(files)){
        file.rename(files[[file]], paste("plate", file, files[[file]], sep="_"))
        plates <- c(plates, paste("plate", file, plates[[plate]], sep="_"))
      }
    }
  }
  if(length(rem)>0){
    plates <- plates[-rem]
  }

  plate_results <- rep(NA, times=length(plates)) %>% as.list()
  names(plate_results) <- plates

  setwd(input_path)
  writeLines("QIAcuity data analysis started", con=paste(input_path, "Log.txt", sep="\\"))

  # for each experiment, do:
  for(plate in 1:length(plates)){
    setwd(output_path)

    readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c(plates[[plate]])) %>%
      writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
    setwd(input_path)

    # find all files with experiment name
    plate_files <- dir()[grep(plates[[plate]], dir())]
    #if csv files among them, analyze those files (assumption: .zip has been unzipped!)
    csv_files <- plate_files[grep(".csv", plate_files)]

    # if no csv files, unzip any zip-files
    if(length(csv_files)==0){
      zip_file <- plate_files[grep(".zip", plate_files)]
      if(length(zip_file)==0){
        stop("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.")
        writeLines("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.\n", con=paste(input_path, "Log.txt", sep="\\"))
      }
      dirname <- paste(output_path, "\\plate ", plate, "_", plates[[plate]], sep="")
      dir.create(dirname)
      unzip(zip_file, exdir=getwd(), overwrite = FALSE)

      setwd(dirname)
      csv_files <- dir()[grep(paste(plates[[plate]] %>% gsub("plate_\\d_", "", .), ".+", ".csv",sep=""), dir())]
      data_list <- csv_files %>% lapply(., read.data.QIAcuity)

    }else{
      data_list <- csv_files %>% lapply(., read.data.QIAcuity)
    }

    # gather data and find baseline of all samples, calculate corrected RFU values
    raw_data <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))

    for(channel in 1:length(data_list)){
      #extract RFU data from raw csv files
      channel_data <- data_list[[channel]]
      channel_data <- channel_data %>% data.frame() %>% mutate(RFU=RFU %>% as.numeric()) %>%
                          group_by(Well) %>% mutate(Partition=Partition %>% as.numeric())

      #put those RFU values into one dataframe (that contains the data of all channels)
      raw_data <- channel_data %>% select(Well, Sample, Partition, RFU) %>% merge(., raw_data, all=TRUE, by=c("Well", "Sample", "Partition"))
      colnames(raw_data)[which(colnames(raw_data)=="RFU")] <- channel_data %>% pull(Channel) %>% unique()
    }

    channels <- raw_data %>% select(!Well) %>% select(!Partition) %>% select(!Sample) %>% colnames()

    coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0))

    wells <- raw_data %>% pull(Well) %>% unique()

    thresholds <- data.frame(channel=channels, baseline=numeric(length=length(channels)), crosstalk = numeric(length=length(channels)), competition=numeric(length=length(channels)))

    maxima <- raw_data %>% group_by(Well) %>% na.omit() %>% summarize_at(channels, get_quantile, probs=0.999)

    pc_wells <- maxima %>% mutate(across(.cols=all_of(channels), .fns=function(x){return(x-min(x))})) %>%
                    mutate(across(.cols=all_of(channels), .fns=function(x){return(x/max(x))}))  %>% mutate(minimum = apply(across(all_of(channels)), 1, min)) %>%
                      mutate(median = apply(across(all_of(channels)), 1, median)) %>% filter(median == max(median) | minimum == max(minimum) | median == min(median)) %>% pull(Well) %>% unique()

    # do baseline correction for subtraction of background signal and removal of artifacts due to position effects
    baseline_corrected <- try(baseline_correction(raw_data=raw_data, smooth=smooth, channels=channels, coupled_channels = coupled_channels, pc_wells=pc_wells))

    if(class(baseline_corrected) =="try-error"){
      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Error: baseline correction failed. Message:", geterrmessage())) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
      baseline_corrected <- raw_data
    }

    #check for correlation between channels to prevent later errors!
    maxima_baseline_c <- baseline_corrected %>% group_by(Well) %>% na.omit() %>% summarize_at(channels, get_quantile, probs=0.999)
    tp_data <- baseline_corrected %>% filter(Well %in% pc_wells)
    maxima_channels <- tp_data %>% select(any_of(channels)) %>% apply(., 2, max, na.rm=TRUE)
    sd_channels <-  tp_data %>% select(any_of(channels)) %>% apply(., 2, sd, na.rm=TRUE)
    for(channel in 1:length(channels)){

      tp <- tp_data %>% select(any_of(channels))  %>% find_turnpoints(., variable=channels[[channel]])
      classified_turnpoints <- classify_peaks(tp, channel_maxima = maxima_channels , variable=channels[[channel]], intensities =tp_data %>% select(any_of(channels)))
      true_pos <- classified_turnpoints %>% filter(tp.peak ==TRUE)
      if(nrow(true_pos)==0){
        true_pos <- classified_turnpoints %>% filter(tp.peaks == TRUE) %>% filter(d.x == max(d.x))
      }
      lower_lim <- classified_turnpoints %>% filter(d.x<true_pos %>% pull(d.x)) %>% filter(d.x==max(d.x)) %>% pull(d.x)
      upper_lim <- classified_turnpoints %>% filter(d.x>true_pos %>% pull(d.x)) %>% filter(d.x==min(d.x)) %>% pull(d.x)

      peak_data <- tp_data %>% filter(!!sym(channels[[channel]])>lower_lim) %>% filter(!!sym(channels[[channel]])<upper_lim) %>%
                    summarize_at(channels[-channel], median)

      peak_data <- peak_data + sd_channels[-channel] -maxima_channels[-channel]
      peak_data <- peak_data[which(peak_data > 0)]
      if(length(peak_data)>0){
        for(peak in 1:length(peak_data)){
          pair <- data.frame(ch1=c(channels[channel], names(peak_data)[peak]), ch2=c(names(peak_data)[peak], channels[channel]))
          coupled_channels <- coupled_channels %>% rbind(., pair)
        }
      }
    }

    coupled_channels <- coupled_channels %>% unique()

    #define thresholds for each channel that separate positive and negative partitions
    temp <- recalculate_thresholds(data=baseline_corrected, thresholds=thresholds, step="baseline", coupled_channels = coupled_channels)

    baseline_corrected_dichot <- temp[[1]]
    thresholds <- temp[[2]]

    # crosstalk correction using a linear model
    temp <- try(crosstalk_correction(baseline_data = baseline_corrected, baseline_data_dichot = baseline_corrected_dichot, pc_wells = pc_wells, coupled_channels = coupled_channels, channels=channels, thresholds=thresholds))

    if(class(temp) =="try-error"){
      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Error: crosstalk correction failed. Message:", geterrmessage())) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
      crosstalk_corrected <- baseline_corrected
      crosstalk_analysis <- ""
    }else{
      crosstalk_corrected <- temp[[1]]
      crosstalk_analysis <- temp[[2]]
    }

    # calculate thresholds for new data
    temp <- recalculate_thresholds(data=crosstalk_corrected, thresholds=thresholds, step="crosstalk", coupled_channels = coupled_channels)
    crosstalk_corrected_dichot <- temp[[1]]
    thresholds <- temp[[2]]

    # perform competition correction
    temp <- try(competition_correction(crosstalk_corrected=crosstalk_corrected, channels=channels, pc_wells=pc_wells, thresholds=thresholds))

    if(class(temp) =="try-error"){
      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Error: competition correction failed. Message:", geterrmessage())) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
      competition_corrected <- crosstalk_corrected
      competition_analysis <- ""
    }else{
      competition_corrected <- temp[[2]]
      competition_analysis <- temp[[1]]
    }

    temp <- recalculate_thresholds(data=competition_corrected, thresholds=thresholds, step="competition", coupled_channels = coupled_channels, min_dist=0.5)
    competition_corrected_dichot <- temp[[1]]
    thresholds <- temp[[2]]

    snr_vec <- c()

    # calculate SNRs for each channel
    for(curr_channel in 1:length(channels)){
      pos <- competition_corrected %>% filter(!!sym(channels[[curr_channel]]) >= thresholds %>% filter(channel== channels[[curr_channel]]) %>% pull("competition")) %>% pull(channels[[curr_channel]])
      neg <- competition_corrected %>% filter(!!sym(channels[[curr_channel]]) < thresholds %>% filter(channel== channels[[curr_channel]]) %>% pull("competition")) %>% pull(channels[[curr_channel]])

      signal <- median(pos) - median(neg)
      noise <- mean(c(rep(sd(pos), times=length(pos)), rep(sd(neg), times=length(neg))))*2 + median(neg)

      SNR <- signal^2/noise^2
      snr_vec[[curr_channel]] <- SNR
      names(snr_vec)[[curr_channel]] <- channels[[curr_channel]]
    }

    #if(readLines("../Log.txt") %>% unlist() %>% paste(collapse = "") %>% nchar() == 0 ){
    readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Analysis completed.")) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
    #}

    output <- list(raw_data = raw_data, baseline_corrected = baseline_corrected, crosstalk_corrected=crosstalk_corrected, competition_corrected=competition_corrected, crosstalk_analysis=crosstalk_analysis, competition_analysis=competition_analysis, thresholds=thresholds)

    # calculate volume

    vol_by_well = data_list %>% bind_rows() %>% mutate(`Is invalid` = ifelse(`Is invalid` == "0", "valid", "invalid")) %>% group_by(Well, `Cycled volume`, Channel, `Is invalid`) %>%
                    summarize(N_total = n()) %>% pivot_wider(names_from= `Is invalid`, values_from=N_total) %>% mutate(`Cycled volume` = as.numeric(`Cycled volume`)*(valid/(invalid + valid)))

    output$Volumes <- vol_by_well
    output$snr <- snr_vec
    plate_results[[plates[[plate]]]] <- output
  }

  return(plate_results)
}


#' @title Generate output of QIAcuity Analyisis
#' @description  This function summarizes and visualizes the results of one or multiple QIAcuity experiments as.
#' @param results The output of the QIAcuityAnalysis function.
#' @param detailed_plots Set to TRUE to generate 2D scatterplots at every correction step
#' @param scatterplots_1d Set to TRUE to generate 1D scatterplots at every correction step
#' @param output_path path to where files and plots should be saved. Will create directory for the experiment there.
#' @return NULL. Generates output files: .csv files for the results per partition, .xlsx files for summaries at each correction step. 2D scatterplots for all channel combinations and 1D scatterplots for each channel at every step. Also puts the crosstalk and competition analyses into .xlsx files.

summarizeQIAcuityResult <- function(results, detailed_plots=TRUE, scatterplots_1d=TRUE, output_path=getwd()){
  plate_summaries <- rep(NA, times=length(results)) %>% as.list()
  names(plate_summaries) <- names(results)
  corrections = c("raw_data", "baseline_corrected", "crosstalk_corrected", "competition_corrected")

  setwd(output_path)

  for(plate in 1:length(results)){
    dirname <- paste(output_path, "\\plate ", plate, "_", names(results)[[plate]], sep="")
    dir.create(dirname)
    setwd(dirname)
    n_channels <- results[[plate]]$thresholds %>% nrow()

    for(correction_step in 1:length(corrections)){
      assign("last.warning", NULL, envir = baseenv())
      well_order <- results[[plate]][[correction_step]] %>% select(Well) %>% unique() %>%
        mutate(letter=Well %>% gsub("[[:digit:]]", "", .)) %>% mutate(digit=Well %>%
            gsub("[[:alpha:]]", "", .)) %>% arrange(digit, letter) %>% pull(Well)

      write.csv(results[[plate]][[correction_step]], file=paste(corrections[[correction_step]], ".csv", sep=""))

      channel_summaries <- data.frame(Sample=character(length=0), Well=character(length=0),
                                      negative=numeric(length = 0), positive = numeric(length = 0),
                                      Channel = character(length=0), poisson_corrected_targets = numeric(length=0))
      for(channel1 in 1:n_channels){
        if(detailed_plots==TRUE){
          for(channel2 in 1:n_channels){
            if(!channel1>=channel2){
              if(n_channels>2){
                p = channel2-1
                if(p==channel1){
                  p = channel2-2
                }
                if(p==0){
                  p = channel2+1
                }
                scatterplot <- results[[plate]][[corrections[[correction_step]]]] %>%
                  dplyr::select(Well | ends_with(results[[plate]]$thresholds$channel[c(channel2, channel1, p)])) %>%
                     na.omit() %>% ggplot(aes(x =!!sym(results[[plate]]$thresholds$channel[[channel2]]), y =!!sym(results[[plate]]$thresholds$channel[[channel1]]), color=!!sym(results[[plate]]$thresholds$channel[[p]]))) + geom_point()
                if(!correction_step==1){
                  scatterplot <- scatterplot + geom_hline(yintercept=results[[plate]]$thresholds[channel1, corrections[[correction_step]] %>%
                                                            gsub("_corrected", "", .)], color="red") + geom_vline(xintercept=results[[plate]]$thresholds[channel2, corrections[[correction_step]] %>% gsub("_corrected", "", .)], color="red")

                }
                plotname <- paste(names(plate_summaries)[[plate]], corrections[[correction_step]], "2dscatter_", paste(results[[plate]]$thresholds$channel[c(channel2, channel1, p)], collapse=""), ".png", sep="" )
                ggsave(scatterplot, filename = plotname)
              }
            }
          }
        }
        if(scatterplots_1d==TRUE){
          plot <- results[[plate]][[corrections[[correction_step]]]]  %>%
            dplyr::select(Well | Partition | ends_with(results[[plate]]$thresholds$channel[channel1])) %>% na.omit() %>% mutate(Partition = Partition %>% as.numeric()) %>% mutate(across(Well, ~ factor(., levels=c(well_order)))) %>% ggplot(aes(x=Partition, y =!!sym(results[[plate]]$thresholds$channel[[channel1]]))) + geom_jitter(size=0.5) + scale_x_discrete(limits=well_order) + facet_wrap(vars(Well), nrow =2) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

          wells <- results[[plate]][[corrections[[correction_step]]]] %>% pull(Well) %>% unique()

          if(!correction_step==1){
            plot <- plot + geom_hline(yintercept=results[[plate]]$thresholds[channel1, corrections[[correction_step]] %>% gsub("_corrected", "", .)], color="red")
          }

          plotname <- paste(names(plate_summaries)[[plate]], "_", corrections[[correction_step]], "_", results[[plate]]$thresholds$channel[[channel1]], ".png", sep="")
          ggsave(plot, filename = plotname, width=14, height=7)
        }


        temp_res <- results
        if(correction_step>1){
          temp_res[[plate]][[corrections[[correction_step]]]]  <- temp_res[[plate]][[corrections[[correction_step]]]] %>% mutate(new = ifelse(temp_res[[plate]][[corrections[[correction_step]]]][temp_res[[plate]]$thresholds$channel[[channel1]]]>temp_res[[plate]]$thresholds[channel1, corrections[[correction_step]] %>% gsub("_corrected", "", .)], 1, 0)) %>% select(!temp_res[[plate]]$thresholds$channel[[channel1]]) %>% rename(!!temp_res[[plate]]$thresholds$channel[[channel1]] := new)

          negs <- temp_res[[plate]][[corrections[[correction_step]]]] %>% filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])==0) %>%
            group_by(Well, Sample) %>% count()
          pos <-  temp_res[[plate]][[corrections[[correction_step]]]] %>% filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])==1) %>%
            group_by(Well, Sample) %>% count()

          summary <- merge(negs, pos, all=TRUE, by=c("Well", "Sample")) %>% rename(positive=n.y, negative=n.x) %>%
            mutate(Channel=temp_res[[plate]]$thresholds$channel[[channel1]])

          summary <- summary %>% mutate(positive= replace(positive, is.na(positive), 0)) %>%
            mutate(poisson_corrected_targets = (-1* (positive + negative)*log(negative/(positive+negative))) %>% round(digits=2))

          mean_RFU_pos <- results[[plate]][[corrections[[correction_step]]]] %>%
            filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])>results[[plate]]$thresholds[channel1, corrections[[correction_step]] %>% gsub("_corrected", "", .)]) %>% pull(temp_res[[plate]]$thresholds$channel[[channel1]]) %>% mean()
          sd_RFU_pos <- results[[plate]][[corrections[[correction_step]]]] %>%
            filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])>results[[plate]]$thresholds[channel1, corrections[[correction_step]] %>% gsub("_corrected", "", .)]) %>% pull(temp_res[[plate]]$thresholds$channel[[channel1]]) %>% sd()

          zscore_RFU_pos <- results[[plate]][[corrections[[correction_step]]]] %>%
                                filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])>results[[plate]]$thresholds[channel1, corrections[[correction_step]] %>%
                                      gsub("_corrected", "", .)]) %>% mutate(z_score = (!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])-mean_RFU_pos)/sd_RFU_pos) %>%
                                          group_by(Well) %>% summarize_at("z_score", mean)
          zscore_RFU_pos <- zscore_RFU_pos %>% mutate(Warning = ifelse(abs(z_score)>2, "Warning: mean RFU value of positive partitions further than 2 standard deviations away from global mean. Potential artifact.", NA))

          summary <- full_join(summary, zscore_RFU_pos) %>% rename("mean_z_score"=z_score)
          channel_summaries <- bind_rows(channel_summaries, summary)
        }
      }
      if(!correction_step==1){

        channel_summaries <- channel_summaries %>% mutate(positive=replace(positive, is.na(positive), 0)) %>%
                                mutate(negative=replace(negative, is.na(negative), 0)) %>%
                                  mutate(valid=positive + negative)

        mm_vol <- 0

        n_part <- channel_summaries %>% pull(valid)

        if(n_part %>% median(., na.rm = TRUE) > 10000){
          mm_vol <- 40/26000
        }else{
          mm_vol <- 12/8500
        }

        poisson.CI <- function(x){return(poisson.test(x)$conf.int %>% round(digits=1) %>% paste(., collapse=" - "))}

        channel_summaries <- merge(channel_summaries, results[[plate]]$Volumes %>%
                                     select(Well, `Cycled volume`, Channel), by=c("Well", "Channel")) %>% mutate(factor = (!!sym("Cycled volume") %>% as.numeric())/(mm_vol*as.numeric(valid)))

        channel_summaries <- channel_summaries %>% mutate(CI=lapply(.$positive, poisson.CI)) %>%
                                    relocate(any_of(c("Sample", "Well", "Channel", "negative", "positive", "CI")))

        channel_summaries <- channel_summaries %>% mutate(Warning = ifelse(lapply(.$CI %>% as.character() %>% strsplit(split=" - "), min)<1 & .$positive > 0, "Warning: confidence interval contains values < 1, positive result may not be reliable.", NA) %>% paste(channel_summaries$Warning, sep=", ") %>% gsub("NA, NA", NA, .) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .))

        channel_summaries <- channel_summaries %>% mutate(poisson_corrected_targets = (-1* (positive/factor + negative/factor)*log(negative/factor/(positive/factor+negative/factor))) %>% round(digits=2)) %>% mutate()

        channel_summaries <- channel_summaries %>% arrange(factor(Well, levels = well_order)) %>% select(!factor) %>% select(!'Cycled volume')

        excel_output <- createWorkbook()
        addWorksheet(excel_output, "dPCR results")

        negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")

        writeData(excel_output, "dPCR results", channel_summaries)
        index_warning <- which(colnames(channel_summaries)=="Warning")
        conditionalFormatting(excel_output, "dPCR results",
                              cols = index_warning,
                              rows = 2:nrow(channel_summaries), type = "contains",
                              rule = "Warning", style = negStyle
        )

        saveWorkbook(excel_output, overwrite=TRUE, file = paste("result_", corrections[[correction_step]], ".xlsx", sep=""))
      }
    }
    write.xlsx(results[[plate]]$crosstalk_analysis, file="crosstalk_analysis.xlsx")
    write.xlsx(results[[plate]]$competition_analysis, file="competition_analysis.xlsx")
    write.xlsx(results[[plate]]$thresholds, file="thresholds.xlsx")
    write.xlsx(results[[plate]]$snr %>% data.frame(), file="snr_competition_corrected.xlsx")
    setwd("..")
  }
}
