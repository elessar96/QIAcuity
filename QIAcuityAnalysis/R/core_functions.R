
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

  if(!require("DescTools", quietly=TRUE))
    install.packages("DescTools")

  library(pastecs)
  library(DescTools)
  library(tidyr)
  library(tidyverse)
  library(ggplot2)

  library(openxlsx)
  library(gtools)

  library(dplyr)

  return(TRUE)
}



#' @title Analyze output of a QIAcuity experiment
#' @description This function takes the files exported from the QIAcuity software for one or more experiments and analyzes them in an automated way.
#' @param input_path path to where the .zip or .csv files exported from the QIAcuity software are saved. No files are written here. The function can handle multiple experiments saved in the same directory. It will always analyze all of them and create separate results for them.
#' @param output_path path to where any results and plots should be saved. Directory for each experiment will be created there.
#' @param noise_suppression_strict Whether groups of fewer than 10 partitions should be considered in the analysis or not. Set to FALSE if the only positive sample on a plate contains very few positive partitions, or to salvage results from a run without a positive sample. Default: TRUE
#' @return A list containing a list of dataframes for each experiment, which contain detailed information about the analysis, for instance the raw or processed data at each step of the analysis, thresholds used at each of the steps, estimates of crosstalk and competition strength etc.
#' @details This function returns a list with an element for each experiment that was in the directory set by path. In each element there are the following data frames: The data for all channels in raw data, after baseline correction, cross talk correction and competition correction; cross talk and competition analyses that give an estimate of the magnitude of the effects the different reactions have on each other; thresholds gives the thresholds for each channel at each step; Volumes contains information about the total volume contained in the partitions of each well.


QIAcuityAnalysis <- function(input_path,
                             output_path,
                             noise_suppression_strict = TRUE,
                             perform_crosstalk_correction = TRUE,
                             perform_competition_correction = TRUE){

  if(setup_conducted ==FALSE){
    setup_conducted <- setup()
  }

  files <- dir(path = input_path)
  # find all experiments in directory
  plates <- gregexec("^.*_RFU_img.*", dir(path = input_path)) %>% regmatches(dir(path = input_path), .) %>% unlist() %>% gsub("\\.csv|\\.xlsx|\\.zip", "",. ) %>% gsub("RFU_img\\d_[[:alpha:]]+_.+", "", .)  %>% unique()# %>% gsub("_RFU_img", "", .) %>% unique()
  plates <- plates[!plates=="character(0)"]

  rem <- c()
  for(plate in 1:length(plates)){
    files <- dir(path = input_path)[grep(plates[[plate]], dir(path = input_path))]
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


  writeLines("QIAcuity data analysis started", con=paste(input_path, "Log.txt", sep="\\"))

  # for each experiment, do:
  for(plate in 1:length(plates)){

    readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c(plates[[plate]])) %>%
      writeLines(., con=paste(input_path, "Log.txt", sep="\\"))

    # find all files with experiment name
    plate_files <- dir(path = input_path)[grep(plates[[plate]], dir(path = input_path))]

    plate_date <- plate_files %>%
      gregexec("[[:alpha:]]_\\d+_\\d+_\\d+_", .)
    plate_date <- plate_files %>%
      substr(., start = unlist(plate_date), stop = unlist(plate_date) + unlist(attr(plate_date[[1]], "match.length")) -2) %>%
      gsub("[[:alpha:]]", "", .) %>%
      gsub("^_", "", .) %>%
      gsub("_", "-", .) %>%
      as.Date(., format = "%d-%m-%Y")

    #if csv files among them, analyze those files (assumption: .zip has been unzipped!)
    csv_files <- plate_files[grep(".csv", plate_files)]
    csv_files <- csv_files[grep("RFU_img\\d+_[[:alpha:]]+_\\d+", csv_files)]

    # if no csv files, unzip any zip-files
    if(length(csv_files)==0){
      zip_file <- plate_files[grep(".zip", plate_files)]
      if(length(zip_file)==0){
        stop("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.")
        writeLines("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.\n", con=paste(input_path, "Log.txt", sep="\\"))
      }
      dirname <- paste(output_path, "\\plate ", plate, "_", plates[[plate]], sep="")
      dir.create(dirname)
      unzip(paste(input_path, "\\", zip_file, sep=""), exdir=dirname, overwrite = FALSE)

      csv_files <- dir(path = dirname)[grep(paste(plates[[plate]] %>% gsub("plate_\\d_", "", .), ".+", ".csv",sep=""), dir(path = dirname))]
      data_list <- paste(dirname, "\\", csv_files, sep="") %>% lapply(., read.data.QIAcuity)

    }else{
      data_list <- paste(output_path, "\\", csv_files, sep="") %>% lapply(., read.data.QIAcuity)
    }

    # gather data and find baseline of all samples, calculate corrected RFU values
    raw_data <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))

    reference_data <- NULL

    for(channel in 1:length(data_list)){
      #extract RFU data from raw csv files
      channel_data <- data_list[[channel]]
      channel_data <- channel_data %>% data.frame() %>% mutate(RFU=RFU %>% as.numeric()) %>%
                          group_by(Well) %>% mutate(Partition=Partition %>% as.numeric())
      if(!"Channel" %in% colnames(channel_data)){
        reference_data <- channel_data
      }else{
        #put those RFU values into one dataframe (that contains the data of all channels)
        raw_data <- channel_data %>%
          select(Well, Sample, Partition, RFU) %>%
          merge(., raw_data, all=TRUE, by=c("Well", "Sample", "Partition"))
        colnames(raw_data)[which(colnames(raw_data)=="RFU")] <- channel_data %>% pull(Channel) %>% unique()
      }

    }

    channels <- raw_data %>% select(!Well) %>% select(!Partition) %>% select(!Sample) %>% colnames()

    coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0))

    wells <- raw_data %>% pull(Well) %>% unique()

    thresholds <- data.frame(channel=channels, raw=numeric(length=length(channels)), baseline=numeric(length=length(channels)), crosstalk = numeric(length=length(channels)), competition=numeric(length=length(channels)))

    maxima <- raw_data %>% group_by(Well) %>% na.omit() %>% summarize_at(channels, get_quantile, probs=0.9999)

    rel_maxima <- maxima %>% mutate(across(.cols=all_of(channels), .fns=function(x){return(x-min(x))})) %>%
      mutate(across(.cols=all_of(channels), .fns=function(x){return(x/max(x))}))  %>% mutate(minimum = apply(across(all_of(channels)), 1, min)) %>%
      mutate(median = apply(across(all_of(channels)), 1, median)) %>%
      mutate(average = apply(across(all_of(channels)), 1, mean))

    pc_wells <- rel_maxima %>% filter(minimum == max(minimum)) %>% pull(Well) %>% unique()
    min_fluorescence <- rel_maxima %>%
      filter(Well %in% pc_wells) %>%
      select(any_of(channels)) %>%
      apply(., 2, max) %>%
      min()

    if(min_fluorescence < 0.6){

      maxima <- rel_maxima %>%
        filter(Well %in% pc_wells) %>%
        select(any_of(channels)) %>%
        t() %>%
        apply(., 1, max, na.rm = TRUE)

      channel_tests <- channels[which(maxima < 0.7)]

      for(channel_pc_detection in 1:length(channel_tests)){
        new_pc <- rel_maxima %>%
          filter(!!sym(channel_tests[[channel_pc_detection]]) == max(!!sym(channel_tests[[channel_pc_detection]]))) %>%
          pull(Well)
        pc_wells <- c(pc_wells, new_pc) %>% unique()
      }

      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append("Error: No valid positive control found. Positive wells for each channel will be searched on the plate to create composite data instead.") %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
    }

    raw_thresh <- try(recalculate_thresholds(raw_data, thresholds, step="raw", coupled_channels = coupled_channels, pc_wells = pc_wells, noise_suppression_strict = noise_suppression_strict))

    if(!class(raw_thresh) =="try-error"){
      thresholds <- raw_thresh[[2]]
      raw_data_dichot <- raw_thresh[[1]]
    }
  #check for link between channels to prevent later errors!

    pc_data <-  raw_data %>% filter(Well %in% pc_wells) %>% na.omit()

    for(channel1 in 1:length(channels)){
      for(channel2 in 1:length(channels)){
        if(!channel1==channel2){
          pos_rate_ch1 <- pc_data %>% filter(!!sym(channels[[channel1]]) > thresholds %>% filter(channel == channels[[channel1]]) %>% pull("raw"))
          pos_rate_ch1 = nrow(pos_rate_ch1)/nrow(pc_data)

          sp_rate_ch1 <- pc_data %>%
            filter(!!sym(channels[[channel1]]) > thresholds %>% filter(channel == channels[[channel1]]) %>% pull("raw")) %>%
            filter(!!sym(channels[[channel2]]) < thresholds %>% filter(channel == channels[[channel2]]) %>% pull("raw"))

          sp_rate_ch1 = nrow(sp_rate_ch1)/nrow(pc_data)

          pos_rate_ch2 <- pc_data %>% filter(!!sym(channels[[channel2]]) > thresholds %>% filter(channel == channels[[channel2]]) %>% pull("raw"))
          pos_rate_ch2 = nrow(pos_rate_ch2)/nrow(pc_data)

          sp_rate_ch2 <- pc_data %>%
            filter(!!sym(channels[[channel2]]) > thresholds %>% filter(channel == channels[[channel2]]) %>% pull("raw")) %>%
            filter(!!sym(channels[[channel1]]) < thresholds %>% filter(channel == channels[[channel1]]) %>% pull("raw"))

          sp_rate_ch2 = nrow(sp_rate_ch2)/nrow(pc_data)

          dp_rate <- pc_data %>% filter(!!sym(channels[[channel2]]) > thresholds %>% filter(channel == channels[[channel2]]) %>% pull("raw")) %>% filter(!!sym(channels[[channel1]]) > thresholds %>% filter(channel == channels[[channel1]]) %>% pull("raw"))
          dp_rate = nrow(dp_rate)/nrow(pc_data)

          max_expected_dp_rate = sqrt(pos_rate_ch1 * pos_rate_ch2)
          ch1ch2_corr <- CCC( x = pc_data %>% pull(!!sym(channels[[channel1]])),
                              y = pc_data %>% pull(!!sym(channels[[channel2]])))$rho.c$est %>% abs()
          if(dp_rate > max_expected_dp_rate | ch1ch2_corr > 0.5){
            coupled_channels <- coupled_channels %>% rbind(., data.frame(ch1=channels[[channel1]], ch2=channels[[channel2]]))
          }


        }
      }
    }
    coupled_channels <- coupled_channels %>% unique()

    #if(nrow(coupled_channels) > length(channels)^2)

    # do baseline correction for subtraction of background signal and removal of artifacts due to position effects
    baseline_corrected <- try(baseline_correction(raw_data=raw_data, smooth=TRUE, channels=channels, coupled_channels = coupled_channels, pc_wells=pc_wells, thresholds = thresholds))

    if(class(baseline_corrected) =="try-error"){
      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Error: baseline correction failed. Message:", geterrmessage())) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
      baseline_corrected <- raw_data
    }

    #define thresholds for each channel that separate positive and negative partitions
    temp <- recalculate_thresholds(data=baseline_corrected, thresholds=thresholds, step="baseline", coupled_channels = coupled_channels, pc_wells=pc_wells, noise_suppression_strict = noise_suppression_strict)

    baseline_corrected_dichot <- temp[[1]]
    thresholds <- temp[[2]]

    if(perform_crosstalk_correction){
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
    }else{
      readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Crosstalk correction skipped.")) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))

      crosstalk_corrected <- baseline_corrected
      crosstalk_analysis <- ""
    }

    # calculate thresholds for new data
    temp <- recalculate_thresholds(data=crosstalk_corrected, thresholds=thresholds, step="crosstalk", coupled_channels = coupled_channels, pc_wells=pc_wells, noise_suppression_strict = noise_suppression_strict)
    crosstalk_corrected_dichot <- temp[[1]]
    thresholds <- temp[[2]]

    # perform competition correction
    if(perform_competition_correction){
      temp <- try(competition_correction(crosstalk_corrected=crosstalk_corrected, channels=channels, pc_wells=pc_wells, thresholds=thresholds, step="competition"))

      if(class(temp) =="try-error"){
        readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Error: competition correction failed. Message:", geterrmessage())) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))
        competition_corrected <- crosstalk_corrected
        competition_analysis <- ""
      }else{
        competition_corrected <- temp[[2]]
        competition_analysis <- temp[[1]]
      }
    }else{
      competition_corrected <- crosstalk_corrected
      competition_analysis <- ""
    }
    temp <- recalculate_thresholds(data=competition_corrected, thresholds=thresholds, step="competition", coupled_channels = coupled_channels, min_dist=0.5, pc_wells=pc_wells, mid=TRUE)
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

    # find artifacts based on z-score in relation to negative partitions
    # identify partitions that deviate from baseline in multiple channels
    # join RFU values and thresholds
    RFU_thresh_joined <- competition_corrected %>%
      pivot_longer(., cols = channels, names_to = "channel", values_to = "RFU") %>%
      full_join(., thresholds %>%
                  select(channel, competition) %>%
                  dplyr::rename(Threshold = competition))
    # identify mean and sd values of negative partitions per well, then calculate z-score in relation to negative partitions for all partitions
    RFU_thresh_joined <- RFU_thresh_joined %>%
      filter(RFU < Threshold) %>%
      group_by(Well, channel) %>%
      summarize(mean_RFU = mean(RFU), sd_RFU = sd(RFU)) %>%
      full_join(RFU_thresh_joined) %>%
      mutate(z_score_neg = abs(RFU - mean_RFU)/sd_RFU)
    # summarize by finding median and minimum z-score across all channels for each partition
    z_score_neg <- RFU_thresh_joined %>%
      ungroup() %>%
      group_by(Well, Sample, Partition) %>%
      summarize(median_z_score = median(z_score_neg, na.rm = TRUE)) %>%
      mutate(median_z_score_outlier = ifelse(median_z_score > 4  & is.finite(median_z_score),1, 0))

    # find out which positive partitions are affected!
    z_score_neg_outliers <- RFU_thresh_joined %>%
      filter(RFU > Threshold) %>%
      right_join(., z_score_neg %>% filter(median_z_score_outlier ==1)) %>%
      filter(!is.na(channel)) %>%
      select(Well, Sample, Partition, channel, median_z_score)

    # find artefacts based on unlikely multiple occupancy!

    # identify channels that are coupled to each other even after crosstalk correction
    summary_multiple_occupancy_pc <- competition_corrected_dichot %>%
      filter(Well %in% pc_wells) %>%
      group_by(pick(channels)) %>%
      summarize(n = n()) %>%
      ungroup() %>%
      mutate(fraction = n/sum(n)) %>%
      rowwise() %>%
      mutate(sum_positives = sum(c_across(channels)))

    sp_fractions <- summary_multiple_occupancy_pc %>%
      filter(sum_positives == 1) %>%
      pivot_longer(cols = channels, names_to = "Channel", values_to = "status") %>%
      filter(status == 1) %>%
      select(Channel, fraction)
    if(length(channels)>1){
    permutations_channels <- gtools::permutations(n = length(channels), r = 2, v = channels)
    colnames(permutations_channels) <- c("channel1", "channel2")
    expected_frequencies <- permutations_channels %>%
      data.frame() %>%
      full_join(., sp_fractions %>% dplyr::rename(channel1 = "Channel", fracch1 = "fraction")) %>%
      full_join(., sp_fractions %>% dplyr::rename(channel2 = "Channel", fracch2 = "fraction")) %>%
      mutate(expected_val = fracch1 * fracch2) %>%
      select(channel1, channel2, expected_val) %>%
      mutate(status = 1) %>%
      pivot_wider(names_from = channel1, values_from = status) %>%
      filter(!is.na(expected_val)) %>%
      filter(!is.na(channel2)) %>%
      apply(., 1, function(x){x[x["channel2"]] <- 1
                              return(x)}) %>%
      t() %>%
      data.frame() %>%
      select(!channel2) %>%
      mutate(across(any_of(channels), function(x){replace(x, is.na(x), 0)})) %>%
      mutate(across(channels, function(x){x %>% trimws() %>% as.numeric() %>% return()})) %>%
      select(expected_val, any_of(channels))

    ratios_channel_combinations <- summary_multiple_occupancy_pc %>%
      filter(sum_positives == 2) %>%
      left_join(., expected_frequencies) %>%
      unique() %>%
      mutate(ratio = as.numeric(fraction)/as.numeric(expected_val)) %>%
      select(channels, ratio)

    ratios_channel_combinations <- ratios_channel_combinations %>%
      pivot_longer(., channels, names_to = "channel1", values_to = "status_ch1") %>%
      full_join(ratios_channel_combinations) %>%
      filter(status_ch1 == 1) %>%
      pivot_longer(., channels, names_to = "channel2", values_to = "status_ch2") %>%
      filter(status_ch2 != 0) %>%
      filter(channel1 != channel2) %>%
      select(channel1, channel2, ratio)

    # get # of valid partitions per well
    n_by_well <- data_list %>%
      bind_rows() %>%
      mutate(`Is invalid` = ifelse(`Is invalid` == "0", "valid", "invalid")) %>%
      group_by(Well, `Cycled volume`, Channel, `Is invalid`) %>%
      summarize(N_total = n()) %>%
      filter(`Is invalid` == "valid") %>%
      ungroup() %>%
      select(Well, Channel, N_total) %>%
      ungroup() %>%
      group_by(Well) %>%
      summarize(N_total = median(N_total))

    positive_rates <- competition_corrected_dichot %>%
      pivot_longer(., cols = all_of(channels), names_to ="channel") %>%
      filter(!is.na(value)) %>%
      group_by(Well, channel, value) %>%
      summarize(n = n()) %>%
      pivot_wider(., id_cols = c("Well", "channel"), names_from = "value", values_from = "n") %>%
      mutate(`1`= replace(`1`, is.na(`1`), 0)) %>%
      mutate(`0`= replace(`0`, is.na(`0`), 0)) %>%
      mutate(positive_rate = `1`/(`1` + `0`)) %>%
      select(!`1`) %>%
      select(!`0`)

    n_double_pos <- competition_corrected_dichot %>%
      group_by(across(all_of(c("Well", channels)))) %>% summarize(n=n()) %>%
      na.omit()

    multiple_occupancy <- data.frame(Well = character(), channel1 =character(), channel2 = character(), n =numeric())

    for(ch1 in 1:length(channels)){
      for(ch2 in 1:length(channels)){
        if(ch1 < ch2){
          n_double_pos_temp <- n_double_pos %>%
            ungroup() %>%
            filter(!!sym(channels[[ch1]]) == 1 & !!sym(channels[[ch2]]) == 1) %>%
            group_by(Well) %>%
            summarize(n = sum(n)) %>%
            mutate(channel1 = channels[[ch1]]) %>%
            mutate(channel2 = channels[[ch2]])
          multiple_occupancy <- multiple_occupancy %>% full_join(n_double_pos_temp)
        }
      }
    }

    multiple_occupancy <- multiple_occupancy

    perm <- permutations(n = length(channels), r = 2, v = channels)
    colnames(perm) <- c("channel1", "channel2")
    # TASK: SUMMARY OF NUMBER OF DOUBLE POSITIVE PARTITIONS ACROSS CHANNEL COMBINATIONS!
    test_likelihood <- perm %>%
      data.frame() %>%
      full_join(., positive_rates %>% dplyr::rename(channel1 = channel)) %>%
      full_join(ratios_channel_combinations) %>%
      mutate(ratio = replace(ratio, is.na(ratio), 1)) %>%
      full_join(., positive_rates %>% dplyr::rename(channel2 = channel) %>% dplyr::rename(positive_rate2 = positive_rate)) %>% full_join(n_by_well) %>%
      mutate(n_double_pos_expected = N_total * positive_rate * positive_rate2 * ratio) %>%
      full_join(multiple_occupancy, .) %>%
      filter(!is.na(n)) %>%
      mutate(p.value = (.) %>% mutate(across(everything(), as.numeric)) %>% apply(., MARGIN = 1, FUN = function(x){poisson.test(x = x[[4]] %>% round(), T=x[[8]], r = x[[5]] * x[[6]] * x[[7]], alternative = "greater")$p.value %>% return()}))

    associated_likelihoods <- competition_corrected_dichot %>%
      pivot_longer(., cols = all_of(channels), values_to = "status", names_to = "channel") %>%
      filter(!is.na(status)) %>%
      full_join(positive_rates) %>%
      mutate(negative_rate = 1 - positive_rate) %>%
      mutate(likelihood = positive_rate %>% replace(., which(status == 0), negative_rate[which(status == 0)])) %>%
      pivot_wider(., id_cols = c(Well, Partition), names_from = channel, values_from = likelihood)

    likelihoods <- associated_likelihoods %>%
      select(all_of(channels)) %>%
      mutate(likelihood = reduce(., `*`)) %>%
      select(likelihood) %>%
      cbind(associated_likelihoods, .)

    likelihoods <- likelihoods %>%
      full_join(., n_by_well, by = "Well") %>%
      mutate(expected_value = likelihood * N_total) %>%
      mutate(outlier = ifelse(expected_value < 0.05/24, TRUE, FALSE))

    multiple_occupancy_outliers <- likelihoods %>%
      filter(outlier) %>%
      select(Well, Partition, outlier, likelihood, expected_value) %>%
      left_join(., competition_corrected_dichot, by = c("Well", "Partition"))
  }
  else{
    multiple_occupancy_outliers <- data.frame()
  }
    # add info about positive partitions with only an aberrant z-score
    RFU_summary_pc_wells <- competition_corrected %>%
      filter(Well %in% pc_wells) %>%
      pivot_longer(., cols = channels, names_to = "channel", values_to = "RFU") %>%
      full_join(thresholds %>% select(channel, competition) %>% dplyr::rename(Threshold = "competition")) %>%
      filter(RFU > Threshold) %>%
      group_by(channel) %>%
      summarize(mean = mean(RFU, na.rm = TRUE), sd = sd(RFU, na.rm = TRUE))

    z_score_pos_outliers <- competition_corrected %>%
      pivot_longer(., cols = channels, names_to = "channel", values_to = "RFU") %>%
      full_join(thresholds %>% select(channel, competition) %>% dplyr::rename(Threshold = "competition")) %>%
      filter(RFU > Threshold) %>%
      full_join(RFU_summary_pc_wells) %>%
      mutate(z_score = (RFU - mean)/sd) %>%
      filter(abs(z_score) > 2)

    readLines(con=paste(input_path, "Log.txt", sep="\\")) %>% append(c("Analysis completed.")) %>% writeLines(., con=paste(input_path, "Log.txt", sep="\\"))


    output <- list(raw_data = raw_data, baseline_corrected = baseline_corrected, crosstalk_corrected=crosstalk_corrected, competition_corrected=competition_corrected, crosstalk_analysis=crosstalk_analysis, competition_analysis=competition_analysis, thresholds=thresholds)

    # calculate volume

    vol_by_well = data_list %>%
      bind_rows() %>%
      mutate(`Is invalid` = ifelse(`Is invalid` == "0", "valid", "invalid")) %>%
      group_by(Well, `Cycled volume`, Channel, `Is invalid`) %>%
      summarize(N_total = n()) %>%
      pivot_wider(names_from= `Is invalid`, values_from=N_total) %>%
      mutate(`Cycled volume` = as.numeric(`Cycled volume`)*(valid/(invalid + valid)))

    output$Volumes <- vol_by_well
    output$snr <- snr_vec
    output$date <- plate_date %>% unique() %>% na.omit()
    output$multiple_occupancy_outliers <- multiple_occupancy_outliers
    output$correlated_fluorescence_outliers <- z_score_neg_outliers
    output$aberrant_fluorescence_outliers <- z_score_pos_outliers
    plate_results[[plates[[plate]]]] <- output
  }

  return(plate_results)
}


#' @title Save output of QIAcuity Analyisis
#' @description  This function summarizes and visualizes the results of one or multiple QIAcuity experiments.
#' @param results A list of dataframes as returned by the function QIAcuityAnalysis
#' @param output_path path to where files and plots should be saved. Will create directory for the experiment there.
#' @param scatterplots_2d Set to TRUE to generate 2D scatterplots at every selected correction step. Default: TRUE
#' @param scatterplots_1d Set to TRUE to generate 1D scatterplots at every selected correction step. Default: TRUE
#' @param save_raw_data Whether plots and summaries based on raw data should be generated and saved. Default: TRUE
#' @param save_baseline_corr Whether plots and summaries based on baseline corrected data should be generated and saved. Default: TRUE
#' @param save_crosstalk_corr Whether plots and summaries based on crosstalk corrected data should be generated and saved. Default: TRUE
#' @param save_competition_corr Whether plots and summaries based on competition corrected data should be generated and saved. Default: TRUE
#' @param save_fluorescence_data Whether fluorescence data should be saved for all selected correction steps. Size: up to 70 MB per file. Default: TRUE
#' @returns NULL. Generates output files in the selected directory:
#' @returns    **Result summary:**  \[experiment name\]\_ result\_\[correction step\]\_corrected.xlsx <br />
#' Contains summary of dPCR results at this correction step. Columns: <br />
#' **Sample:** Sample ID <br />
#' **Well:** Well on plate <br />
#' **Channel:** Data from this channel is summarized in this line <br />
#' **valid:** Number of valid partitions in this well in this channel <br />
#' **negative:** Number of negative partitions in this well in this channel <br />
#' **positive:** Number of positive partitions in this well in this channel <br />
#' **CI:** 95% confidence interval of the quantification <br />
#' **poisson_corrected_targets:** Total number of copies of the DNA target in this channel in this well <br />
#' **mean_z_score:** Based on a positive sample. Z-score quantifies, by how many standard deviations the fluorescence intensity of a partition deviates from the mean fluorescence intensity of positive partitions in a positive control. Values lower than -2 usually indicate that a partition is background noise, while values greater than 2 can indicate imaging artifacts. <br />
#' **Warning:** Warnings are printed if there is any indication of artifacts or uncertain quantification <br />
#' @returns    **2d scatterplots:** \[experiment name\]\_ \[correction step\]\_\[channel x axis\]\[channel y axis\]\[channel color\].png
#' @returns    **1d scatterplots:** \[experiment name\]\_ \[correction step\]\_\[channel\].png
#' @returns    **thresholds:** \[experiment name]\_ thresholds.xlsx
#' @returns    **crosstalk analysis:** \[experiment name]\_ crosstalk\_analysis.xlsx
#' @returns    **competition analysis:** \[experiment name]\_ competition\_analysis.xlsx
#' @returns    **Signal-to-noise ratio analysis:** \[experiment name]\_snr_\[correction step\].xlsx
#' @returns    **Fluorescence data:** \[experiment name\]\_\[correction step\].csv <br />
#' Contains the fluorescence data. Each partition has a defined position on the plate \(overall partition ID, Well and ID of the partition within the well\), sample name, and fluorescence intensities in all the channels.

saveQIAcuityResult <- function(results, output_path, scatterplots_2d =TRUE, scatterplots_1d=TRUE,  save_raw_data = TRUE, save_baseline_corr = TRUE, save_crosstalk_corr = TRUE, save_competition_corr = TRUE, save_fluorescence_data  =TRUE){
  plate_summaries <- rep(NA, times=length(results)) %>% as.list()
  names(plate_summaries) <- names(results)
  corrections = c("raw_data", "baseline_corrected", "crosstalk_corrected", "competition_corrected")
  save_data = c(save_raw_data, save_baseline_corr, save_crosstalk_corr, save_competition_corr)


  for(plate in 1:length(results)){
    #dirname <- paste(output_path, "\\plate ", plate, "_", names(results)[[plate]], sep="")
    #dir.create(dirname)
    dirname <- output_path
    channels <- results[[plate]]$thresholds %>% pull(channel)
    n_channels <- results[[plate]]$thresholds %>% nrow()

    for(correction_step in 1:length(corrections)){
      if(save_data[[correction_step]] == FALSE){
        next()
      }
      assign("last.warning", NULL, envir = baseenv())
      well_order <- results[[plate]][[correction_step]] %>% select(Well) %>% unique() %>%
        mutate(letter=Well %>% gsub("[[:digit:]]", "", .)) %>% mutate(digit=Well %>%
            gsub("[[:alpha:]]", "", .)) %>% arrange(digit, letter) %>% pull(Well)
      if(save_fluorescence_data == TRUE){
        write.csv(results[[plate]][[correction_step]], file=paste(dirname, "\\", names(results)[[plate]], "_", corrections[[correction_step]], ".csv", sep=""))
      }
      channel_summaries <- data.frame(Sample=character(length=0), Well=character(length=0),
                                      negative=numeric(length = 0), positive = numeric(length = 0),
                                      Channel = character(length=0), poisson_corrected_targets = numeric(length=0))
      for(channel1 in 1:n_channels){
        if(scatterplots_2d==TRUE){
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
                plotname <- paste(names(results)[[plate]], corrections[[correction_step]], paste(results[[plate]]$thresholds$channel[c(channel2, channel1, p)], collapse=""), ".png", sep="" ) %>% gsub("[[:space:]]", "", .) %>% gsub("corrected", "corr", .) %>% paste(dirname, "\\", ., sep="")

                ggsave(scatterplot, filename = plotname, width=4, height=4)
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

          plotname <- paste(names(results)[[plate]], "_", corrections[[correction_step]], "_", results[[plate]]$thresholds$channel[[channel1]], ".png", sep="") %>% gsub("corrected", "corr", .) %>% paste(dirname, "\\", ., sep="")
          ggsave(plot, filename = plotname, width=14, height=7)
        }


        temp_res <- results
        if(correction_step>1){

          temp_res[[plate]][[corrections[[correction_step]]]]  <- temp_res[[plate]][[corrections[[correction_step]]]] %>% mutate(new = ifelse(temp_res[[plate]][[corrections[[correction_step]]]][temp_res[[plate]]$thresholds$channel[[channel1]]]>temp_res[[plate]]$thresholds[channel1, corrections[[correction_step]] %>% gsub("_corrected", "", .)], 1, 0)) %>% select(!temp_res[[plate]]$thresholds$channel[[channel1]]) %>% dplyr::rename(!!temp_res[[plate]]$thresholds$channel[[channel1]] := new)

          negs <- temp_res[[plate]][[corrections[[correction_step]]]] %>% filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])==0) %>%
            group_by(Well, Sample) %>% count()
          pos <-  temp_res[[plate]][[corrections[[correction_step]]]] %>% filter(!!sym(temp_res[[plate]]$thresholds$channel[[channel1]])==1) %>%
            group_by(Well, Sample) %>% count()

          summary <- merge(negs, pos, all=TRUE, by=c("Well", "Sample")) %>% dplyr::rename(positive=n.y, negative=n.x) %>%
            mutate(Channel=temp_res[[plate]]$thresholds$channel[[channel1]])

          summary <- summary %>% mutate(positive= replace(positive, is.na(positive), 0))

          #summary <- full_join(summary, zscore_RFU_pos) %>% dplyr::rename("mean_z_score"=z_score)
          channel_summaries <- bind_rows(channel_summaries, summary)


        }
      }
      if(!correction_step==1){

        channel_summaries <- channel_summaries %>% mutate(positive=replace(positive, is.na(positive), 0)) %>%
                                mutate(negative=replace(negative, is.na(negative), 0)) %>%
                                  mutate(valid=positive + negative) %>%
          mutate(poisson_corrected_targets = (-1* (positive + negative)*log(negative/(positive+negative))) %>% round(digits=2))


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

        channel_summaries <- channel_summaries %>% mutate(CI=lapply(.$positive, poisson.CI), .after=poisson_corrected_targets)

        channel_summaries <- channel_summaries %>% mutate(poisson_corrected_targets = (-1* (positive/factor + negative/factor)*log(negative/factor/(positive/factor+negative/factor))) %>% round(digits=2)) %>% mutate()

        if(length(results[[plate]]$date) == 0){
          results[[plate]]$date = NA
        }


        channel_summaries <- channel_summaries %>%
          mutate(min = CI %>% as.character() %>% unlist() %>% strsplit(" - ") %>% lapply(., as.numeric) %>% lapply(., min) %>% unlist()) %>%
          mutate(max = CI %>% as.character() %>% unlist() %>% strsplit(" - ") %>% lapply(., as.numeric) %>% lapply(., max) %>% unlist()) %>%
          mutate(min = (-1* (min/factor + (valid-min)/factor)*log((valid-min)/factor/(min/factor+(valid-min)/factor))) %>% round(digits=2)) %>%
          mutate(max = (-1* (max/factor + (valid-max)/factor)*log((valid-max)/factor/(max/factor+(valid-max)/factor))) %>% round(digits=2)) %>%
          mutate(max = replace(max %>% as.character(), max %>% as.character() =="NaN", "Inf")) %>%
          mutate(CI = paste(min, "-", max)) %>%
          select(!any_of(c("min", "max"))) %>%
          relocate(any_of(c("Sample", "Well", "Channel", "valid",  "negative", "positive", "poisson_corrected_targets", "CI"))) %>%
          dplyr::rename(`CI 95%` = CI) %>%
          mutate(Date = results[[plate]]$date) %>%
          mutate(Experiment = names(results)[[plate]])


        channel_summaries <- channel_summaries %>%
          arrange(factor(Well, levels = well_order)) %>%
          select(!factor)


        channel_summaries <- channel_summaries %>% mutate(Warning = ifelse(lapply(.$CI %>% as.character() %>% strsplit(split=" - "), min)<1 & .$positive > 0, "Warning: confidence interval contains values < 1, positive result may not be reliable.", NA) %>% paste(channel_summaries$Warning, sep=", ") %>% gsub("NA, NA", NA, .) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .))

        # add info about outliers
        if(correction_step == 4){
          # first, add info about partitions with aberrant z-score only
          if(length(channels) > 1){
          multiple_occupancy_outliers <- results[[plate]][["multiple_occupancy_outliers"]]

          multiple_occupancy_outliers_summarized <- multiple_occupancy_outliers %>%
            pivot_longer(., values_to = "status", names_to = "Channel", cols = channels) %>%
            filter(status==1) %>%
            group_by(Sample, Well, Channel) %>%
            summarize(`Multiple occupancy outliers` = n(), `Maximum expected value of outliers` = max(expected_value))
          }else{
            multiple_occupancy_outliers <- data.frame(Partition = numeric(), Well = character(), Sample = character())
            multiple_occupancy_outliers_summarized <- channel_summaries %>% select(Sample, Well, Channel) %>% mutate(`Multiple occupancy outliers` = 0)
          }
          channel_summaries <- full_join(channel_summaries, multiple_occupancy_outliers_summarized, by = c("Sample", "Well", "Channel")) %>%
            mutate(`Multiple occupancy outliers` = replace(`Multiple occupancy outliers`, is.na(`Multiple occupancy outliers`), 0)) %>%
            #relocate(., `Multiple occupancy outliers`, `Maximum expected value of outliers`, .after = mean_z_score) %>%
            mutate(Warning_outliers = ifelse(`Multiple occupancy outliers` > 0, "Warning: partitions with unusual multiple occupation, potential artifact.", NA)) %>%
            mutate(Warning = paste(Warning, Warning_outliers) %>% gsub("NA$|^NA", "", .)) %>%
            select(!Warning_outliers)

          correlated_fluorescence_outliers <- results[[plate]]$correlated_fluorescence_outliers %>%
            anti_join(., multiple_occupancy_outliers, by = c("Partition", "Well", "Sample")) %>%
            group_by(Well, Sample, channel) %>%
            summarize(`Correlated fluorescence outliers` = n(), median_z_score_mean = mean(median_z_score))

          channel_summaries <- left_join(channel_summaries, correlated_fluorescence_outliers %>% dplyr::rename(Channel = "channel"), by = c("Well", "Sample", "Channel")) %>%
            mutate(`Correlated fluorescence outliers` = replace(`Correlated fluorescence outliers`, is.na(`Correlated fluorescence outliers`), 0)) %>%
            mutate(Warning_outliers = ifelse(`Correlated fluorescence outliers` > 0, "Warning: partitions with increased fluorescence in all channels, potential artifact.", NA)) %>%
            mutate(Warning = paste(Warning, Warning_outliers) %>% gsub("NA$|^NA", "", .)) %>%
            select(!Warning_outliers)

          z_score_outliers <- results[[plate]]$aberrant_fluorescence_outliers %>%
            anti_join(., multiple_occupancy_outliers, by = c("Partition", "Well", "Sample")) %>%
            anti_join(., results[[plate]]$correlated_fluorescence_outliers, by = c("Partition", "Well", "Sample")) %>%
            group_by(Well, Sample, channel) %>%
            summarize(`Aberrant fluorescence outliers` = n(), z_score_mean_of_outliers = mean(z_score))

          channel_summaries <- left_join(channel_summaries, z_score_outliers %>% dplyr::rename(Channel = "channel"), by = c("Well", "Sample", "Channel")) %>%
            mutate(`Aberrant fluorescence outliers` = replace(`Aberrant fluorescence outliers`, is.na(`Aberrant fluorescence outliers`), 0)) %>%
            mutate(expected_outliers = positive*(2*0.021)) %>%
            mutate(`Aberrant fluorescence outliers` = `Aberrant fluorescence outliers`-round(expected_outliers, digits = 0)) %>%
            mutate(`Aberrant fluorescence outliers` = replace(`Aberrant fluorescence outliers`, `Aberrant fluorescence outliers` < 0, 0)) %>%
            mutate(`Aberrant fluorescence outliers` = replace(`Aberrant fluorescence outliers`, z_score_mean_of_outliers > -2, 0)) %>%
            mutate(Warning_outliers = ifelse(`Aberrant fluorescence outliers` > expected_outliers, "Warning: more partitions with extreme fluorescence values in this channel than expected, potential artifact.", NA)) %>%
            mutate(Warning = paste(Warning, Warning_outliers) %>% gsub("NA$|^NA", "", .)) %>%
            select(!Warning_outliers) %>%
            select(!expected_outliers)
        }



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

        saveWorkbook(excel_output, overwrite=FALSE, file = paste(dirname, "\\", names(results)[[plate]], "result_", corrections[[correction_step]], ".xlsx", sep=""))
      }
    }
    write.xlsx(results[[plate]]$crosstalk_analysis, file= paste(dirname, "\\", names(results)[[plate]], "_crosstalk_analysis.xlsx", sep = ""))
    write.xlsx(results[[plate]]$competition_analysis, file= paste(dirname, "\\", names(results)[[plate]], "_competition_analysis.xlsx", sep = ""))
    write.xlsx(results[[plate]]$thresholds, file= paste(dirname, "\\", names(results)[[plate]], "_thresholds.xlsx", sep = ""))
    write.xlsx(results[[plate]]$snr %>% data.frame(), file= paste(dirname, "\\", names(results)[[plate]], "_snr_competition_corrected.xlsx", sep = ""))
  }
}
