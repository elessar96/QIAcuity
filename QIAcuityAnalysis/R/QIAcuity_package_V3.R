setup_conducted <- FALSE

#' @title Read QIAcuity data
#' @description  This function reads in raw data .csv-files from the QIAcuity system
#' @param filename Path to the location of the QIAcuity .csv-file to be read.
#' @return An object of class dataframe of the RFU values from the .csv-file

read.data.QIAcuity <- function(filename){
  out <- read.csv(filename, header=FALSE)
  colnames(out) <- out[2,]
  out <- out[c(-1,-2),]

  return(out)
}

#' @title Get quantile
#' @description This function gives just the numeric value of a defined quantile in a vector of data.
#' @param input Input vector of data that can be forced to numeric.
#' @param probs quantile to be returned
#' @return A numeric value.
#'
get_quantile <- function(input, probs=0.25){
  quants <- input %>% as.numeric() %>% na.omit() %>% quantile(., probs=probs)
  return(quants[[1]])
}

check_NAs <- function(vec){
  vec_tab <- vec %>% unlist() %>% is.na() %>% table()
  if(!"TRUE" %in% names(vec_tab)){
    vec_tab <- c(vec_tab, "TRUE"=0)
  }

  if(vec_tab["TRUE"]<sum(vec_tab)*0.99){
    return(TRUE)
  }else{
    return(FALSE)
  }

}

#' @title Find turnpoints in a distribution
#' @description This function detects turnpoints in a distribution
#' @param data A dataframe or matrix only with numeric RFU values. Columns: Channels, Rows: partitions. Important:
#' @param variable A character specifying the column of data to be examined
#' @return A dataframe with turnpoints

find_turnpoints <- function(data, variable, n_points=30){
  # filter out non-numeric data
  data <- data %>% mutate_if(is.factor, as.numeric) %>% mutate_if(is.character, as.numeric) %>% select_if(check_NAs) %>% na.omit()

  # fit density curve to data
  d <- density(x=data %>% pull(variable), bw="sj", n=n_points)

  # find turnpoints
  ts_y <- ts(d$y)

  tp <- turnpoints(ts_y)
  tp_df <- data.frame(d.y=tp$points, d.x= d$x[tp$pos], tp$peaks, tp$pits) %>% right_join(., data.frame(d$x, d$y))

  # if first turn point is not a peak, prevent loss of a potential peak, and add pits at beginning and end if necessary:

  if(!tp$firstispeak){
    omitted_points <- tp_df %>% filter(d.x < (tp_df %>% filter(tp.pits==TRUE) %>% filter(d.x==min(d.x)) %>% pull(d.x))) %>% filter(d.y>0)
    if(nrow(omitted_points)>0){
      tp_df$tp.peaks[which(tp_df$d.x == omitted_points %>% filter(d.y==max(d.y)) %>% pull(d.x))] <- TRUE
      tp_df <- data.frame(d.y=0, d.x = (tp_df %>% filter(tp.peaks==TRUE) %>% pull(d.x) %>% min())-10, tp.peaks=FALSE, tp.pits=TRUE) %>% rbind(., tp_df)
    }
  }else{
    tp_df <- data.frame(d.y=0, d.x = (tp_df %>% filter(tp.peaks==TRUE) %>% pull(d.x) %>% min())-10, tp.peaks=FALSE, tp.pits=TRUE) %>% rbind(., tp_df)
  }

  if(tp_df %>% filter(d.x==max(d.x)) %>% pull(tp.peaks) == FALSE){
    tp_df <- data.frame(d.y=0, d.x = (tp_df %>% filter(tp.peaks==TRUE) %>% pull(d.x) %>% max())+10, tp.peaks=FALSE, tp.pits=TRUE) %>% rbind(., tp_df)
  }

  tp_df <- tp_df %>% filter(tp.peaks==TRUE | tp.pits == TRUE) %>% arrange(d.x)

  #get rid of multiple following points counted as same pit/peak!
  for(row in 1:(nrow(tp_df)-1)){
    if(tp_df$tp.peaks[[row]] == tp_df$tp.peaks[[row+1]] ){
      tp_df$d.x[c(row, row+1)] <- mean(tp_df$d.x[c(row, row+1)])
      tp_df$d.y[c(row, row+1)] <- mean(tp_df$d.y[c(row, row+1)])
    }
  }

  tp_df <- tp_df %>% unique()

  return(tp_df)
}

#' @title Classify peaks previously detected by find_turnpoints()
#' @description This function classifies turnpoints in a distribution
#' @param turnpoints Output of find_turnpoints
#' @param intensities Fluorescence intensity values in a data frame
#' @param variable Which column in intensities is the channel to be examined
#' @param channel_maxima A named list giving the maximum values expected for each of the channels found in intensities
#' @param reference_peaks A dataframe of previously classified peaks to guide this classification; optional.
#' @return A dataframe with turnpoints

classify_peaks <- function(turnpoints, intensities, variable, channel_maxima, reference_peaks=NULL){
  # find peaks and pits
  peaks <- turnpoints %>% filter(tp.peaks) %>% arrange(d.x)
  pits <- turnpoints %>% filter(tp.pits)

  peak_info <- mutate(intensities, peak=cut(intensities %>% pull(variable) %>% as.numeric(),breaks=pits$d.x, labels=paste("peak", c(1:nrow(peaks))))) %>% na.omit()

  #normalize intensity data
  other_channels <- colnames(intensities)
  for(chan in 1: length(other_channels)){
    peak_info[,other_channels[[chan]]] <- peak_info[,other_channels[[chan]]]/channel_maxima[other_channels[[chan]]]
  }
  peak_info <- peak_info %>% na.omit() %>% group_by(peak) %>% summarize_at(.vars = colnames(intensities), .funs= c("min"=min, "median"=median, "max"=max)) %>% na.omit() %>% mutate_if(is.numeric, round, digits=1)

  if(ncol(peak_info)==4){
    peak_info <- peak_info %>% rename(., !!sym(paste(variable, "_median", sep="")):=median, !!sym(paste(variable, "_min", sep="")):=min, !!sym(paste(variable, "_max", sep="")):=max)# %>% rename(!!sym(paste(variable, "median")) := median) %>% rename(!!sym(paste(variable, "max")) := max)
  }

  # add position and height info
  peak_indices <- peak_info %>% pull(peak) %>% gsub("peak ", "", .) %>% as.numeric()
  peak_info <- peak_info %>% mutate(position=peaks$d.x[peak_indices]) %>% mutate(height=peaks$d.y[peak_indices])

  #summarize crosstalk from other channels
  peak_info <- peak_info %>% mutate(max_others = peak_info %>% select(!starts_with(variable)) %>%  select(ends_with("median")) %>% apply(., 1, max))

  # rank intensity in current channel
  peak_info <- peak_info %>% mutate(median_current = peak_info %>% select(starts_with(variable)) %>% select(ends_with("median")) %>% unlist())

  if(!is.null(reference_peaks)){
    # if reference TN and TP peaks are provided, make decision purely on peak positions!
    if(nrow(peak_info)>2){
      ref_height <- peak_info %>% mutate(rank_height = rank(-height)) %>% filter(rank_height ==2 ) %>% pull(height)
      peak_info <- peak_info %>% filter(height>(ref_height/5))
    }

    peak_info <- peak_info %>% mutate(dist.tn=abs(position - reference_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x))) %>% mutate(dist.tp=abs(position - reference_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x)))
    peak_info <- peak_info %>% mutate(tn.peak.dist = ifelse(dist.tn == dist.tn %>% abs() %>% min(), 1, 0)) %>% mutate(tp.peak = ifelse(dist.tp == dist.tp %>% abs() %>% min(), TRUE, FALSE))


    tn_candidates <- peak_info %>% mutate(rank_height = rank(height)) %>% mutate(rank=tn.peak.dist*2 + rank_height) %>% filter(rank == max(rank))
    tp_candidates <- peak_info %>% filter(tp.peak==TRUE) %>% filter(dist.tn > abs(reference_peaks[1,]$d.x - reference_peaks[2,]$d.x)/4) %>% arrange(peak) %>% mutate(rank_dist.tn =rank(-dist.tn)) %>% mutate(rank_dist.tp =rank(dist.tp)) %>% mutate(max_others = rank(max_others)) %>% mutate(rank = max_others + rank_dist.tp) %>% filter(rank == min(rank))

    #if the same peak is classified as both TN and TP, decide based on distance to reference peaks:
    if(nrow(tn_candidates)==1 & nrow(tp_candidates)==1){
      if(tn_candidates$peak[[1]] == tp_candidates$peak[[1]]){
        if(tn_candidates$dist.tn > tn_candidates$dist.tp){
          tn_candidates <- tn_candidates[-1,]
        }else{
          tp_candidates <- tp_candidates[-1,]
        }
      }

    }

  }else{
    #find candidates based on height, position and crosstalk intensity
    if(nrow(peak_info)>2){
      tn_candidates <- peak_info %>% filter(height > max(height)/100)
      if(length(other_channels)==0){
        tn_candidates <- tn_candidates %>% arrange(peak) %>% mutate(position = rank(abs(as.numeric(position)))) %>% mutate(rank = position) %>% filter(rank == min(rank))
      }else{
        ct_detection <- tn_candidates %>% mutate(across(.cols=matches("median"), .fns=function(x){return(x-min(x))}))
        ct_detection <- ct_detection %>% mutate(max_others = ct_detection %>% select(!starts_with(variable)) %>% select(matches("median")) %>% apply(., MARGIN=1, FUN=max))

        non_ct_peaks <- ct_detection %>% filter(max_others < min(max_others)+0.5*sd(max_others))
        if(nrow(non_ct_peaks)==0){
          non_ct_peaks <- tn_candidates %>% filter(max_others <= min(max_others) + 0.2)
        }
        tn_candidates <- non_ct_peaks %>% arrange(peak)  %>% mutate(max_others = rank(max_others)) %>% mutate(position = rank(abs(as.numeric(position)))) %>% mutate(rank = max_others + position) %>% filter(rank == min(rank)) %>% filter(position==min(position))
      }
      tp_candidates <- peak_info %>% filter(position > (peak_info %>% filter(peak %in% tn_candidates$peak) %>% pull(position) %>% max())) %>% filter(!peak %in% c(tn_candidates %>% filter(position==min(position)) %>% pull(peak))) %>% arrange(peak) %>% mutate(pos = peaks$d.x[.$peak]) %>% mutate(max_others = rank(max_others)) %>% mutate(position = rank(-position)) %>% mutate(rank = max_others + position) %>% filter(rank == min(rank))
    }else{
      tn_candidates <- peak_info %>% filter(position==min(position))
      tp_candidates <- peak_info %>% filter(position==max(position))
    }
  }

  tn_peak <- peaks[tn_candidates %>% filter(max_others==min(max_others)) %>% filter(position==min(position)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]

  tp_peak <- peaks[tp_candidates %>% filter(position==min(position)) %>% filter(max_others==min(max_others)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]

  if(nrow(tn_peak)==1 & nrow(tp_peak)==1){
    if(tn_peak$d.x==tp_peak$d.x){
      tn_peak <- peaks[tn_candidates %>% filter(position==min(position)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
      tp_peak <- peaks[tp_candidates %>% filter(position==max(position)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
    }
  }
  #define true positive and true negative peaks

  other_peaks <- peak_info %>% filter(!peak %in% (tp_candidates %>% pull(peak))) %>% filter(!peak %in% (tn_candidates %>% pull(peak)))

  peaks <- mutate(peaks, crosstalk.peak = FALSE)

  #if other peaks are found - figure out whether these are due to crosstalk!
  if(nrow(other_peaks)>0){
    crosstalk_peaks <- list()

    for(opeak in 1:nrow(other_peaks)){

      channel_candidates <- other_peaks[opeak,] %>% select(ends_with("median"))
      channel_crosstalk <- colnames(channel_candidates)[which(channel_candidates == (other_peaks[opeak,] %>% pull(max_others)))] %>% gsub("_median", "", .)

      if((other_peaks[opeak,] %>% pull(max_others)) > (other_peaks[opeak,] %>% pull(!!sym(paste(variable, "_median", sep=""))))){
        peaks$crosstalk.peak[other_peaks[opeak,] %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric()] <- TRUE
      }#else{
        #print(paste("Crosstalk from", variable, "to", channel_crosstalk[[1]]))
      #}
    }
  }
  if(nrow(tp_peak)>0){
    peaks <- peaks %>% mutate(tp.peak = ifelse(d.x==tp_peak$d.x, TRUE, FALSE))
  }else{
    if(!tn_peak$d.y == peaks %>% filter(crosstalk.peak ==FALSE)  %>% pull(d.y) %>% max()){
      tn_peak <- peaks %>% filter(crosstalk.peak ==FALSE) %>% filter(d.y==max(d.y))
    }
    peaks <- peaks %>% mutate(tp.peak = FALSE)
  }
  if(nrow(tn_peak)>0){
    peaks <- peaks %>% mutate(tn.peak = ifelse(d.x==tn_peak$d.x, TRUE, FALSE))
  }else{
    peaks <- peaks %>% mutate(tn.peak = FALSE)
  }


  tp_df_mod <- full_join(peaks, pits %>% mutate(tp.peak=FALSE) %>% mutate(tn.peak=FALSE) %>% mutate(crosstalk.peak=FALSE)) %>% arrange(d.x)

  return(tp_df_mod)
  }



#' @title Define a threshold based on mixture modeling
#' @description This function defines a threshold that best divides a population of data points into two separate groups based on mixture modeling.
#' @param input_data An object of class dataframe containing fluorescence intensities fom ddPCR partitions.
#' @param variable  Name of the variable in the dataframe to be pulled
#' @param min_dist Determines the minimum distance between threshold and negative partitions, always as a fraction of the distance between true positive and true negative populations.
#' @param references Reference peaks for guiding classification
#' @param pc_data Positive control data as a reference to calculate maximum expected values for each channel
#' @return A numeric value giving the optimal threshold for distinguishing the true positive and true negative populations in the data.

density_threshold <- function(input_data, variable, min_dist=0.2, references=NULL, pc_data =NULL){

  data_subset <- input_data %>% select(!any_of(c("Well", "Sample", "Partition"))) %>% na.omit()

  # find turnpoints
  tp_df <- find_turnpoints(data_subset, variable=variable)

  if(is.null(pc_data)){
    pc_data <- data_subset
  }

  pc_data <- pc_data %>% select(!any_of(c("Well", "Sample", "Partition")))

  channel_maxima <- pc_data %>% apply(., 2, max, na.rm=TRUE)

  # find true positives and negatives
  classified_turnpoints <- classify_peaks(tp_df, intensities = data_subset , channel_maxima = channel_maxima, variable=variable, reference_peaks = references)

  #determine position of threshold
  positions <- data.frame(tp=NA, tn=NA, ct=NA)
  positions$tp[[1]] <- classified_turnpoints %>% filter(tp.peak == TRUE) %>% pull(d.x) %>% max()
  positions$tn[[1]] <- classified_turnpoints %>% filter(tn.peak == TRUE) %>% pull(d.x) %>% max()

  pot_crosstalk_peaks <- classified_turnpoints
  if(is.numeric(positions$tp)){
    pot_crosstalk_peaks <- pot_crosstalk_peaks %>% filter(d.x < positions$tp)
  }
  if(is.numeric(positions$tn)){

    pot_crosstalk_peaks <- pot_crosstalk_peaks %>% filter(d.x > positions$tn)

  }
  positions$ct[[1]] <-   pot_crosstalk_peaks %>% filter(crosstalk.peak == TRUE) %>% pull(d.x) %>% max()

  min_thresh <- positions$tn[[1]]

  if(is.finite(min_thresh)){
    if(positions$ct > min_thresh & is.finite(positions$ct)){
      if(is.finite(positions$tp)){
        if(positions$tp < positions$ct){
          min_thresh <- positions$tn
        }else{
          min_thresh <- positions$ct
        }

      }else{
        min_thresh <- positions$ct
      }

    }
  }else{
    min_thresh <- classified_turnpoints %>% filter(tp.pits==TRUE) %>% min()
  }

  threshold <- classified_turnpoints %>% filter(tp.pits == TRUE) %>% filter(d.x > min_thresh) %>% filter(d.x == min(d.x)) %>% pull(d.x)

  if(is.finite(positions$tp) & is.finite(positions$tn)){
    dist <- (threshold-positions$tn)/(positions$tp-positions$tn)
    if(dist < min_dist){
      threshold <- positions$tn + (positions$tp-positions$tn)*min_dist
    }
  }

  return(threshold)

}

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
  for(i in 1:length(plates)){
    files <- dir()[grep(plates[[i]], dir())]
    files <- files[grep(".zip", files)]
    if(length(files) > 1){
      rem <- c(rem, i)
      for(j in 1:length(files)){
        file.rename(files[[j]], paste("plate", j, files[[j]], sep="_"))
        plates <- c(plates, paste("plate", j, plates[[i]], sep="_"))
      }
    }
  }
  if(length(rem)>0){
    plates <- plates[-rem]
  }

  plate_results <- rep(NA, times=length(plates)) %>% as.list()
  names(plate_results) <- plates
  # for each experiment, do:
  for(k in 1:length(plates)){
    print(plates[k])
    setwd(input_path)
    # find all files with experiment name
    plate_files <- dir()[grep(plates[[k]], dir())]
    #if csv files among them, analyze those files (assumption: .zip has been unzipped!)
    csv_files <- plate_files[grep(".csv", plate_files)]

    # if no csv files, unzip any zip-files
    if(length(csv_files)==0){
      zip_file <- plate_files[grep(".zip", plate_files)]
      if(length(zip_file)==0){
        stop("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.")
      }
      dirname <- paste(output_path, "\\plate ", k, "_", plates[[k]], sep="")
      dir.create(dirname)
      unzip(zip_file, exdir=dirname, overwrite = FALSE)

      setwd(dirname)
      csv_files <- dir()[grep(paste(plates[[k]] %>% gsub("plate_\\d_", "", .), ".+", ".csv",sep=""), dir())]
      dataframes1 <- csv_files %>% lapply(., read.data.QIAcuity)

    }else{
      dataframes1 <- csv_files %>% lapply(., read.data.QIAcuity)
    }

    # gather data and find baseline of all samples, calculate corrected RFU values
    raw_data <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))

    for(channel in 1:length(dataframes1)){
      #extract RFU data from raw csv files
      temp <- dataframes1[[channel]]
      temp <- temp %>% data.frame() %>% mutate(RFU=RFU %>% as.numeric()) %>% group_by(Well) %>% mutate(Partition=Partition %>% as.numeric())

      #put those RFU values into one dataframe (that contains the data of all channels)
      raw_data <- temp %>% select(Well, Sample, Partition, RFU) %>% merge(., raw_data, all=TRUE, by=c("Well", "Sample", "Partition"))
      colnames(raw_data)[which(colnames(raw_data)=="RFU")] <- temp %>% pull(Channel) %>% unique()
    }

    channels <- raw_data %>% select(!Well) %>% select(!Partition) %>% select(!Sample) %>% colnames()

    coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0))

    wells <- raw_data %>% pull(Well) %>% unique()

    thresholds <- data.frame(channel=channels, baseline=numeric(length=length(channels)), crosstalk = numeric(length=length(channels)), competition=numeric(length=length(channels)))

    maxima <- raw_data %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), get_quantile, probs=0.999)

    pc_wells <- maxima %>% mutate(across(.cols=channels, .fns=function(x){return(x-min(x))})) %>% mutate(across(.cols=channels, .fns=function(x){return(x/max(x))}))  %>% mutate(minimum = apply(across(channels), 1, min)) %>% mutate(median = apply(across(channels), 1, median)) %>% filter(median == max(median) | minimum == max(minimum)) %>% pull(Well) %>% unique()

    #pc_wells <- maxima %>% mutate(across(.cols=channels, .fns=rank)) %>% mutate(rank= select(., any_of(channels)) %>% apply(., 1, max)) %>% filter(rank==max(rank)) %>% pull(Well) %>% c(pc_wells, .) %>% unique()

    baseline_corrected <- raw_data

    for(channel in 1:length(channels)){
      # find true positive and negative peaks in PC
      current_channel <- channels[[channel]]

      tp_data <- raw_data %>% filter(Well %in% pc_wells) %>% select(!any_of(c("Well", "Partition", "Sample"))) %>% na.omit()

      if(current_channel %in% coupled_channels){
        exclude <- coupled_channels %>% filter(ch1 == current_channel | ch2 ==current_channel) %>% unlist() %>% unique()
        exclude <- exclude[which(!exclude==current_channel)]

        tp_data <- tp_data %>% select(!any_of(exclude))
      }else{
        exclude <- c()
      }

      turnpoints <- tp_data  %>% find_turnpoints(., variable=current_channel)
      classified_turnpoints <- classify_peaks(turnpoints = turnpoints, intensities = tp_data, variable = current_channel, channel_maxima = tp_data %>% apply(., 2, get_quantile, probs=0.999))

      ref_peaks <- classified_turnpoints %>% filter(tp.peak == TRUE | tn.peak == TRUE)
      if((ref_peaks %>% filter(tp.peak==TRUE) %>%nrow())==0){
        tp_artif <- classified_turnpoints %>% filter(tp.peaks==TRUE) %>% filter(d.x > ref_peaks %>% filter(tn.peak=TRUE) %>% pull(d.x)) %>% filter(d.x==min(d.x)) %>% mutate(tp.peak=TRUE)
        ref_peaks <- ref_peaks %>% rbind(., tp_artif)
      }

      # calculate baseline for each well
      baseline <- data.frame(Well=character(length=0), Partition=character(length=0), baseline= numeric(length=0))

      for(well in 1:length(wells)){
        # try to find threshold for each well guided by positive controls, limited by crosstalk signals
        well_data <- raw_data %>% na.omit() %>% select(!any_of(exclude)) %>% filter(Well==wells[[well]])

        turnpoints_well <- well_data %>% select(!any_of(c("Well", "Sample", "Partition"))) %>% na.omit() %>% find_turnpoints(., variable=current_channel, n_points=15)
        classified_turnpoints_well <- classify_peaks(turnpoints = turnpoints_well, intensities = well_data %>% select(!any_of(c("Well", "Sample", "Partition"))), variable = current_channel, channel_maxima = tp_data %>% apply(., 2, get_quantile, probs=0.999), reference_peaks = ref_peaks)

        tn_peak_well <- classified_turnpoints_well %>% filter(tn.peak == TRUE)

        if(nrow(tn_peak_well)>0){
          upper_lim <- classified_turnpoints_well %>% filter(d.x > tn_peak_well %>% pull(d.x)) %>% pull(d.x) %>% min()
          lower_lim <- classified_turnpoints_well %>% filter(d.x < tn_peak_well %>% pull(d.x)) %>% pull(d.x) %>% max()
          datapoints <- well_data %>% filter(!!sym(current_channel) > lower_lim) %>% filter(!!sym(current_channel) < upper_lim) %>% nrow()
        }else{
          datapoints <- 0
        }

        if(nrow(tn_peak_well) == 0 | datapoints<(nrow(well_data)/100)){
          baseline <- data.frame(Well=well_data$Well, Partition = well_data$Partition, baseline = ref_peaks %>% filter(tn.peak == TRUE) %>% pull(d.x)) %>% rbind(baseline, .)
          warning(paste("Well", wells[[well]], ", Channel",  current_channel, ". Baseline estimation potentially failed. Please inspect data."))
        }else{
          if(smooth==TRUE){
              baseline <- smooth_data(well_data, upper_lim = upper_lim, lower_lim= lower_lim, current_channel=current_channel) %>% select(Well, Partition, baseline) %>% rbind(baseline, .)
            }else{
              baseline_well <- well_data %>% filter(!!sym(current_channel) > lower_lim) %>% filter(!!sym(current_channel) < upper_lim) %>% select(any_of(current_channel)) %>% unlist() %>% median()
              baseline <- data.frame(Well=well_data$Well, Partition = well_data$Partition, baseline = baseline_well) %>% rbind(baseline, .)
            }
          }
      }
      baseline_corrected <- baseline_corrected %>% select(any_of(c("Well", "Partition", "Sample", current_channel))) %>% merge(., baseline, by=c("Well", "Partition")) %>% mutate(!!sym(current_channel) := !!sym(current_channel) - baseline) %>% select(!baseline) %>% full_join(baseline_corrected %>% select(!any_of(current_channel)))
    }

    #check for correlation between channels to prevent later errors!
    maxima_baseline_c <- baseline_corrected %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), get_quantile, probs=0.999)
    tp_data <- baseline_corrected %>% filter(Well %in% pc_wells)
    maxima_channels <- tp_data %>% select(any_of(channels)) %>% apply(., 2, max, na.rm=TRUE)
    sd_channels <-  tp_data %>% select(any_of(channels)) %>% apply(., 2, sd, na.rm=TRUE)
    for(i in 1:length(channels)){

      tp <- tp_data %>% select(any_of(channels))  %>% find_turnpoints(., variable=channels[[i]])
      classified_turnpoints <- classify_peaks(tp, channel_maxima = maxima_channels , variable=channels[[i]], intensities =tp_data %>% select(any_of(channels)))
      true_pos <- classified_turnpoints %>% filter(tp.peak ==TRUE)
      lower_lim <- classified_turnpoints %>% filter(d.x<true_pos %>% pull(d.x)) %>% filter(d.x==max(d.x)) %>% pull(d.x)
      upper_lim <- classified_turnpoints %>% filter(d.x>true_pos %>% pull(d.x)) %>% filter(d.x==min(d.x)) %>% pull(d.x)

      peak_data <- tp_data %>% filter(!!sym(channels[[i]])>lower_lim) %>% filter(!!sym(channels[[i]])<upper_lim) %>% summarize_at(all_of(channels[-i]), median)
      peak_data <- peak_data + sd_channels[-i] -maxima_channels[-i]
      peak_data <- peak_data[which(peak_data > 0)]
      if(length(peak_data)>0){
        for(j in 1:length(peak_data)){
          pair <- data.frame(ch1=c(channels[i], names(peak_data)[j]), ch2=c(names(peak_data)[j], channels[i]))
          coupled_channels <- coupled_channels %>% rbind(., pair)
        }
      }
    }

    coupled_channels <- coupled_channels %>% unique()

    #define thresholds for each channel that separate positive and negative samples

    baseline_corrected_dichot <- baseline_corrected

    maxima <- baseline_corrected %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), get_quantile, probs=0.999)

    for(i in 1:length(channels)){
      #if(channels[[i]] %in% (coupled_channels %>% unlist())){
      #  omit <- coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]]) %>% unlist() %>% unique()
      #  comparison_channels <- channels[which(!channels %in% omit)]

      #}else{
      comparison_channels <- channels[-i]
      #}

      threshold <- baseline_corrected %>% filter(Well %in% pc_wells) %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.2)

      thresholds[which(thresholds$channel ==channels[[i]]), "baseline"]<- threshold

      # make dichotomized dataset that just gives info on whether a partition is positive or not
      baseline_corrected_dichot[,channels[[i]]] <- ifelse(baseline_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }

    # crosstalk correction using a linear mode
    crosstalk_corrected <- baseline_corrected
    crosstalk_analysis = data.frame(ch1=character(length=0), ch2=character(length=0), cross_talk=numeric(length=0))

    for(i in 1:length(channels)){
      for(j in 1:length(channels)){
        if(!i==j){
          #pc_wells_channel <- baseline_corrected %>% group_by(Well) %>% summarize(max_current=get_quantile(!!sym(channels[[i]]), probs=0.999)) %>% filter(max_current > max(max_current)-2*sd(max_current)) %>% pull(Well)
          pc_data <- crosstalk_corrected %>% filter(Well %in% pc_wells)
          excluded_channels <- coupled_channels %>% filter(ch1 %in% channels[c(i,j)] | ch2 %in% channels[c(i,j)]) %>% unlist() %>% unique()
          filtered_points <- baseline_corrected_dichot %>% filter(Well %in% pc_wells)
          filtered_points$sum <- filtered_points %>% select(any_of(channels)) %>% select(!any_of(channels[c(i,j)])) %>% select(!any_of(excluded_channels)) %>% apply(., 1, sum)
          filtered_points <- filtered_points %>% filter(sum==0)

          thr_ch2 <- thresholds %>% filter(channel == channels[[j]]) %>% pull(baseline)
          thr_ch1 <- thresholds %>% filter(channel == channels[[i]]) %>% pull(baseline)

          pc_data <- filtered_points %>% select(Well, Sample, Partition) %>% left_join(., pc_data)

          turnpoints_channel <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% select(any_of(channels[[i]])) %>% find_turnpoints(variable = channels[[i]])

          tn_peak <- turnpoints_channel %>% mutate(height=rank(-d.y)) %>% filter(height < 3) %>% filter(d.x==min(d.x))

          thr_ch1_ct <- turnpoints_channel %>% filter(d.x > tn_peak$d.x) %>% pull(d.x) %>% min()

          sp_points <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% filter(!!sym(channels[[i]]) < thr_ch1_ct)
          sp_points_ch2 <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2) %>% filter(!!sym(channels[[i]]) > thr_ch1_ct)

          dp_points <- pc_data %>% filter(!!sym(channels[[i]]) > thr_ch1_ct) %>% filter(!!sym(channels[[j]]) > thr_ch2)

          expected_dps <- nrow(pc_data) * ((nrow(sp_points) + nrow(dp_points))/nrow(pc_data)) * ((nrow(sp_points_ch2) + nrow(dp_points))/nrow(pc_data))

          dn_points <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2)  %>% filter(!!sym(channels[[i]]) < thr_ch1)

          if(nrow(dp_points)>expected_dps*100){
            coupled_channels <- data.frame(ch1=channels[[i]], ch2=channels[[j]]) %>% rbind(coupled_channels)
          }

          if(nrow(sp_points)<100){
            print("Warning: fewer than 100 positive points found. Aborting crosstalk calculation")
            next()
          }

          dn_med <- c(dn_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), dn_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          sp_med <- c(sp_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), sp_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))

          slope <- (sp_med[[1]] - dn_med[[1]])/(sp_med[[2]] - dn_med[[2]])
          yintercept <- dn_med[[1]] - slope*dn_med[[2]]

          crosstalk <- (crosstalk_corrected %>% pull(channels[[j]]))*slope + yintercept #predict(model_crosstalk, baseline_corrected)*(1-(corr_factor-1)^2)
          crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=slope) %>% bind_rows(., crosstalk_analysis)

          crosstalk_corrected[,channels[[i]]] <- crosstalk_corrected[,channels[[i]]] - crosstalk
        }
      }
    }

    crosstalk_corrected_dichot <- crosstalk_corrected

    maxima <- crosstalk_corrected %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), get_quantile, probs=0.999)

    for(i in 1:length(channels)){
      threshold <- crosstalk_corrected %>% filter(Well %in% pc_wells) %>% mutate(Sample=Well) %>% select(any_of(c(channels[[i]]))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.2)

      thresholds[which(thresholds$channel ==channels[[i]]), "crosstalk"]<- threshold

      # make dichotomized dataset that just gives info on whether a partition is positive or not
      crosstalk_corrected_dichot[,channels[[i]]] <- ifelse(crosstalk_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }

    # perform competition correction
    competition_corrected <- crosstalk_corrected
    competition_analysis = data.frame(ch1=character(length=0), ch2=character(length=0), competition=numeric(length=0))

    for(i in 1:length(channels)){
      for(j in 1:length(channels)){
        if(!i==j){
          pc_data <- competition_corrected %>% filter(Well %in% pc_wells)

          thr_ch1 <- thresholds %>% filter(channel==channels[[i]]) %>% pull(crosstalk)
          thr_ch2 <- thresholds %>% filter(channel==channels[[j]]) %>% pull(crosstalk)

          sp_points <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2) %>% filter(!!sym(channels[[i]]) > thr_ch1)
          sp_points_ch2 <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% filter(!!sym(channels[[i]]) < thr_ch1)

          dp_points <- pc_data %>% filter(!!sym(channels[[i]]) > thr_ch1) %>% filter(!!sym(channels[[j]]) > thr_ch2)

          expected_dps <- nrow(pc_data) * ((nrow(sp_points) + nrow(dp_points))/nrow(pc_data)) * ((nrow(sp_points_ch2) + nrow(dp_points))/nrow(pc_data))

          dn_points <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2)  %>% filter(!!sym(channels[[i]]) < thr_ch1)

          if(nrow(dp_points)<10){
            print("Warning: fewer than 10 double positive points found. Aborting competition calculation")
            next()
          }

          dp_med <- c(dp_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), dp_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          sp_med <- c(sp_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), sp_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          dn_med <- c(dn_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), dn_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))

          slope <- (sp_med[[1]] - dp_med[[1]])/(sp_med[[2]] - dp_med[[2]])

          max_ch1 <- max(dp_points %>% pull(channels[[i]]))

          competition <- ((competition_corrected %>% pull(channels[[j]]))*slope)* ((crosstalk_corrected %>% pull(channels[[i]])) /max_ch1) #predict(model_crosstalk, baseline_corrected)*(1-(corr_factor-1)^2)
          competition_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], competition=slope) %>% bind_rows(., competition_analysis)

          competition_corrected[,channels[[i]]] <- competition_corrected[,channels[[i]]] - competition

        }
      }
    }

    # re-calculate thresholds
    competition_corrected_dichot <- competition_corrected

    maxima <- competition_corrected %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), get_quantile, probs=0.999)

    for(i in 1:length(channels)){
      # define which channels the data should be compared to
      #if(channels[[i]] %in% (coupled_channels %>% unlist())){
      #  omit <- coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]]) %>% unlist() %>% unique()
      #  comparison_channels <- channels[which(!channels %in% omit)]
#
      #}else{
      #  comparison_channels <- channels[-i]
  #    }

      #calculate new threshold for this channel
     # threshold <- competition_corrected %>% filter(Well %in% pc_wells) %>% mutate(Sample=Well) %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.5)
      threshold <- competition_corrected %>% filter(Well %in% pc_wells) %>% mutate(Sample=Well) %>% select(any_of(c(channels[[i]]))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.4)

      thresholds[which(thresholds$channel ==channels[[i]]), "competition"]<- threshold

      # make dichotomized dataset that just gives info on whether a partition is positive or not
      competition_corrected_dichot[,channels[[i]]] <- ifelse(competition_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }

    snr_vec <- c()

    # calculate SNRs for each channel
    for(i in 1:length(channels)){
      pos <- competition_corrected %>% filter(!!sym(channels[[i]]) >= thresholds %>% filter(channel== channels[[i]]) %>% pull("competition")) %>% pull(channels[[i]])
      neg <- competition_corrected %>% filter(!!sym(channels[[i]]) < thresholds %>% filter(channel== channels[[i]]) %>% pull("competition")) %>% pull(channels[[i]])

      signal <- median(pos) - median(neg)
      noise <- mean(c(rep(sd(pos), times=length(pos)), rep(sd(neg), times=length(neg))))*2 + median(neg)

      SNR <- signal^2/noise^2
      snr_vec[[i]] <- SNR
      names(snr_vec)[[i]] <- channels[[i]]
    }

    output <- list(raw_data = raw_data, baseline_corrected = baseline_corrected, crosstalk_corrected=crosstalk_corrected, competition_corrected=competition_corrected, crosstalk_analysis=crosstalk_analysis, competition_analysis=competition_analysis, thresholds=thresholds)

    # calculate volume
    vol_by_well = dataframes1 %>% bind_rows() %>% group_by(Well, Channel) %>% select(Well, starts_with("Cycled"), Channel)%>% unique()
    output$Volumes <- vol_by_well
    output$snr <- snr_vec
    plate_results[[plates[[k]]]] <- output
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

  for(i in 1:length(results)){
    dirname <- paste(output_path, "\\plate ", i, "_", names(results)[[i]], sep="")
    dir.create(dirname)
    setwd(dirname)
    n_channels <- results[[i]]$thresholds %>% nrow()

    for(j in 1:length(corrections)){
      well_order <- results[[i]][[j]] %>% select(Well) %>% unique() %>% mutate(letter=Well %>% gsub("[[:digit:]]", "", .)) %>% mutate(digit=Well %>% gsub("[[:alpha:]]", "", .)) %>% arrange(digit, letter) %>% pull(Well)

      write.csv(results[[i]][[j]], file=paste(corrections[[j]], ".csv", sep=""))

      channel_summaries <- data.frame(Sample=character(length=0), Well=character(length=0), negative=numeric(length = 0), positive = numeric(length = 0), Channel = character(length=0), poisson_corrected_targets = numeric(length=0))
      for(k in 1:n_channels){
        if(detailed_plots==TRUE){
          for(m in 1:n_channels){
            if(!k>=m){
              if(n_channels>2){
                p = m-1
                if(p==k){
                  p = m-2
                }
                if(p==0){
                  p = m+1
                }
                scatterplot <- results[[i]][[corrections[[j]]]] %>% dplyr::select(Well | ends_with(results[[i]]$thresholds$channel[c(m, k, p)])) %>% na.omit() %>% ggplot(aes_string(x =results[[i]]$thresholds$channel[[m]], y =results[[i]]$thresholds$channel[[k]], color=results[[i]]$thresholds$channel[[p]])) + geom_point()
                if(!j==1){
                  scatterplot <- scatterplot + geom_hline(yintercept=results[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)], color="red") + geom_vline(xintercept=results[[i]]$thresholds[m, corrections[[j]] %>% gsub("_corrected", "", .)], color="red")

                }
                plotname <- paste(names(plate_summaries)[[i]], corrections[[j]], "2dscatter_", paste(results[[i]]$thresholds$channel[c(m, k, p)], collapse=""), ".png", sep="" )
                ggsave(scatterplot, filename = plotname)
              }
            }
          }
        }
        if(scatterplots_1d==TRUE){
          plot <- results[[i]][[corrections[[j]]]]  %>% dplyr::select(Well | Partition | ends_with(results[[i]]$thresholds$channel[k])) %>% na.omit() %>% mutate(Partition = Partition %>% as.numeric()) %>% mutate(across(Well, factor, levels=c(well_order))) %>% ggplot(aes_string(x="Partition", y =results[[i]]$thresholds$channel[[k]])) + geom_jitter(size=0.5) + scale_x_discrete(limits=well_order) + facet_wrap(vars(Well), nrow =2) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

          wells <- results[[i]][[corrections[[j]]]] %>% pull(Well) %>% unique()

          if(!j==1){
            plot <- plot + geom_hline(yintercept=results[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)], color="red")
          }

          plotname <- paste(names(plate_summaries)[[i]], "_", corrections[[j]], "_", results[[i]]$thresholds$channel[[k]], ".png", sep="")
          ggsave(plot, filename = plotname, width=14, height=7)
        }


        temp_res <- results
        if(j>1){
          temp_res[[i]][[corrections[[j]]]]  <- temp_res[[i]][[corrections[[j]]]] %>% mutate(new = ifelse(temp_res[[i]][[corrections[[j]]]][temp_res[[i]]$thresholds$channel[[k]]]>temp_res[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)], 1, 0)) %>% select(!temp_res[[i]]$thresholds$channel[[k]]) %>% rename(!!temp_res[[i]]$thresholds$channel[[k]] := new)

          negs <- temp_res[[i]][[corrections[[j]]]] %>% filter(!!sym(temp_res[[i]]$thresholds$channel[[k]])==0) %>% group_by(Well, Sample) %>% count()
          pos <-  temp_res[[i]][[corrections[[j]]]] %>% filter(!!sym(temp_res[[i]]$thresholds$channel[[k]])==1) %>% group_by(Well, Sample) %>% count()

          summary <- merge(negs, pos, all=TRUE, by=c("Well", "Sample")) %>% rename(positive=n.y, negative=n.x) %>% mutate(Channel=temp_res[[i]]$thresholds$channel[[k]])

          summary <- summary %>% mutate(positive= replace(positive, is.na(positive), 0)) %>% mutate(poisson_corrected_targets = (-1* (positive + negative)*log(negative/(positive+negative))) %>% round(digits=2))

          mean_RFU_pos <- results[[i]][[corrections[[j]]]] %>% filter(!!sym(temp_res[[i]]$thresholds$channel[[k]])>results[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)]) %>% pull(temp_res[[i]]$thresholds$channel[[k]]) %>% mean()
          sd_RFU_pos <- results[[i]][[corrections[[j]]]] %>% filter(!!sym(temp_res[[i]]$thresholds$channel[[k]])>results[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)]) %>% pull(temp_res[[i]]$thresholds$channel[[k]]) %>% sd()

          zscore_RFU_pos <- results[[i]][[corrections[[j]]]] %>% filter(!!sym(temp_res[[i]]$thresholds$channel[[k]])>results[[i]]$thresholds[k, corrections[[j]] %>% gsub("_corrected", "", .)]) %>% mutate(z_score = abs(!!sym(temp_res[[i]]$thresholds$channel[[k]])-mean_RFU_pos)/sd_RFU_pos) %>% group_by(Well) %>% summarize_at("z_score", mean)

          zscore_RFU_pos <- zscore_RFU_pos %>% mutate(Warning = ifelse(z_score>2, "Warning: mean RFU value of positive partitions further than 2 standard deviations away from global mean. Potential artifact.", NA))

          summary <- full_join(summary, zscore_RFU_pos) %>% rename("mean_z_score"=z_score)
          channel_summaries <- bind_rows(channel_summaries, summary)
        }
    }
    if(!j==1){

      channel_summaries <- channel_summaries %>% mutate(positive=replace(positive, is.na(positive), 0)) %>% mutate(negative=replace(negative, is.na(negative), 0)) %>% mutate(valid=positive + negative)

      mm_vol <- 0

      n_part <- channel_summaries %>% pull(valid)

      if(n_part %>% median(., na.rm = TRUE) > 10000){
        mm_vol <- 40/26000
      }else{
        mm_vol <- 12/8500
      }


      poisson.CI <- function(x){return(poisson.test(x)$conf.int %>% round(digits=1) %>% paste(., collapse=" - "))}

      channel_summaries <- merge(channel_summaries, results[[i]]$Volumes, by=c("Well", "Channel")) %>% mutate(factor = (!!sym("Cycled volume") %>% as.numeric())/(mm_vol*valid))


      channel_summaries <- channel_summaries %>% mutate(CI=lapply(.$positive, poisson.CI)) %>% relocate(any_of(c("Sample", "Well", "Channel", "negative", "positive", "CI")))

      channel_summaries <- channel_summaries %>% mutate(Warning = ifelse(lapply(.$CI %>% as.character() %>% strsplit(split=" - "), min)<1 & .$positive > 0, "Warning: confidence interval contains values < 1, positive result may not be reliable.", NA) %>% paste(channel_summaries$Warning, sep=", ") %>% gsub("NA, NA", NA, .) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .))


      channel_summaries <- channel_summaries %>% mutate(poisson_corrected_targets = (poisson_corrected_targets/factor) %>% round(digits=2)) %>% mutate()

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

      saveWorkbook(excel_output, overwrite=TRUE, file = paste("result_", corrections[[j]], ".xlsx", sep=""))
    }

    }
    write.xlsx(results[[i]]$crosstalk_analysis, file="crosstalk_analysis.xlsx")
    write.xlsx(results[[i]]$competition_analysis, file="competition_analysis.xlsx")
    write.xlsx(results[[i]]$thresholds, file="thresholds.xlsx")
    write.xlsx(results[[i]]$snr %>% data.frame(), file="snr_competition_corrected.xlsx")
    setwd("..")
  }
}

smooth_data <- function(data, current_channel, upper_lim, lower_lim){

  per_cluster <- list()

  breakpoint_data <- data %>% filter(!!sym(current_channel)>lower_lim) %>% filter(!!sym(current_channel)<upper_lim)

  frm <- formula(paste(current_channel, "~ s(Partition, bs='cs')"))

  model <- mgcv::gam(frm, data=breakpoint_data)

  breakpoint_data$predicted <- model$fitted.values
  breakpoint_data <- breakpoint_data %>% arrange(Partition)
  breakpoint_data$first_deriv <- (breakpoint_data$predicted - c(breakpoint_data$predicted[2:nrow(breakpoint_data)], NA))/(breakpoint_data$Partition - c(breakpoint_data$Partition[2:nrow(breakpoint_data)], NA))

  tps <- turnpoints(breakpoint_data$first_deriv)
  break_ID <- tps$tppos %>% na.omit()
  breakpoint_pos <- breakpoint_data$Partition[break_ID]

  breaks <- c(0, breakpoint_pos, data %>% pull(Partition) %>% max()) %>% unique()

  for(p in 2:length(breaks)){
    data_baseline_est <- breakpoint_data %>% filter(Partition > breaks[p-1]) %>% filter(Partition < breaks[p])

    if(nrow(data_baseline_est)>10){
      frm <- formula(paste(current_channel, "~ poly(Partition, 2)"))

      model_baseline <- data_baseline_est %>% lm(frm, .)

      pred <- predict(model_baseline, data %>% filter(Partition > breaks[p-1]) %>% filter(Partition <= breaks[p]))
      pred <- data.frame(Partition=names(pred) %>% as.numeric(), RFU_temp = pred)

      tempPerWell <- data %>% filter(Partition > breaks[p-1]) %>% filter(Partition <= breaks[p]) %>% mutate(!!sym(current_channel) := !!sym(current_channel) - pred$RFU_temp) %>% mutate(baseline=pred$RFU_temp)
      per_cluster[[p]] <- tempPerWell
    }else{
      pred <- breakpoint_data %>% select(all_of(current_channel)) %>% unlist() %>% median()
      tempPerWell <- data %>% filter(Partition > breaks[p-1]) %>% filter(Partition <= breaks[p]) %>% mutate(!!sym(current_channel) := !!sym(current_channel) - pred) %>% mutate(baseline=pred)
      per_cluster[[p]] <- tempPerWell
      }
  }
  output <- bind_rows(per_cluster)

  return(output)
  }
