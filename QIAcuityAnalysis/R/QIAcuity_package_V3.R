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
  
  if((vec_tab["TRUE"]/sum(vec_tab)) < 0.5){
    return(TRUE)
  }else{
    return(FALSE)
  }
  
}

#' @title Find turnpoints in a distribution
#' @description This function detects turnpoints in a distribution
#' @param data Test
#' @return A dataframe with turnpoints

find_turnpoints <- function(data, variable=NULL){
  # filter out non-numeric data
  data <- data %>% mutate_if(is.factor, as.numeric) %>% mutate_if(is.character, as.numeric) %>% select_if(check_NAs) %>% na.omit()
  
  # fit density curve to data
  d <- density(x=data %>% pull(variable), bw="sj", n=30)
  
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
#' @param data Test
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
  
  # summarize data
  peak_info <- peak_info %>% na.omit() %>% group_by(peak) %>% summarize_at(.vars = colnames(intensities), .funs= c("min"=min, "median"=median, "max"=max)) %>% na.omit() %>% mutate_if(is.numeric, round, digits=1)
  
  # add position and height info
  peak_indices <- peak_info %>% pull(peak) %>% gsub("peak ", "", .) %>% as.numeric()
  peak_info <- peak_info %>% mutate(position=peaks$d.x[peak_indices]) %>% mutate(height=peaks$d.y[peak_indices])
  
  #summarize crosstalk from other channels
  peak_info <- peak_info %>% mutate(max_others = peak_info %>% select(!starts_with(variable)) %>%  select(ends_with("median")) %>% apply(., 1, max))
  
  # rank intensity in current channel
  peak_info <- peak_info %>% mutate(median_current = peak_info %>% select(starts_with(variable)) %>% select(ends_with("median")) %>% unlist())
  
  if(!is.null(reference_peaks)){
    # if reference TN and TP peaks are provided, make decision purely on peak positions!
    peak_info <- peak_info %>% mutate(dist.tn=abs(position - reference_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x))) %>% mutate(dist.tp=abs(position - reference_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x)))
    peak_info <- peak_info %>% mutate(tn.peak = ifelse(dist.tn == dist.tn %>% abs() %>% min(), TRUE, FALSE)) %>% mutate(tp.peak = ifelse(dist.tp == dist.tp %>% abs() %>% min(), TRUE, FALSE))
    
    tn_candidates <- peak_info %>% arrange(peak) %>% mutate(rank_dist.tn =rank(dist.tn)) %>% mutate(rank_dist.tp =rank(-dist.tp)) %>% mutate(max_others = rank(max_others)) %>% mutate(height=rank(-as.numeric(height))) %>% mutate(rank = height + max_others + rank_dist.tn) %>% filter(rank == min(rank))
    tp_candidates <- peak_info %>% arrange(peak) %>% mutate(rank_dist.tn =rank(-dist.tn)) %>% mutate(rank_dist.tp =rank(dist.tp)) %>% mutate(max_others = rank(max_others)) %>% mutate(rank = max_others + rank_dist.tp) %>% filter(rank == min(rank))
    
    #if the same peak is classified as both TN and TP, decide based on distance to reference peaks: 
    if(nrow(tn_candidates)==1 & nrow(tp_candidates)==1 & tn_candidates$peak[[1]] == tp_candidates$peak[[1]]){
      if(tn_candidates$dist.tn > tn_candidates$dist.tp){
        tn_candidates <- tn_candidates[-1,]
      }else{
        tp_candidates <- tp_candidates[-1,]
      }
    }
  }else{
    #find candidates based on height, positin and crosstalk intensity
    tn_candidates <- peak_info %>% arrange(peak) %>% mutate(max_others = rank(max_others)) %>% mutate(position = rank(abs(as.numeric(position)))) %>% mutate(height=rank(-as.numeric(height))) %>% mutate(rank = max_others + position + height) %>% filter(rank == min(rank))
    tp_candidates <- peak_info %>% filter(!peak %in% c(tn_candidates %>% filter(position==min(position)) %>% pull(peak))) %>% arrange(peak) %>% mutate(pos = peaks$d.x[.$peak]) %>% mutate(max_others = rank(max_others)) %>% mutate(position = rank(-position)) %>% mutate(rank = max_others + position) %>% filter(rank == min(rank))
  }
  
  tn_peak <- peaks[tn_candidates %>% filter(max_others==min(max_others)) %>% filter(position==min(position)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
  
  tp_peak <- peaks[tp_candidates %>% filter(max_others==min(max_others)) %>% filter(position==min(position)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
  
  #define true positive and true negative peaks
    

  other_peaks <- peak_info %>% filter(!peak %in% (tp_candidates %>% pull(peak))) %>% filter(!peak %in% (tn_candidates %>% pull(peak)))
  
  peaks <- mutate(peaks, crosstalk.peak = FALSE)
  
  #if other peaks are found - figure out whether these are due to crosstalk!
  if(nrow(other_peaks)>0){
    crosstalk_peaks <- list()
    
    for(opeak in 1:nrow(other_peaks)){
      
      channel_candidates <- other_peaks[opeak,] %>% select(ends_with("median"))
      channel_crosstalk <- colnames(channel_candidates)[which(channel_candidates == (other_peaks[opeak,] %>% pull(max_others)))] %>% gsub("_median", "", .)
      
      if((other_peaks[opeak,] %>% pull(max_others)) > (other_peaks[opeak,] %>% pull(median_current))){
        peaks$crosstalk.peak[other_peaks[opeak,] %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric()] <- TRUE
      }else{
        print(paste("Crosstalk from", variable, "to", channel_crosstalk[[1]]))
      }
    }
  }
  if(nrow(tp_peak)>0){
    peaks <- peaks %>% mutate(tp.peak = ifelse(d.x==tp_peak$d.x, TRUE, FALSE))
  }
  if(nrow(tn_peak)>0){
    peaks <- peaks %>% mutate(tn.peak = ifelse(d.x==tn_peak$d.x, TRUE, FALSE))
  }
   
  
  tp_df_mod <- full_join(peaks, pits %>% mutate(tp.peak=FALSE) %>% mutate(tn.peak=FALSE) %>% mutate(crosstalk.peak=FALSE)) %>% arrange(d.x)
  
  return(tp_df_mod)
  }



#' @title Define a threshold based on mixture modeling
#' @description This function defines a threshold that best divides a population of data points into two separate groups based on mixture modeling.
#' @param data Either an object of class dataframe (if parameter variable is defined), of a vector that can be forced to numeric
#' @param variable  Name of the variable in the dataframe to be pulled
#' @param stringent If set to TRUE, forces script to set threshold half-way between mu of the two normal distributions. This decreases the number of false positive results.
#' @return An object of class dataframe of the RFU values from the .csv-file 

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
  
  if(!require("ggpubr", quietly=TRUE))
    install.packages("ggpubr")
  
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
  library(ggpubr)
  library(dplyr)
  library(openxlsx)
}

#' @title Analyze output of a QIAcuity experiment
#' @description This function takes the files exported from the QIAcuity software for one or more experiments and analyzes them in an automated way.
#' @param path path to where the .zip or .csv files exported from the QIAcuity software are saved
#' @param coupled_channels Data frame containing any channels the targets of which are on the same chromosome or otherwise to be expecte in the same partitions
#' @return A list containing a list of dataframes for each experiment.
#' @details This function returns a list with an element for each experiment that was in the directory set by path. In each element there are the following data frames: The data for all channels in raw data, after baseline correction, cross talk correction and competition correction; cross talk and competition analyses that give an estimate of the magnitude of the effects the different reactions have on each other; thresholds gives the thresholds for each channel at each step; Volumes contains information about the total volume contained in the partitions of each well.

QIAcuityAnalysis <- function(path, coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0))){ #}, outlier_detection=TRUE, alpha=0.01){
  setwd(path)
  
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
    # find all files with experiment name
    plate_files <- dir()[grep(plates[[k]], dir())]
    #if csv files among them, analyze those files (assumption: .zip has been unzipped!)
    csv_files <- plate_files[grep(".csv", plate_files)]
    dirname <- paste("plate", k)
    dir.create(dirname)
    
    # if no csv files, unzip any zip-files
    if(length(csv_files)==0){
      zip_file <- plate_files[grep(".zip", plate_files)]
      if(length(zip_file)==0){
        stop("Error: neither suitable .csv nor .zip file containing QIAcuity raw data were found.")
      }

      unzip(zip_file, exdir=dirname, overwrite = TRUE)
      setwd(dirname)
      csv_files <- dir()[grep(".csv", dir())]
      dataframes1 <- csv_files %>% lapply(., read.data.QIAcuity)
      setwd("..")
    }else{
      dataframes1 <- csv_files %>% lapply(., read.data.QIAcuity)
    }

    # gather data and find baseline of all samples, calculate corrected RFU values
    raw_data <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))
   
    crosstalk_corrected <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))
    competition_corrected <- data.frame(Well=character(length=0), Partition=numeric(length=0), Sample=character(length=0))
    
    for(channel in 1:length(dataframes1)){
      #extract RFU data from raw csv files
      temp <- dataframes1[[channel]]
      temp <- temp %>% data.frame() %>% mutate(RFU=RFU %>% as.numeric()) %>% group_by(Well) %>% mutate(Partition=Partition %>% as.numeric())
      
      #put those RFU values into one dataframe (that contains the data of all channels)
      raw_data <- temp %>% select(Well, Sample, Partition, RFU) %>% merge(., raw_data, all=TRUE, by=c("Well", "Sample", "Partition"))
      colnames(raw_data)[which(colnames(raw_data)=="RFU")] <- temp %>% pull(Channel) %>% unique()
    }
    
    channels <- raw_data %>% select(!Well) %>% select(!Partition) %>% select(!Sample) %>% colnames()
    
    wells <- raw_data %>% pull(Well) %>% unique()
    
    maxima <- raw_data %>% group_by(Well) %>% na.omit() %>% summarize_at(all_of(channels), max)
    
    pc_wells <- mutate(maxima, overall = (maxima %>% select(!Well) %>% apply(., 1, min))) %>% filter(overall > max(overall) - (sd(overall)/2)) %>% pull(Well)
    
    baseline_corrected <- raw_data
    
    for(channel in 1:length(channels)){
      # find true positive and negative peaks in PC
      current_channel <- channels[[channel]]
      
      tp_data <- raw_data %>% filter(Well %in% pc_wells) %>% select(!any_of(c("Well", "Partition", "Sample"))) %>% na.omit() 
      
      if(current_channel %in% coupled_channels){
        exclude <- coupled_channels %>% filter(ch1 == current_channel | ch2 ==current_channel) %>% unlist() %>% unique()
        exclude <- exclude[which(!exclude==current_channel)]
        
        tp_data <- tp_data %>% select(!any_of(exclude))
      }
      
      turnpoints <- tp_data  %>% find_turnpoints(., variable=current_channel)
      classified_turnpoints <- classify_peaks(turnpoints = turnpoints, intensities = tp_data, variable = current_channel, channel_maxima = tp_data %>% apply(., 2, max))
      
      ref_peaks <- classified_turnpoints %>% filter(tp.peak == TRUE | tn.peak == TRUE)
      
      # calculate baseline for each well
      baseline <- data.frame(Well=wells, baseline= numeric(length=length(wells)))

      for(well in 1:length(wells)){
        # try to find threshold for each well guided by positive controls, limited by crosstalk signals
        well_data <- raw_data %>% na.omit() %>% select(!any_of(exclude)) %>% filter(Well==wells[[well]])
        
        threshold_well <- well_data %>% density_threshold(., variable=current_channel, references = ref_peaks, pc_data = raw_data %>% filter(Well %in% pc_wells))
         
        n_below <- well_data %>% filter(Well==wells[[well]]) %>% filter(!!sym(current_channel)<threshold_well) %>% nrow()
        n_total <- well_data %>% filter(Well==wells[[well]])  %>% nrow()
         
        # if too few negative partitions are present, proper baseline correction is impossible
         if(n_below < n_total/10){
           baseline$baseline[[well]] <- ref_peaks %>% filter(tn.peak == TRUE) %>% pull(d.x)
         }else{
                      
           baseline$baseline[[well]] <- well_data %>% filter(!!sym(current_channel)<threshold_well) %>% pull(current_channel) %>% median()
         }
      }
      
      baseline_corrected <- baseline_corrected %>% select(any_of(c("Well", "Partition", "Sample", current_channel))) %>% merge(., baseline, by=c("Well")) %>% mutate(!!sym(current_channel) := !!sym(current_channel) - baseline) %>% select(!baseline) %>% full_join(baseline_corrected %>% select(!any_of(current_channel)))
      
      
    }

    #check for correlation between channels to prevent later errors!
  
    cors <- raw_data %>% filter(Well %in% pc_wells) %>% select(!Well) %>% select(!Partition) %>% select(!Sample) %>% na.omit() %>% cor()
    cors <- cors>0.4 & cors<1
    
    row <-  (which(cors==TRUE)/ncol(cors)) %>% ceiling(.)
    row[which(row == 0)] <- 1
    col <-  which(cors==TRUE)%%ncol(cors)
    col[which(col == 0)] <- ncol(cors)
    
    pairs <- data.frame(ch1=colnames(cors)[col], ch2=rownames(cors)[row])
    
    for(i in 1:nrow(pairs)){
      pairs[i,] <- pairs[i,][order(pairs[i,])]
    }
    
    pairs <- pairs %>% unique() %>% na.omit() %>% data.frame()
    
    #re-define thresholds for each channel that separate positive and negative samples 

    thresholds <- data.frame(channel=channels, baseline=numeric(length=length(channels)), crosstalk = numeric(length=length(channels)), competition=numeric(length=length(channels)))
    
    baseline_corrected_dichot <- baseline_corrected

    for(i in 1:length(channels)){
      if(channels[[i]] %in% coupled_channels){
        omit <- coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]]) %>% unlist() %>% unique()
        comparison_channels <- channels[which(!channels %in% omit)]
        
      }else{
        comparison_channels <- channels[-i]
      }
      
      threshold <- baseline_corrected %>% filter(Well %in% pc_wells) %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.2)
      
      thresholds[which(thresholds$channel ==channels[[i]]), "baseline"]<- threshold    
      
      # make dichotomized dataset that just gives info on whether a partition is positive or not
      baseline_corrected_dichot[,channels[[i]]] <- ifelse(baseline_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }
    
    # crosstalk correction using a linear model
    
    crosstalk_corrected <- baseline_corrected
    crosstalk_analysis = data.frame(ch1=character(length=0), ch2=character(length=0), cross_talk=numeric(length=0))
    
    for(i in 1:length(channels)){
      
      for(j in 1:length(channels)){
                
        if(!i==j){
          #PROBLEM: if channels are coupled (e.g. targets on same chromosome), the filtering would lead to massive loss of data; not filtering will lead to less clear data!
          excluded=c()
          if(nrow(coupled_channels)>0){
            if(channels[[i]] %in% (coupled_channels %>%  unlist())){
              excluded <- coupled_channels %>% filter(ch1==channels[[i]] | ch2==channels[[i]]) %>% unlist() %>% unique() %>% c(excluded, .)
            }
            if(channels[[j]] %in% (coupled_channels %>% unlist())){
              excluded <- coupled_channels %>% filter(ch1==channels[[j]] | ch2==channels[[j]]) %>% unlist() %>% unique() %>% c(excluded, .)
            }
            
            if(channels[[i]] %in% excluded && channels[[j]] %in% excluded){
              next()
            }
            
            if(channels[[j]] %in% coupled_channels){
              #PROBLEM: if two channels are coupled and exert crosstalk, the script can interpret it as crosstalk from both!
              # therefore, always only analyze the channel that is closer wavelength-wise
              
              channel_order <- c("G", "Y", "O", "R", "C")
              current_channel_pos <- which(channel_order ==channels[[i]])
              comparison_channel_pos <- which(channel_order ==channels[[j]])
              comp_coupled <- coupled_channels %>% filter(ch1==channels[[j]] | ch2 == channels[[j]]) %>% unlist() 
              comp_coupled <- comp_coupled[which(!comp_coupled==channels[[j]])]
              
              coupled_channel_pos <- which(channel_order == comp_coupled)
              
              if(abs(current_channel_pos - comparison_channel_pos)>abs(current_channel_pos - coupled_channel_pos)){
                next()
              }
            }
            
            excluded <- excluded[which(!excluded %in% channels[c(i,j)])] %>% unique()
            
          }
          if(length(excluded)==0){
            excluded=c("PLACEHOLDER")
          }
          
          new_thresh_data <- baseline_corrected_dichot %>% select(!any_of(channels[[i]])) %>% select(!any_of(channels[[j]])) %>% select(!any_of(excluded)) %>% mutate(., sum := !!. %>% select(any_of(channels)) %>% rowSums(., na.rm = TRUE))  %>% filter(sum==0 | is.na(sum)==TRUE)  %>% select(Well, Partition) %>% left_join(., baseline_corrected) %>% na.omit()
          
          crosstalk_threshold <- new_thresh_data %>% select(!any_of(excluded)) %>% filter(Well %in% pc_wells) %>% density_threshold(input_data = ., variable=channels[[i]])
          
          negs <- new_thresh_data %>% filter(!!sym(channels[[i]])<crosstalk_threshold)
          
          sp_points <- negs %>% filter(!!sym(channels[[j]]) > thresholds %>% filter(channel==channels[[j]]) %>% pull(baseline))
          
          if(nrow(sp_points)<100){
            print("Warning: fewer than 100 positive points found. Aborting crosstalk calculation")
            next()
          }
          
          dn_points <- negs %>% filter(!!sym(channels[[j]]) < thresholds %>% filter(channel==channels[[j]]) %>% pull(baseline))
          
          dn_med <- c(dn_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), dn_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          sp_med <- c(sp_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), sp_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          
          slope <- (sp_med[[1]] - dn_med[[1]])/(sp_med[[2]] - dn_med[[2]])
          yintercept <- dn_med[[1]] - slope*dn_med[[2]]
          
          crosstalk <- (baseline_corrected %>% pull(channels[[j]]))*slope + yintercept #predict(model_crosstalk, baseline_corrected)*(1-(corr_factor-1)^2)
          crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=slope) %>% bind_rows(., crosstalk_analysis)
            
          crosstalk_corrected[,channels[[i]]] <- crosstalk_corrected[,channels[[i]]] - crosstalk
        }
      }
    }
    
    crosstalk_corrected_dichot <- crosstalk_corrected
    
    for(i in 1:length(channels)){
      if(channels[[i]] %in% coupled_channels){
        omit <- coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]]) %>% unlist() %>% unique()
        comparison_channels <- channels[which(!channels %in% omit)]
        
      }else{
        comparison_channels <- channels[-i]
      }
      
      threshold <- crosstalk_corrected %>% filter(Well %in% pc_wells)  %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.2)
      
      thresholds[which(thresholds$channel ==channels[[i]]), "crosstalk"]<- threshold    
      
      # make dichotomized dataset that just gives info on whether a partition is positive or not
      crosstalk_corrected_dichot[,channels[[i]]] <- ifelse(crosstalk_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }
    
    
    # perform competition correction
    competition_corrected <- crosstalk_corrected
    competition_analysis <- list()
    
    for(i in 1:length(channels)){
      # for each channel, do:
      current_channel <- channels[[i]]
      # find data from PC qwells
      
      comp_data <- competition_corrected %>% filter(Well %in% pc_wells) %>% select(!any_of(c("Well", "Partition", "Sample"))) %>% na.omit()
      
      # find peaks
      peaks <- find_turnpoints(comp_data, variable = current_channel)
      # classify them
      classified_peaks <- peaks %>% classify_peaks(., comp_data, variable = current_channel, channel_maxima = comp_data %>% apply(., 2, max)) 
      # only examine peaks between TP and TN peak
      rel_tps <- classified_peaks %>% filter(d.x > classified_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x)) %>% filter(d.x < classified_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x))
      if(nrow(rel_tps)>1){
        breaks <- rel_tps %>% filter(tp.pits ==TRUE) %>% pull(d.x)
        
        current_baseline <- classified_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x)
        
        # summarize; find relevant competition-induced peaks!
        peak_summary <- comp_data %>% mutate(peak:= !!sym(current_channel) %>% cut(., breaks=breaks, labels=1:(length(breaks)-1))) %>% group_by(peak) %>% summarize_at(.vars=colnames(comp_data), median) %>% na.omit()
        # find partitions that should
        comp_correction_candidates <- competition_corrected %>%  mutate(peak:= !!sym(current_channel) %>% cut(., breaks=breaks, labels=1:(length(breaks)-1))) %>% na.omit() %>% select(all_of(c("Partition", "Well", "peak"))) %>% left_join(., crosstalk_corrected_dichot)
        
        neg_sd <- crosstalk_corrected %>% filter(R< thresholds[which(thresholds$channel == channels[[i]]), "crosstalk"]) %>% pull(R) %>% sd(., na.rm=TRUE)
        neg_mean <- crosstalk_corrected %>% filter(R< thresholds[which(thresholds$channel == channels[[i]]), "crosstalk"]) %>% pull(R) %>% mean(., na.rm=TRUE)

        comp_correction_candidates <- comp_correction_candidates %>% filter(if_any(ends_with(channels[which(!channels==current_channel)]), ~ . >0))
        
        comp_correction_candidates <- merge(comp_correction_candidates, competition_corrected, by=c("Partition", "Well"))
        
        vars <- comp_correction_candidates  %>% select(matches("^\\w\\.x")) %>% colnames()
        
        comp_correction_candidates <- comp_correction_candidates %>% mutate(z_score = abs(!!sym(paste(channels[[i]], ".y", sep="")) - neg_mean)/neg_sd) %>% filter(z_score >6)
        
        comp_correction_candidates <- comp_correction_candidates %>% group_by(across(vars))
        
        competition_summary <- comp_correction_candidates %>% add_count() %>% summarize_at(.vars=c("n", paste(current_channel,".y", sep="")), median)
        
        competition_summary <- competition_summary %>% mutate(strength = !!sym(paste(current_channel,".y", sep=""))/(classified_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x))) %>% filter(n>10) %>% select(!n)
        
        comp_correction_candidates <- comp_correction_candidates %>% inner_join(., competition_summary %>% select(!ends_with(".y")), by=competition_summary %>% select(ends_with(".x")) %>% colnames()) %>% mutate(!!sym(paste(current_channel,".y", sep="")) := !!sym(paste(current_channel,".y", sep=""))/strength)
        
        comp_corrected_current <- comp_correction_candidates %>% ungroup() %>% select(any_of(c("Partition", "Well", "Sample.y", paste(channels, ".y", sep=""))))
        
        colnames(comp_corrected_current) <- comp_corrected_current %>% colnames() %>% gsub(".y", "", .)
        
        comp_corrected_current <- comp_corrected_current %>% relocate(colnames(competition_corrected))
        
        competition_corrected <- competition_corrected %>% rows_update(., comp_corrected_current, by=c("Partition", "Well"))
        
        competition_analysis[[i]] <- competition_summary
      }
      
    }
    
    competition_analysis <- bind_rows(competition_analysis)
    
    # re-calculate thresholds
    
    competition_corrected_dichot <- competition_corrected
    
    for(i in 1:length(channels)){
      # define which channels the data should be compared to
      if(channels[[i]] %in% coupled_channels){
        omit <- coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]]) %>% unlist() %>% unique()
        comparison_channels <- channels[which(!channels %in% omit)]
        
      }else{
        comparison_channels <- channels[-i]
      }
      
      #calculate new threshold for this channel
      threshold <- competition_corrected %>% filter(Well %in% pc_wells) %>% mutate(Sample=Well) %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = 0.5)
      
      thresholds[which(thresholds$channel ==channels[[i]]), "competition"]<- threshold    
      
      # make dichotomized dataset that just gives info on whether a partition is positive or not
      competition_corrected_dichot[,channels[[i]]] <- ifelse(competition_corrected_dichot[,channels[[i]]] > threshold, 1, 0)
    }
    
    output <- list(raw_data = raw_data, baseline_corrected = baseline_corrected, crosstalk_corrected=crosstalk_corrected, competition_corrected=competition_corrected, crosstalk_analysis=crosstalk_analysis, competition_analysis=competition_analysis, thresholds=thresholds)
    
    # calculate volume
    vol_by_well = dataframes1 %>% bind_rows() %>% group_by(Well, Channel) %>% select(Well, starts_with("Cycled"), Channel)%>% unique()
    output$Volumes <- vol_by_well
    plate_results[[plates[[k]]]] <- output
  }
  
  return(plate_results)
}

#' @title Generate output of QIAcuity Analyisis
#' @description  This function summarizes and visualizes the results of one or multiple QIAcuity experiments as.
#' @param results The output of the QIAcuityAnalysis function.
#' @param detailed_plots Set to TRUE to generate 2D scatterplots at every correction step
#' @param scatterplots_1d Set to TRUE to generate 1D scatterplots at every correction step
#' @return NULL. Generates output files: .csv files for the results per partition, .xlsx files for summaries at each correction step. 2D scatterplots for all channel combinations and 1D scatterplots for each channel at every step. Also puts the crosstalk and competition analyses into .xlsx files. 

summarizeQIAcuityResult <- function(results, detailed_plots=TRUE, scatterplots_1d=TRUE){
  plate_summaries <- rep(NA, times=length(results)) %>% as.list()
  names(plate_summaries) <- names(results)
  corrections = c("raw_data", "baseline_corrected", "crosstalk_corrected", "competition_corrected")
  for(i in 1:length(results)){
    dirname <- paste("plate", i)
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
  
      n_part <- channel_summaries %>% mutate(valid=positive + negative) %>% pull(valid)# %>% median(., na.rm = TRUE)
      
      mm_vol <- 0
      
      if(n_part %>% median(., na.rm = TRUE) > 10000){
        mm_vol <- 40/26000*n_part
      }else{
        mm_vol <- 12/8500*n_part        
      }
      
      poisson.CI <- function(x){return(poisson.test(x)$conf.int %>% round(digits=1) %>% paste(., collapse=" - "))}      
      
      
      
      channel_summaries <- channel_summaries %>% mutate(CI=lapply(.$positive, poisson.CI)) %>% relocate(any_of(c("Sample", "Well", "Channel", "negative", "positive", "CI"))) 
      
      channel_summaries <- channel_summaries %>% mutate(Warning = ifelse(lapply(.$CI %>% as.character() %>% strsplit(split=" - "), min)<1 & .$positive > 0, "Warning: confidence interval contains values < 1, positive result may not be reliable.", NA) %>% paste(channel_summaries$Warning, sep=", ") %>% gsub("NA, NA", NA, .) %>% gsub("NA, ", "", .) %>% gsub(", NA", "", .))
      
      channel_summaries <- merge(channel_summaries, results[[i]]$Volumes) %>% mutate(factor = (!!sym("Cycled volume") %>% as.numeric())/mm_vol)
      
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
    
    setwd("..")
  }
}
