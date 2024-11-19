#' @title Eliminate background signal and artifacts
#' @description  This function eliminates artifacts by searching for break points where sudden changes in fluorescence intensity on the plate occur, and fits a model to each region identified to estimate artifact effects and subtract them from the data.
#' @param data The fluorescence data from QIAcuity
#' @param current_channel Which column the artifact correction should be performed on.
#' @param upper_lim Upper RFU limit of which partitions should be considered as part of the negative group
#' @param lower_lim Lower RFU limit of which partitions should be considered as part of the negative group
#' @return Returns data in the same format, with only the data in column current_channel altered to remove artifacts
#'
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

#' @title Eliminate background signal and artifacts across a whole experiment
#' @description  This function eliminates artifacts and background signal from QIAcuity data from a whole experiment in multiple steps.
#' @param raw_data The raw fluorescence data from QIAcuity
#' @param smooth Boolean; default:TRUE. Whether or not the function should search for break points in the data and eliminate any artifacts that skew fluorescence intensities.
#' @param coupled_channels Defines which channels are highly correlated.
#' @param channels Defines which columns contain fluorescence data.
#' @param pc_wells Defines which wells are the positive controls
#' @return Returns data in the same format as raw_data, with background signal and artifacts removed from data.

baseline_correction <- function(raw_data, smooth=TRUE, coupled_channels, channels, pc_wells, thresholds){
  baseline_corrected <- raw_data
  for(channel in 1:length(channels)){
    # find true positive and negative peaks in PC
    current_channel <- channels[[channel]]
    tp_data <- raw_data %>% filter(Well %in% pc_wells) %>% select(!any_of(c("Well", "Partition", "Sample"))) %>% na.omit()
    wells <- raw_data %>% pull(Well) %>% unique()

    baseline <- data.frame(Well=character(length=0), Partition=character(length=0), baseline= numeric(length=0))

    if(current_channel %in% coupled_channels){
      exclude <- coupled_channels %>% filter(ch1 == current_channel | ch2 ==current_channel) %>% unlist() %>% unique()
      exclude <- exclude[which(!exclude==current_channel)]

      tp_data <- tp_data %>% select(!any_of(exclude))
    }else{
      exclude <- c()
    }

    turnpoints <- tp_data  %>% find_turnpoints(., variable=current_channel)
    classified_turnpoints <- classify_peaks(turnpoints = turnpoints, intensities = tp_data, variable = current_channel, channel_maxima = tp_data %>% apply(., 2, get_quantile, probs=0.999), thresholds=thresholds, step="raw", coupled_channels=coupled_channels)

    channel_maxima <- tp_data %>% apply(., 2, get_quantile, probs=0.999)

    ref_peaks <- classified_turnpoints %>% filter(tp.peak == TRUE | tn.peak == TRUE)
    if((ref_peaks %>% filter(tp.peak==TRUE) %>% nrow())==0){
      tp_artif <- classified_turnpoints %>% filter(tp.peaks==TRUE) %>% filter(d.x > ref_peaks %>% filter(tn.peak=TRUE) %>% pull(d.x)) %>% filter(d.x==min(d.x)) %>% mutate(tp.peak=TRUE)
      ref_peaks <- ref_peaks %>% rbind(., tp_artif)
    }
    if((ref_peaks %>% filter(tn.peak==TRUE) %>% nrow())==0){
      # include three largest peaks as options only
      tn_artif <- classified_turnpoints %>% filter(tp.peaks==TRUE) %>% mutate(d.y.rank = rank(-d.y)) %>% filter(d.y.rank < 4) %>% filter(d.x < ref_peaks %>% filter(tp.peak=TRUE) %>% pull(d.x)) %>% filter(d.x==min(d.x)) %>% mutate(tn.peak=TRUE) %>% select(!d.y.rank)
      ref_peaks <- ref_peaks %>% rbind(., tn_artif)
    }

    # calculate baseline for each well

    for(well in 1:length(wells)){
      # try to find threshold for each well guided by positive controls, limited by crosstalk signals
      well_data <- raw_data %>% na.omit() %>% select(!any_of(exclude)) %>% filter(Well==wells[[well]])

      turnpoints_well <- well_data %>% select(!any_of(c("Well", "Sample", "Partition"))) %>% na.omit() %>% find_turnpoints(., variable=current_channel, n_points=15)
      classified_turnpoints_well <- classify_peaks(turnpoints = turnpoints_well, intensities = well_data %>% select(!any_of(c("Well", "Sample", "Partition"))), variable = current_channel, channel_maxima = channel_maxima, reference_peaks = ref_peaks, thresholds = thresholds, step="raw", coupled_channels = coupled_channels)

      tn_peak_well <- classified_turnpoints_well %>% filter(tn.peak == TRUE)

      if(nrow(tn_peak_well)>0){
        upper_lim <- classified_turnpoints_well %>% filter(d.x > tn_peak_well %>% pull(d.x)) %>% pull(d.x) %>% min()
        lower_lim <- classified_turnpoints_well %>% filter(d.x < tn_peak_well %>% pull(d.x)) %>% pull(d.x) %>% max()
        datapoints <- well_data %>% filter(!!sym(current_channel) > lower_lim) %>% filter(!!sym(current_channel) < upper_lim) %>% nrow()
      }else{
        datapoints <- 0
      }

      if(nrow(tn_peak_well) == 0){
        baseline <- data.frame(Well=well_data$Well, Partition = well_data$Partition, baseline = ref_peaks %>% filter(tn.peak == TRUE) %>% pull(d.x)) %>% rbind(baseline, .)
        warning(paste("Well", wells[[well]], ", Channel",  current_channel, ". Baseline estimation potentially failed. Please inspect data."))
      }else{
        if(datapoints<(nrow(well_data)/10)){
          baseline <- data.frame(Well=well_data$Well, Partition = well_data$Partition, baseline = classified_turnpoints_well %>% filter(tn.peak == TRUE) %>% pull(d.x)) %>% rbind(baseline, .)
          warning(paste("Well", wells[[well]], ", Channel",  current_channel, ". Baseline estimation potentially failed. Please inspect data."))
        }else{
          if(smooth==TRUE){
            baseline <- smooth_data(well_data, upper_lim = upper_lim, lower_lim= lower_lim, current_channel=current_channel) %>% select(Well, Partition, baseline) %>% rbind(baseline, .)
            #baseline %>% filter(Well == wells[[well]]) %>% pull(baseline) %>% mean() %>% print()
          }else{
            baseline_well <- well_data %>% filter(!!sym(current_channel) > lower_lim) %>% filter(!!sym(current_channel) < upper_lim) %>% select(any_of(current_channel)) %>% unlist() %>% median()
            baseline <- data.frame(Well=well_data$Well, Partition = well_data$Partition, baseline = baseline_well) %>% rbind(baseline, .)
          }
        }

      }
    }

    baseline_corrected <- baseline_corrected %>% select(any_of(c("Well", "Partition", "Sample", current_channel))) %>% merge(., baseline, by=c("Well", "Partition")) %>% mutate(!!sym(current_channel) := !!sym(current_channel) - baseline) %>% select(!baseline) %>% full_join(baseline_corrected %>% select(!any_of(current_channel))) %>% select(Well, Partition, Sample, all_of(channels))# %>% return()

    }

  return(baseline_corrected)
}

#' @title Calculate thresholds for input data and dichotomize it
#' @description  This function takes QIAcuity fluorescence data, calculates thresholds for all channels and then classifies each partition based on these thresholds as 0 or 1.
#' @param data The fluorescence data from QIAcuity
#' @param thresholds Data frame where threshold values are stored.
#' @param step Column of the data frame where the data should be stored.
#' @param coupled_channels Defines which channels are highly correlated.
#' @param min_dist Minimum distance of the threshold to negative band relative to difference between negative and positive band (0.2 = 20% of the difference between the negative and positive band is the minimum distance between negative band and threshold)
#' @return Returns a list of (1) the thresholds and (2) the data, transformed to 0 or 1 depending on whether the fluorescence value exceeds the threshold

recalculate_thresholds <- function(data, thresholds, step, coupled_channels, min_dist=0.2, mid=FALSE, pc_wells, noise_suppression_strict =TRUE){
  channels <- thresholds$channel
  #maxima <- data %>% group_by(Well) %>% na.omit() %>% summarize_at(channels, get_quantile, probs=0.999)
  #pc_wells <- maxima %>% mutate(across(.cols=any_of(channels), .fns=function(x){return(x-min(x))})) %>% mutate(across(.cols=any_of(channels), .fns=function(x){return(x/max(x))}))  %>% mutate(minimum = apply(across(channels), 1, min)) %>% mutate(median = apply(across(channels), 1, median)) %>% filter(median == max(median) | minimum == max(minimum) | median == min(median)) %>% pull(Well) %>% unique()
  corrs <- colnames(thresholds)[-1]
  position = which(corrs==step)-1
  if(position ==0){
    position <- 1
  }
  prev_step <- corrs[position]
  if(step %in% c("baseline", "raw")){
    recalc=TRUE
  }else{
    recalc=FALSE
  }

  data_dichot <- data

  for(i in 1:length(channels)){
    if(step=="crosstalk"){
      comparison_channels <- c()
    }else{
      comparison_channels <- channels[-i]
      comparison_channels <- comparison_channels[which(!comparison_channels %in% (coupled_channels  %>% filter(ch1==channels[[i]]) %>% unlist() %>% unique()))]
    }

    threshold <- data %>% filter(Well %in% pc_wells) %>% select(any_of(c(channels[[i]], comparison_channels))) %>%  density_threshold(., variable=channels[[i]], min_dist = min_dist, mid=mid, thresholds = thresholds, step=prev_step, coupled_channels=coupled_channels, recalc = recalc, noise_suppression_strict = noise_suppression_strict)

    thresholds[which(thresholds$channel ==channels[[i]]), step]<- threshold

    # make dichotomized dataset that just gives info on whether a partition is positive or not
    data_dichot[,channels[[i]]] <- ifelse(data[,channels[[i]]] > threshold, 1, 0)
  }

  return(list(data_dichot, thresholds))
}

#' @title Perform crosstalk correction
#' @description  This function takes QIAcuity fluorescence data, calculates the effects of crosstalk between the channels and corrects for them
#' @param baseline_data The baseline corrected fluorescence data from QIAcuity
#' @param baseline_data_dichot The baseline corrected fluorescence data from QIAcuity, dichotomized
#' @param pc_wells Wells that contain the positive controls
#' @param coupled_channels Defines which channels are highly correlated.
#' @param channels Defines which columns in the data contain fluorescence data
#' @param thresholds Data frame with the thresholds for all of the channels
#' @return Returns a list of (1) the crosstalk corrected data and (2) an overview of the crosstalk effects between the different channels

crosstalk_correction <- function(baseline_data, baseline_data_dichot, pc_wells, coupled_channels, channels, thresholds){
  crosstalk_corrected <- baseline_data
  crosstalk_analysis = data.frame(ch1=character(length=0), ch2=character(length=0), cross_talk=numeric(length=0))

  for(i in 1:length(channels)){
    for(j in 1:length(channels)){
      if(!i==j){
        #if(channels[[j]] %in% (coupled_channels %>% filter(ch1==channels[[i]] | ch2 == channels[[i]] )  %>% unlist() %>% unique())){
        #  print("Warning: coupled channels. Aborting crosstalk correction")
        #  crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=NA) %>% bind_rows(., crosstalk_analysis)
        #  next()
        #}

        pc_data <- crosstalk_corrected %>% filter(Well %in% pc_wells)

        filtered_points <- baseline_data_dichot %>% filter(Well %in% pc_wells)
        excluded_channels <- coupled_channels %>% filter(ch1 %in% channels[c(i,j)] | ch2 %in% channels[c(i,j)]) %>% unlist() %>% unique()
        excluded_channels <- excluded_channels[which(!excluded_channels %in% c(channels[[j]], channels[[i]]))]

        filtered_points$sum <- filtered_points %>% select(any_of(channels)) %>% select(!any_of(channels[c(i,j)])) %>% select(!any_of(excluded_channels)) %>% apply(., 1, sum)
        filtered_points <- filtered_points %>% filter(sum <= median(sum, na.rm =TRUE))

        pc_data <- filtered_points %>% select(Well, Sample, Partition) %>% left_join(., pc_data)

        # re-define thresholds to potentially eliminate noise from competition effects
        thresh_calc_data <- pc_data %>% filter(!!sym(channels[[i]]) < thresholds %>% filter(channel == channels[[j]]) %>% pull(baseline))
        if(nrow(thresh_calc_data)>0){
          thresh_temp <- recalculate_thresholds(thresh_calc_data, thresholds, step = "crosstalk", pc_wells = pc_wells, coupled_channels= coupled_channels)
        }else{
          thresh_temp <- list("", thresholds %>% mutate(crosstalk = baseline))
        }

        thr_ch2 <- thresh_temp[[2]] %>% filter(channel == channels[[j]]) %>% pull(crosstalk)
        if(thr_ch2 < 3){
          thr_ch2 <- thresholds %>% filter(channel == channels[[j]]) %>% pull(baseline)
        }
        if(thr_ch2 < 3){
          thr_ch2 <- 3
        }
        thr_ch1 <- thresh_temp[[2]] %>% filter(channel == channels[[i]]) %>% pull(crosstalk)
        if(thr_ch1 < 3){
          thr_ch1 <- thresholds %>% filter(channel == channels[[i]]) %>% pull(baseline)
        }
        if(thr_ch1 < 3){
          thr_ch1 <- 3
        }
        sp_points <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% filter(!!sym(channels[[i]]) < thr_ch1)

        if(nrow(sp_points)<20){
          thr_ch2 <- thresholds %>% filter(channel == channels[[j]]) %>% pull(baseline)
          thr_ch1 <- thresholds %>% filter(channel == channels[[i]]) %>% pull(baseline)
          sp_points <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% filter(!!sym(channels[[i]]) < thr_ch1)
          }

        sp_points_ch2 <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2) %>% filter(!!sym(channels[[i]]) > thr_ch1)

        dp_points <- pc_data %>% filter(!!sym(channels[[i]]) > thr_ch1) %>% filter(!!sym(channels[[j]]) > thr_ch2)

        expected_dps <- nrow(pc_data) * ((nrow(sp_points) + nrow(dp_points))/nrow(pc_data)) * ((nrow(sp_points_ch2) + nrow(dp_points))/nrow(pc_data))

        dn_points <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2)  %>% filter(!!sym(channels[[i]]) < thr_ch1)

        if(nrow(sp_points)<20){
          if(nrow(dp_points)>nrow(sp_points)*3 & nrow(sp_points_ch2 > nrow(sp_points)*3)){
            thr_ch1 <- sp_points_ch2 %>% pull(channels[[i]]) %>% median(na.rm=TRUE)
            sp_points <- pc_data %>% filter(!!sym(channels[[j]]) > thr_ch2) %>% filter(!!sym(channels[[i]]) < thr_ch1)
            sp_points_ch2 <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2) %>% filter(!!sym(channels[[i]]) > thr_ch1)

            dp_points <- pc_data %>% filter(!!sym(channels[[i]]) > thr_ch1) %>% filter(!!sym(channels[[j]]) > thr_ch2)

            expected_dps <- nrow(pc_data) * ((nrow(sp_points) + nrow(dp_points))/nrow(pc_data)) * ((nrow(sp_points_ch2) + nrow(dp_points))/nrow(pc_data))

            dn_points <- pc_data %>% filter(!!sym(channels[[j]]) < thr_ch2)  %>% filter(!!sym(channels[[i]]) < thr_ch1)

          }else{
            print("Warning: fewer than 20 positive points found. Aborting crosstalk calculation")
            crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=NA) %>% bind_rows(., crosstalk_analysis)
            next()
          }
        }

        if(nrow(dn_points)<100){
            print("Warning: fewer than 100 negative points found. Aborting crosstalk calculation")
            crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=NA) %>% bind_rows(., crosstalk_analysis)
            next()
        }

        if(nrow(dp_points)>expected_dps*100){
          coupled_channels <- data.frame(ch1=channels[[i]], ch2=channels[[j]]) %>% rbind(coupled_channels)
        }


          dn_med <- c(dn_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), dn_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          sp_med <- c(sp_points %>% pull(channels[[i]]) %>% median(na.rm=TRUE), sp_points %>% pull(channels[[j]]) %>% median(na.rm=TRUE))
          if(!is.na(sp_med[[1]])  & !is.na(dn_med[[1]])){

          slope <- (sp_med[[1]] - dn_med[[1]])/(sp_med[[2]] - dn_med[[2]])

          if(slope>1){
            break()
          }

          yintercept <- dn_med[[1]] - slope*dn_med[[2]]

          crosstalk <- (crosstalk_corrected %>% pull(channels[[j]]))*slope + yintercept #predict(model_crosstalk, baseline_corrected)*(1-(corr_factor-1)^2)
          crosstalk_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], cross_talk=slope) %>% bind_rows(., crosstalk_analysis)

          crosstalk_corrected[,channels[[i]]] <- crosstalk_corrected[,channels[[i]]] - crosstalk
          }
      }
    }
  }
  return(list(crosstalk_corrected, crosstalk_analysis))
}

#' @title Perform competition correction
#' @description  This function takes QIAcuity fluorescence data, calculates the effects of competition between the different reactions and corrects for them
#' @param crosstalk_corrected The crosstalk corrected fluorescence data from QIAcuity
#' @param pc_wells Wells that contain the positive controls
#' @param channels Defines which columns in the data contain fluorescence data
#' @param thresholds Data frame with the thresholds for all of the channels
#' @return Returns a list of (1) an overview of the competition effects between the different channels and (2) the competition corrected data

competition_correction <- function(crosstalk_corrected, channels, pc_wells, thresholds, step){
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
        competition[which(is.na(competition))] <- 0
        competition_analysis <- data.frame(ch1=channels[[j]], ch2=channels[[i]], competition=slope) %>% bind_rows(., competition_analysis)

        competition_corrected[,channels[[i]]] <- competition_corrected[,channels[[i]]] - competition

      }
    }
  }
  return(list(competition_analysis, competition_corrected))
}
