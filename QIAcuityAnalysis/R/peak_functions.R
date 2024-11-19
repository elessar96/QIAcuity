
#' @title Define a threshold based on mixture modeling
#' @description This function defines a threshold that best divides a population of data points into two separate groups based on mixture modeling.
#' @param input_data An object of class dataframe containing fluorescence intensities fom ddPCR partitions.
#' @param variable  Name of the variable in the dataframe to be pulled
#' @param min_dist Determines the minimum distance between threshold and negative partitions, always as a fraction of the distance between true positive and true negative populations.
#' @param references Reference peaks for guiding classification
#' @param pc_data Positive control data as a reference to calculate maximum expected values for each channel
#' @return A numeric value giving the optimal threshold for distinguishing the true positive and true negative populations in the data.

density_threshold <- function(input_data, variable, min_dist=0.2, references=NULL, pc_data =NULL, mid=FALSE, thresholds, step, coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0)), recalc = FALSE, noise_suppression_strict =TRUE){

  data_subset <- input_data %>% select(!any_of(c("Well", "Sample", "Partition"))) %>% na.omit()

  # find turnpoints
  tp_df <- find_turnpoints(data_subset, variable=variable)

  if(is.null(pc_data)){
    pc_data <- data_subset
  }

  pc_data <- pc_data %>% select(!any_of(c("Well", "Sample", "Partition")))

  channel_maxima <- pc_data %>% apply(., 2, get_quantile, probs=0.999)
  channel_minima <- pc_data %>% apply(., 2, get_quantile, probs=0.001)
  if((channel_maxima/channel_minima) %>% min() <1.5){
    channel_maxima <- pc_data %>% apply(., 2, get_quantile, probs=1-(5/26000))
  }
  # find true positives and negatives
  classified_turnpoints <- classify_peaks(tp_df, intensities = data_subset , channel_maxima = channel_maxima, variable=variable, reference_peaks = references, thresholds = thresholds, step = step, coupled_channels=coupled_channels, recalc, noise_suppression_strict = noise_suppression_strict)

  #determine position of threshold
  positions <- data.frame(tp=NA, tn=NA, ct=NA)
  positions$tp[[1]] <- classified_turnpoints %>% filter(tp.peak == TRUE) %>% pull(d.x) %>% max()
  positions$tn[[1]] <- classified_turnpoints %>% filter(tn.peak == TRUE) %>% pull(d.x) %>% max()

  if(!is.finite(positions$tp)){
    positions$tp <- classified_turnpoints %>%
      filter(tp.peaks == TRUE & tn.peak == FALSE & !ct_peak ==TRUE) %>%
      pull(d.x) %>%
      max()
  }
  if(!is.finite(positions$tn)){
    positions$tn <- classified_turnpoints %>%
      filter(tp.peaks == TRUE & tp.peak == FALSE & !ct_peak ==TRUE) %>%
      pull(d.x) %>%
      min()
  }

  pot_crosstalk_peaks <- classified_turnpoints
  if(is.numeric(positions$tp)){
    pot_crosstalk_peaks <- pot_crosstalk_peaks %>% filter(d.x < positions$tp)
  }
  if(is.numeric(positions$tn)){

    pot_crosstalk_peaks <- pot_crosstalk_peaks %>% filter(d.x > positions$tn)

  }
  positions$ct[[1]] <-   pot_crosstalk_peaks %>% filter(ct_peak == TRUE) %>% pull(d.x) %>% max()

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

  if(mid==TRUE){
    if(is.finite(positions$tp)){
      threshold <- positions$tn + (positions$tp-positions$tn)*0.5
    }else{
      threshold <- positions$tn + ((classified_turnpoints %>% filter(tp.peaks == TRUE) %>% filter(d.x > positions$tn) %>% filter(d.y > max(d.y)/3) %>% filter(d.x == max(d.x, na.rm=TRUE)) %>% pull(d.x)) + positions$tn)*0.5
    }

  }

  return(threshold)

}


#' @title Classify peaks previously detected by find_turnpoints
#' @description This function classifies turnpoints in a distribution
#' @param turnpoints Output of find_turnpoints
#' @param intensities Fluorescence intensity values in a data frame
#' @param variable Which column in intensities is the channel to be examined
#' @param channel_maxima A named list giving the maximum values expected for each of the channels found in intensities
#' @param reference_peaks A dataframe of previously classified peaks to guide this classification; optional.
#' @return A dataframe with turnpoints


classify_peaks <- function(turnpoints, intensities, variable, channel_maxima, reference_peaks=NULL, thresholds, step, coupled_channels=data.frame(ch1=character(length=0), ch2=character(length=0)), recalc = FALSE, noise_suppression_strict = TRUE){
  # find peaks and pits
  peaks <- turnpoints %>% filter(tp.peaks) %>% arrange(d.x)
  pits <- turnpoints %>% filter(tp.pits)

  peak_info <- mutate(intensities, peak=cut(intensities %>% pull(variable) %>% as.numeric(),breaks=pits$d.x, labels=paste("peak", c(1:nrow(peaks))))) %>% na.omit()

  #normalize intensity data
  channels <- colnames(intensities)
  for(chan in 1: length(channels)){
    peak_info[,channels[[chan]]] <- peak_info[,channels[[chan]]]/channel_maxima[channels[[chan]]]
  }
  if(noise_suppression_strict ==TRUE){

    size_threshold <- peak_info %>%
      na.omit() %>%
      group_by(peak) %>%
      summarize(across(colnames(intensities), function(x){median(x, na.rm=TRUE) %>% return()}), n=n()) %>%
      na.omit() %>%
      mutate_if(is.numeric, round, digits=1) %>%
      filter(n < max(n)) %>%
      filter(n == max(n)) %>%
      pull(n) %>%
      unique()/100

    if(length(size_threshold) ==0){
      size_threshold <- 10
    }else{
      if(size_threshold < 10 | is.na(size_threshold)){
        size_threshold <- 10
      }
    }
  }else{
    size_threshold = 0
  }

  peak_info <- peak_info %>% na.omit() %>% group_by(peak) %>%
      summarize(across(colnames(intensities), function(x){median(x, na.rm=TRUE) %>% return()}), n=n()) %>%
      na.omit() %>% mutate_if(is.numeric, round, digits=1) %>% filter(n > size_threshold) %>% select(!n)

  colnames(peak_info) <- paste(colnames(peak_info), "_median", sep="") %>% gsub("peak_median", "peak", .)

  # add position and height info
  peak_indices <- peak_info %>% pull(peak) %>% gsub("peak ", "", .) %>% as.numeric()
  peak_info <- peak_info %>% mutate(position=peaks$d.x[peak_indices]) %>% mutate(height=peaks$d.y[peak_indices])

  #summarize crosstalk from other channels
  peak_info <- peak_info %>% mutate(across(matches("median"), .fns = function(x){return(x-min(abs(x), na.rm=TRUE))})) %>%
    mutate(max_others = select(., !starts_with(variable)) %>%
             select(ends_with("median")) %>% apply(., 1, max))
  # rank intensity in current channel
  peak_info <- peak_info %>% mutate(median_current = peak_info %>% select(starts_with(variable)) %>% select(ends_with("median")) %>% unlist())

  # identify crosstalk peaks based on number of partitions positive in other channels
  if(recalc){
    # for first estimate of cross-talk, use simpler method
    ct_threshold <- peak_info %>%
      mutate(rank_others = rank(max_others)) %>%
      filter(rank_others < 3) %>%
      pull(max_others) %>%
      max()+ 0.1

    non_ct_peaks <- peak_info %>% filter(max_others < ct_threshold )

    if(nrow(non_ct_peaks)<2){
      non_ct_peaks <- peak_info %>% filter(max_others <= min(max_others) + 0.2)
    }

    peak_info <- peak_info %>% mutate(ct_peak = ifelse(.$peak %in% non_ct_peaks$peak, FALSE, TRUE))
  }else{
    # more precise method requires thresholds which are caluclated later
    ct_info <- data.frame(peak=numeric(length=0), channel_from=character(length=0), channel_to=character(length=0))

    all_partitions <- mutate(intensities, peak=cut(intensities %>% pull(variable) %>% as.numeric(),breaks=pits$d.x, labels=paste("peak", c(1:nrow(peaks))))) %>% na.omit()

    for(comp_channel in 1:length(channels)){
      if(channels[[comp_channel]] == variable | channels[[comp_channel]] %in% (coupled_channels %>% filter(ch1 == variable) %>% unlist() %>% unique())){
        next
      }else{
        expected_positive_ratio <- cut(intensities %>% pull(channels[[comp_channel]]), breaks= c(min(intensities %>% pull(channels[[comp_channel]]), na.rm=TRUE), thresholds %>% filter(channel == channels[[comp_channel]]) %>% pull(step) %>% as.numeric(), max(intensities %>% pull(channels[[comp_channel]]), na.rm=TRUE) ), labels=c(0, 1)) %>% na.omit() %>% table()
        expected_positive_ratio <- expected_positive_ratio["1"]/expected_positive_ratio["0"]

        peak_list_comp <- data.frame(peak=numeric(length=0), pos=numeric(length=0), neg=numeric(length=0))
        for(curr_peak in 1:nrow(peak_info)){
          peak_partitions <- all_partitions %>% filter(peak == paste("peak", curr_peak))

          n_pos <- peak_partitions %>% filter(peak==paste("peak", curr_peak)) %>% na.omit() %>% filter(!!sym(channels[[comp_channel]]) > thresholds %>% filter(channel==channels[[comp_channel]]) %>% pull(step)) %>% nrow()
          n_neg <- peak_partitions %>% filter(peak==paste("peak", curr_peak)) %>% na.omit() %>% filter(!!sym(channels[[comp_channel]]) <= thresholds %>% filter(channel==channels[[comp_channel]]) %>% pull(step)) %>% nrow()

          peak_list_comp <- full_join(peak_list_comp, data.frame(peak=curr_peak, pos=n_pos, neg=n_neg))
        }

        #peak_list_comp <- peak_list_comp %>% mutate(pos = replace(pos, pos==0 | is.na(pos), 1)) %>% mutate(neg = replace(neg, neg==0 | is.na(neg), 1))
        peak_list_comp <- peak_list_comp %>% mutate(rel=pos/neg)

        if(nrow(peak_list_comp)>1){
          ct_peaks <- peak_list_comp %>% filter(rel > 2*(expected_positive_ratio))
        }else{
          ct_peaks <- peak_list_comp[-1,]
        }

        if(nrow(ct_peaks)>0){
          ct_info_add <- data.frame(peak=ct_peaks %>% pull(peak), channel_from = channels[[comp_channel]], channel_to = variable)
          ct_info <- ct_info_add %>% full_join(., ct_info)

        }
      }
    }
    peak_info <- peak_info %>% mutate(ct_peak = ifelse(.$peak %in% paste("peak", ct_info %>% pull(peak)), TRUE, FALSE))
  }

  join_peaks = TRUE
  while(join_peaks){
    join_peaks = FALSE
    pits_temp <- pits
    for(pit in nrow(pits_temp):1){
      associated_peak_high <- peak_info %>% filter(position > pits_temp$d.x[[pit]]) %>% filter(position == min(position))
      associated_peak_low <- peak_info %>% filter(position < pits_temp$d.x[[pit]]) %>% filter(position == max(position))

      join_assoc_peaks <-(c(associated_peak_high$height, associated_peak_low$height) %>% max(., na.rm=TRUE))<2*pits_temp$d.y[[pit]] & !associated_peak_high$ct_peak & !associated_peak_low$ct_peak

      if(length(join_assoc_peaks)>0){
        if(!is.na(join_assoc_peaks)){
        if(join_assoc_peaks){
          joined_peak_means = associated_peak_high %>%
            full_join(associated_peak_low) %>%
            summarize_at(.vars = c(colnames(peak_info)[grep("_median", colnames(peak_info))], "max_others", "median_current"), .funs = mean)
          joined_peak_maxs = associated_peak_high %>%
            full_join(associated_peak_low) %>%
            summarize_at(.vars = c("height"), .funs = max)
          joined_peak <- cbind(joined_peak_maxs, joined_peak_means)
          joined_peak$peak <- paste(associated_peak_high$peak, associated_peak_low$peak)
          joined_peak$position <- associated_peak_high$position*(associated_peak_high$height/(associated_peak_low$height + associated_peak_high$height)) + associated_peak_low$position*(associated_peak_low$height/(associated_peak_low$height + associated_peak_high$height))
          joined_peak$ct_peak = FALSE
          pits_temp <- pits_temp[-pit,]
          peaks <- peaks %>%
            filter(d.x != associated_peak_low %>% pull(position)) %>%
            filter(d.x != associated_peak_high %>% pull(position))

          peaks <- peaks %>%
            full_join(joined_peak %>% select(height, position) %>% rename(d.y = height, d.x = position) %>% mutate(tp.peaks = TRUE, tp.pits = FALSE)) %>%
            arrange(d.x)

          peak_info <- peak_info %>%
            filter(!peak == associated_peak_high$peak & !peak ==associated_peak_low$peak) %>%
            full_join(., joined_peak)
          join_peaks = TRUE
        }
        }
      }
    }
    pits <- pits_temp
  }

  # define upper and lower limits as estimates for positive and negative band positions
  if(!is.null(reference_peaks)){
    upper_lim <- reference_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x)
    lower_lim <- reference_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x)
  }else{
    upper_lim <-intensities %>% pull(variable) %>% get_quantile(., probs=0.9999)
    lower_lim <-intensities %>% pull(variable) %>% get_quantile(., probs=0.0001)
  }

  # calculate peak distances to positive and negative bands
  peak_info <- peak_info %>% mutate(dist_up = abs(position - upper_lim)) %>% mutate(dist_low = abs(position - lower_lim))

  # consider only bands that are close to the lower limit and are not crosstalk peaks for the true negative band
  tn_candidates <- peak_info %>% filter(!ct_peak) %>% filter(dist_up > dist_low)

  # if no band fulfills conditions,  pick the a band that has as low fluorescence values as possible while being as large as possible
  if(is.null(reference_peaks)){
    if(nrow(tn_candidates)==0){
      tn_candidates <- peak_info %>% filter(!ct_peak) %>% mutate(rank=rank(-dist_low)*2 + rank(-max_others) + rank(height)) %>% filter(rank == max(rank))
      if(nrow(tn_candidates)==0){
        tn_candidates <- peak_info %>% mutate(rank=rank(-dist_low)*2 + rank(-max_others) + rank(height) + rank(-median_current)) %>% filter(rank == max(rank))
      }
    }
  }

  # if multiple bands are still candidates for tn: pick lowest of the candidates
  if(nrow(tn_candidates)>1){
    tn_peak <- tn_candidates %>%
      mutate(rank_height = rank(height)) %>% mutate(rank=rank(-dist_low)*2 + rank_height) %>% filter(rank == max(rank)) %>% filter(dist_low == min(dist_low))
  }else{
    tn_peak <- tn_candidates
  }

  # for the true positive band: consider only bands that are not crosstalk peaks and have a high fluorescence intensity
  tp_candidates <- peak_info %>% filter(!peak %in% tn_peak$peak) %>% filter(!ct_peak) %>% filter(dist_up < dist_low) %>%
    mutate(max_others_rank = rank(max_others)) %>% mutate(height_rank = rank(-height)) %>% mutate(dist_up_rank = rank(dist_up)) %>% mutate(rank = max_others_rank + 2*dist_up_rank + 0.5*height_rank) %>%
    filter(rank == min(rank, na.rm=TRUE))

  if(nrow(tp_candidates)>1){
    tp_peak <- tp_candidates %>%
      mutate(rank_height = rank(height)) %>% mutate(rank=rank(-dist_up)*2 + rank_height) %>% filter(rank == max(rank))
  }else{
    tp_peak <- tp_candidates
  }

  #if the same peak is classified as both TN and TP, decide based on distance to reference peaks:
  if(nrow(tn_peak)==1 & nrow(tp_peak)==1){
    if(tn_peak$peak[[1]] == tp_peak$peak[[1]]){
      if(tn_peak$dist_low > tp_peak$dist_up){
        tn_peak <- tn_peak[-1,]
      }else{
        tp_peak <- tp_peak[-1,]
      }
    }
  }

  tn_peak <- tn_peak %>% select(height, position) %>% rename(d.y = height, d.x = position) %>% mutate(tp.peaks = TRUE, tp.pits = FALSE)
  tp_peak <- tp_peak %>% select(height, position) %>% rename(d.y = height, d.x = position) %>% mutate(tp.peaks = TRUE, tp.pits = FALSE)

  if(nrow(tn_peak)==1 & nrow(tp_peak)==1){
    if(tn_peak$d.x==tp_peak$d.x){
      tn_peak <- peaks[tn_candidates %>% filter(position==min(position, na.rm=TRUE)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
      tp_peak <- peaks[tp_candidates %>% filter(position==max(position, na.rm=TRUE)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
    }
  }

  peaks <- peaks %>%
    mutate(peak=1:nrow(peaks)) %>%
    full_join(peak_info %>% rename(d.x = "position") %>%
    select(d.x, ct_peak)) %>%
    select(!peak) %>%
    filter(!is.na(d.x))

  if(nrow(tp_peak)>0){
    peaks <- peaks %>% mutate(tp.peak = ifelse(d.x==tp_peak$d.x, TRUE, FALSE))
  }else{
    peaks <- peaks %>% mutate(tp.peak = FALSE)
  }

  if(nrow(tn_peak)>0){
    peaks <- peaks %>% mutate(tn.peak = ifelse(d.x==tn_peak$d.x, TRUE, FALSE))
  }else{
    peaks <- peaks %>% mutate(tn.peak = FALSE)
  }

  tp_df_mod <- full_join(peaks, pits %>% mutate(tp.peak=FALSE) %>% mutate(tn.peak=FALSE) %>% mutate(ct_peak=FALSE)) %>% arrange(d.x)

  return(tp_df_mod)
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
