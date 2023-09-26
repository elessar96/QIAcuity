
#' @title Define a threshold based on mixture modeling
#' @description This function defines a threshold that best divides a population of data points into two separate groups based on mixture modeling.
#' @param input_data An object of class dataframe containing fluorescence intensities fom ddPCR partitions.
#' @param variable  Name of the variable in the dataframe to be pulled
#' @param min_dist Determines the minimum distance between threshold and negative partitions, always as a fraction of the distance between true positive and true negative populations.
#' @param references Reference peaks for guiding classification
#' @param pc_data Positive control data as a reference to calculate maximum expected values for each channel
#' @return A numeric value giving the optimal threshold for distinguishing the true positive and true negative populations in the data.

density_threshold <- function(input_data, variable, min_dist=0.2, references=NULL, pc_data =NULL, mid=FALSE){

  data_subset <- input_data %>% select(!any_of(c("Well", "Sample", "Partition"))) %>% na.omit()

  # find turnpoints
  tp_df <- find_turnpoints(data_subset, variable=variable)

  if(is.null(pc_data)){
    pc_data <- data_subset
  }

  pc_data <- pc_data %>% select(!any_of(c("Well", "Sample", "Partition")))

  channel_maxima <- pc_data %>% apply(., 2, get_quantile, probs=0.999)

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

  if(mid==TRUE){
    threshold <- positions$tn + (positions$tp-positions$tn)*0.5
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
  peak_info <- peak_info %>% na.omit() %>% group_by(peak) %>%
    summarize_at(.vars = colnames(intensities), .funs= c("min"=min, "median"=median, "max"=max), na.rm=TRUE) %>%
    na.omit() %>% mutate_if(is.numeric, round, digits=1)

  if(ncol(peak_info)==4){
    peak_info <- peak_info %>% rename(., !!sym(paste(variable, "_median", sep="")):=median, !!sym(paste(variable, "_min", sep="")):=min, !!sym(paste(variable, "_max", sep="")):=max)
  }

  # add position and height info
  peak_indices <- peak_info %>% pull(peak) %>% gsub("peak ", "", .) %>% as.numeric()
  peak_info <- peak_info %>% mutate(position=peaks$d.x[peak_indices]) %>% mutate(height=peaks$d.y[peak_indices])

  #summarize crosstalk from other channels
  peak_info <- peak_info %>% mutate(max_others = peak_info %>% select(!starts_with(variable)) %>%
                                      select(ends_with("median")) %>% apply(., 1, max))

  # rank intensity in current channel
  peak_info <- peak_info %>% mutate(median_current = peak_info %>% select(starts_with(variable)) %>% select(ends_with("median")) %>% unlist())

  if(!is.null(reference_peaks)){
    # if reference TN and TP peaks are provided, make decision purely on peak positions!
    if(nrow(peak_info)>2){
      ref_height <- peak_info %>% mutate(rank_height = rank(-height)) %>% filter(rank_height ==2 ) %>% pull(height)
      peak_info <- peak_info %>% filter(height>(ref_height/5))
    }

    peak_info <- peak_info %>% mutate(dist.tn=abs(position - reference_peaks %>% filter(tn.peak==TRUE) %>% pull(d.x))) %>%
                      mutate(dist.tp=abs(position - reference_peaks %>% filter(tp.peak==TRUE) %>% pull(d.x)))

    peak_info <- peak_info %>% mutate(tn.peak.dist = ifelse(dist.tn == dist.tn %>% abs() %>% min(., na.rm=TRUE), 1, 0)) %>%
                      mutate(tp.peak = ifelse(dist.tp == dist.tp %>% abs() %>% min(., na.rm=TRUE), TRUE, FALSE))

    tn_candidates <- peak_info %>% filter(height > max(height)/5000) %>% filter(dist.tp > dist.tn)

    if(nrow(tn_candidates)>1){
      tn_candidates <- tn_candidates %>% filter(max_others < max_others %>% median()+0.2) %>%
        mutate(rank_height = rank(height)) %>% mutate(rank=tn.peak.dist*2 + rank_height) %>% filter(rank == max(rank))
    }

    tp_candidates <- peak_info %>% filter(dist.tp < dist.tn) %>% filter(tp.peak==TRUE) %>%
      filter(dist.tn > abs(reference_peaks[1,]$d.x - reference_peaks[2,]$d.x)/4) %>% arrange(peak) %>%
      mutate(rank_dist.tn =rank(-dist.tn)) %>% mutate(rank_dist.tp =rank(dist.tp)) %>% mutate(max_others = rank(max_others)) %>%
      mutate(rank = max_others + rank_dist.tp) %>% filter(rank == min(rank, na.rm=TRUE))

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
        tn_candidates <- tn_candidates %>% arrange(peak) %>% mutate(position = rank(abs(as.numeric(position)))) %>% mutate(rank = position) %>%
          filter(rank == min(rank, na.rm=TRUE))
      }else{
        ct_detection <- tn_candidates %>% mutate(across(.cols=matches("median"), .fns=function(x){return(x-min(x, na.rm=TRUE))}))
        ct_detection <- ct_detection %>% mutate(max_others = ct_detection %>% select(!starts_with(variable)) %>% select(matches("median")) %>%
                                                  apply(., MARGIN=1, FUN=max))

        non_ct_peaks <- ct_detection %>% filter(max_others < min(max_others, na.rm=TRUE)+0.5*sd(max_others))
        if(nrow(non_ct_peaks)==0){
          non_ct_peaks <- tn_candidates %>% filter(max_others <= min(max_others, na.rm = TRUE) + 0.2)
        }
        tn_candidates <- non_ct_peaks %>% arrange(peak)  %>% mutate(max_others = rank(max_others)) %>%
          mutate(position = rank(abs(as.numeric(position)))) %>% mutate(rank = max_others + position) %>% filter(rank == min(rank, na.rm=TRUE)) %>%
          filter(position==min(position, na.rm=TRUE))
      }
      tp_candidates <- peak_info %>% filter(height > height %>% max()/100) %>% filter(position > (peak_info %>% filter(peak %in% tn_candidates$peak) %>%
                              filter(height > height %>% max()/30)  %>% pull(position) %>% max())) %>% filter(!peak %in% c(tn_candidates %>%
                              filter(position==min(position, na.rm=TRUE)) %>% pull(peak))) %>% arrange(peak) %>% mutate(pos = peaks$d.x[.$peak]) %>%
                              mutate(max_others = rank(max_others)) %>% mutate(position = rank(-position)) %>% mutate(rank = max_others + position) %>%
                              filter(rank == min(rank, na.rm=TRUE))
    }else{
      tn_candidates <- peak_info %>% filter(position==min(position, na.rm=TRUE))
      tp_candidates <- peak_info %>% filter(position==max(position, na.rm=TRUE))
    }
  }

  tn_peak <- peaks[tn_candidates %>% filter(max_others==min(max_others, na.rm=TRUE)) %>% filter(position==min(position, na.rm=TRUE)) %>%
                     pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]

  tp_peak <- peaks[tp_candidates %>% filter(position==min(position, na.rm=TRUE)) %>% filter(max_others==min(max_others, na.rm=TRUE)) %>%
                     pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]

  if(nrow(tn_peak)==1 & nrow(tp_peak)==1){
    if(tn_peak$d.x==tp_peak$d.x){
      tn_peak <- peaks[tn_candidates %>% filter(position==min(position, na.rm=TRUE)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
      tp_peak <- peaks[tp_candidates %>% filter(position==max(position, na.rm=TRUE)) %>% pull(peak) %>% as.character() %>% gsub("peak ", "", .) %>% as.numeric(),]
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
    if(!tn_peak$d.y == peaks %>% filter(crosstalk.peak ==FALSE)  %>% pull(d.y) %>% max(., na.rm=TRUE)){
      tn_peak <- peaks %>% filter(crosstalk.peak ==FALSE) %>% filter(d.y==max(d.y, na.rm=TRUE))
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
