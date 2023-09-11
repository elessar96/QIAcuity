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