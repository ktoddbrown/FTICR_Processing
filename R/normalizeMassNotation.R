#' Standardize mass notation
#' 
#' Force C, H, O, N, S, and P for the headers and calculate various ratios and indecies of interest.
#'
#' @param data.df data frame with samples and elements
#' @param sampleIDs an array of strings with the sample column names
#' @param elementKey a list with names that match C, H, O, N, S, P with the associated column names. Note that if you want to count C+C13 then set list(C=c('C', 'C13')).
#' @param verbose a boolean to flag lots of output that is hopefully informative
#'
#' @return data.frame with indecies and ratios of interest plus the standardize element names and sampleIDs
#' @export
normalizeMassNotation <- function(data.df, sampleIDs,
                                  elementKey = list(C='C', H='H', O='O', N='N', S='S', P='P'),
                                  verbose=FALSE){
  assert_that(class(data.df) %in% 'data.frame')
  ##Force expected naming conventions
#   for(elementStr in unlist(elementKey)){
#     if(sum(names(data.df) %in% elementKey[[elementStr]]) == 0){
#       stop('column name in elementKey does not match csv table: ', elementKey[elementStr])
#     }
#     if(sum(names(data.df) %in% elementKey[[elementStr]]) > 1){
#       if(verbose) cat('merging', elementStr,'counts; ')
#       data.df[[elementStr]] <- rowSums(data.df[,names(data.df) %in% elementKey[[elementStr]]])
#     }else{
#       if(verbose) cat('forcing',elementStr,'header name; ')
#       data.df[[elementStr]] <- data.df[,names(data.df) %in% elementKey[[elementStr]]]
#     }
#   }
  if(verbose) cat('\n')
  if(verbose) print(head(data.df))
  
  ##Only return expected convention
  data.df <- data.df[,c(sampleIDs, unlist(elementKey))]
  
  ###Make metrics
  if(all(c('C', 'H', 'O') %in% names(data.df))){
    if(verbose) cat('calculating: OtoC, HtoC\n')
    data.df$OtoC <- data.df$O/data.df$C
    data.df$HtoC <- data.df$H/data.df$C
    
    if(all(c('N', 'S', 'P') %in% names(data.df))){
      if(verbose) cat('calculating: DBE, AI, AImod\n')
      #Koch, B. P. and Dittmar, T.: From mass to structure: an aromaticity index for high-resolution mass data of natural organic matter, Rapid Commun. Mass Spectrom., 20(5), 926–932, doi:10.1002/rcm.2386, 2006.
      data.df$DBE <- 1+ 0.5*(2*data.df$C - data.df$H + data.df$N + data.df$P)
      data.df$AI <- (1 + data.df$C - data.df$O - data.df$S - data.df$H*0.5)/(data.df$C - data.df$O - data.df$S - data.df$N - data.df$P)
      data.df$AImod <- (1 + data.df$C - data.df$O*0.5 - data.df$S - data.df$H*0.5)/(data.df$C - data.df$O*0.5 - data.df$S - data.df$N - data.df$P)
      if(verbose) print(head(data.df))
    }else{
      warning('N, S, and P not found in elementKey. Not calculating DBE, AI, or AImod.')
    }
  }else{
    warning('O, H, and C not found in elementKey. Not calculating OtoC, HtoC, DBE, AI, or AImod.')
  }
  
  return(data.df)
}