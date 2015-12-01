#' Read in characteristics of the mass to charge ratios
#' 
#'
#' @param fileIn a string identifying the csv file to read our data from
#' @param massHeader an array of strings or single string naming the mass to charge header
#' @param elementKey a list of all the column names to be counted for each element
#' @param verbose boolean flag for verbose outputs
#'
#' @return a data frame with the mass, orginal elements [if names different from standard], 
#' element count, and metrics [if appropreate elements are present]
#' @examples
#' ##massHeader can have extra entries that will be ignored
#' #readMass(fileIn='testfile.csv', massHeader=c('m.z', 'mass', 'nonsence'))
#' 
#' ##elementKey can contain isotopes which will be included in the final element count
#' #readMass(fileIn='testfile.csv', massHeader='m.z', elementKey=list(C=c('C13', 'C'),
#' #H='H', O='O', N=c('N', 'N15'), S='S', P='P'))
#' 
#' @import assertthat
#'
#' @export
readMass <-  function(fileIn, massHeader = 'Mass',
                       elementKey = list(C='C', H='H', O='O', N='N', S='S', P='P'),
                       verbose=FALSE){
  assert_that(file.exists(fileIn))
  
  if(verbose) cat('reading in from file [', fileIn, ']\n')
  header.df <- read.csv(fileIn, nrows=1)
  
  elementHeader <- unlist(elementKey)
  ###Create mass to charge ratio table for all masses
  
  if(verbose) cat('reading in mass characteristics for colms [',elementHeader,']\n')
  massCols <- as.character(names(header.df) %in% c(massHeader, elementHeader))
  massCols[names(header.df) %in% c(massHeader, elementHeader)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, elementHeader)] <- 'NULL'
  massCharacteristic <- read.csv(fileIn, colClasses=massCols)
  if(verbose) print(head(massCharacteristic))
  
  ##Force expected naming conventions
  for(elementStr in names(elementKey)){
    if(sum(names(massCharacteristic) %in% elementKey[[elementStr]]) == 0){
      stop('column name in elementKey does not match csv table: ', elementKey[elementStr])
    }
    if(sum(names(massCharacteristic) %in% elementKey[[elementStr]]) > 1){
      if(verbose) cat('merging', elementStr,'counts \n')
      massCharacteristic[[elementStr]] <- rowSums(massCharacteristic[,names(massCharacteristic) %in% elementKey[[elementStr]]])
    }else{
      if(verbose) cat('forcing',elementStr,'header name \n')
      massCharacteristic[[elementStr]] <- massCharacteristic[,names(massCharacteristic) %in% elementKey[[elementStr]]]
    }
  }

  if(verbose) print(head(massCharacteristic))
  
  ###Make metrics
  if(all(c('C', 'H', 'O') %in% names(elementKey))){
    if(verbose) cat('calculating: OtoC, HtoC\n')
    massCharacteristic$OtoC <- massCharacteristic$O/massCharacteristic$C
    massCharacteristic$HtoC <- massCharacteristic$H/massCharacteristic$C
    
    if(all(c('N', 'S', 'P') %in% names(elementKey))){
      if(verbose) cat('calculating: DBE, AI, AImod')
      #Koch, B. P. and Dittmar, T.: From mass to structure: an aromaticity index for high-resolution mass data of natural organic matter, Rapid Commun. Mass Spectrom., 20(5), 926–932, doi:10.1002/rcm.2386, 2006.
      massCharacteristic$DBE <- 1+ 0.5*(2*massCharacteristic$C - massCharacteristic$H + massCharacteristic$N + massCharacteristic$P)
      massCharacteristic$AI <- (1 + massCharacteristic$C - massCharacteristic$O - massCharacteristic$S - massCharacteristic$H*0.5)/(massCharacteristic$C - massCharacteristic$O - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
      massCharacteristic$AImod <- (1 + massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$H*0.5)/(massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
      if(verbose) print(head(massCharacteristic))
    }else{
      warning('N, S, and P not found in elementKey. Not calculating DBE, AI, or AImod.')
    }
  }else{
    warning('O, H, and C not found in elementKey. Not calculating OtoC, HtoC, DBE, AI, or AImod.')
  }
  
  return(massCharacteristic)
}