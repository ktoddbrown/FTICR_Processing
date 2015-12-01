#' Read in characteristics of the mass to charge ratios
#' 
#' We assume that the elemental headers in the file are: 'H', 'O', 'N', 'S', 'P'. Carbon headers can be passed as arguments 
#'
#' @param fileIn a string identifying the csv file to read our data from
#' @param massHeader an array of strings or single string naming the mass to charge header
#' @param carbonHeader an array of strings or single string naming the carbon headers
#' @param verbose boolean flag for verbose outputs
#'
#' @return a data frame with the long table format with the sample ID, mass ID, and intensity
#' @export
#' @import assertthat
#'
readMass <-  function(fileIn='data/smallTest.csv',
                       massHeader = c('Mass','m.z'),
                       carbonHeader = c('C', 'X13C', 'C13')[1],
                       verbose=FALSE){
  assert_that(file.exists(fileIn))
  
  if(verbose) cat('reading in from file [', fileIn, ']\n')
  header.df <- read.csv(fileIn, nrows=1)
  
  elementHeader <- c(carbonHeader, 'H','O', 'N', 'S', 'P')
  ###Create mass to charge ratio table for all masses
  
  if(verbose) cat('reading in mass characteristics\n')
  massCols <- as.character(names(header.df) %in% c(massHeader, elementHeader))
  massCols[names(header.df) %in% c(massHeader, elementHeader)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, elementHeader)] <- 'NULL'
  massCharacteristic <- read.csv(fileIn, colClasses=massCols)
  
  if(verbose) cat('calculating: OtoC, HtoC, AImod\n')
  massCharacteristic$OtoC <- massCharacteristic$O/sum(massCharacteristic[names(massCharacteristic) %in% carbonHeader])
  massCharacteristic$HtoC <- massCharacteristic$H/sum(massCharacteristic[names(massCharacteristic) %in% carbonHeader])
  
  #Koch, B. P. and Dittmar, T.: From mass to structure: an aromaticity index for high-resolution mass data of natural organic matter, Rapid Commun. Mass Spectrom., 20(5), 926–932, doi:10.1002/rcm.2386, 2006.
  massCharacteristic$DBE <- 1+ 0.5*(2*sum(massCharacteristic[names(massCharacteristic) %in%  carbonHeader]) - massCharacteristic$H + massCharacteristic$N + massCharacteristic$P)
  massCharacteristic$AI <- (1 + sum(massCharacteristic[names(massCharacteristic) %in%  carbonHeader]) - massCharacteristic$O - massCharacteristic$S - massCharacteristic$H*0.5)/(sum(massCharacteristic[names(massCharacteristic) %in%  carbonHeader]) - massCharacteristic$O - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
  massCharacteristic$AImod <- (1 + sum(massCharacteristic[names(massCharacteristic) %in%  carbonHeader]) - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$H*0.5)/(sum(massCharacteristic[names(massCharacteristic) %in%  carbonHeader]) - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
  if(verbose) print(head(massCharacteristic))
  
  return(massCharacteristic)
}