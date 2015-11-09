#' Read in characteristics of the mass to charge ratios
#'
#' @param fileIn a string identifying the csv file to read our data from
#' @param massHeader an array of strings or single string naming the mass to charge header
#' @param elementHeader an array of strings or single string naming the element header
#' @param verbose boolean flag for verbose outputs
#'
#' @return a data frame with the long table format with the sample ID, mass ID, and intensity
#' @export
#' @import assertthat
#'
#' @examples
readMass <-  function(fileIn='data/smallTest.csv',
                       massHeader = c('Mass','m.z'),
                       elementHeader = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13'),
                       verbose=FALSE){
  assert_that(file.exists(fileIn))
  
  if(verbose) cat('reading in from file [', fileIn, ']\n')
  header.df <- read.csv(fileIn, nrows=1)
  
  ###Create mass to charge ratio table for all masses
  
  if(verbose) cat('reading in mass characteristics\n')
  massCols <- as.character(names(header.df) %in% c(massHeader, elementHeader))
  massCols[names(header.df) %in% c(massHeader, elementHeader)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, elementHeader)] <- 'NULL'
  massCharacteristic <- read.csv(fileIn, colClasses=massCols)
  
  if(verbose) cat('calculating: OtoC, HtoC, AImod\n')
  massCharacteristic$OtoC <- massCharacteristic$O/massCharacteristic$C
  massCharacteristic$HtoC <- massCharacteristic$H/massCharacteristic$C
  massCharacteristic$AImod <- (1 + massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$H*0.5)/(massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
  if(verbose) print(head(massCharacteristic))
  
  return(massCharacteristic)
}