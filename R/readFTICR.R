#' Read in chuncks of an FT-ICR-MS file
#'
#' @param fileIn a string identifying the csv file to read our FT-ICR-MS data from
#' @param massHeader an array of strings or single string naming the mass to charge header
#' @param elementHeader an array of strings or single string naming the element header
#' @param sampleRegStr a regular expression matching the sample headers (and ONLY the sample headers)
#' @param samplesToRead an array of indecies or regular expression or array of names of the samples to be read in
#' @param loadMassCharacteristics a boolean flag to append the mass characteristics to the samples
#' @param verbose boolean flag for verbose outputs
#'
#' @return a data frame with the long table format with the sample ID, mass ID, and intensity
#' @export
#' @import reshape2 plyr assertthat
#'
#' @examples
readFTICR <-  function(fileIn='data/smallTest.csv',
                       massHeader = c('Mass','m.z'),
                       elementHeader = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13'),
                       sampleRegStr = '(X.out)|(^X\\d+$)|(std)|(IntCal_)',
                       samplesToRead=1:2, loadMassCharacteristics=FALSE,
                       verbose=TRUE){
  assert_that(file.exists(fileIn))
  
  if(verbose) cat('reading in from file [', fileIn, ']\n')
  header.df <- read.csv(fileIn, nrows=1)
  
  ###Create mass to charge ratio table for all masses
  if(loadMassCharacteristics){
    massCharacteristic <- readMass(fileIn, massHeader, elementHeader, verbose)
  }
  
  ###Create the mass counts for each sample
  compounds.df <- data.frame()
  #Flag the mass header to be read in
  colToReadIn <- names(header.df) %in% massHeader
  if(is.numeric(samplesToRead)){
    sampleCols <- grepl(sampleRegStr, names(header.df))
    #Flag the N columns to read in
    colToReadIn[sampleCols][samplesToRead] <- TRUE
  }else if(is.character(samplesToRead)){
    if(length(samplesToRead) == 1){
      colToReadIn <- colToReadIn | grepl(samplesToRead, names(header.df))    
    }else{
      colToReadIn <- colToReadIn | names(header.df) %in% samplesToRead
    }
  }else{
    stop('Badly defined samples: ', samplesToRead)
  }
  
  colReadArr <- as.character(colToReadIn)
  colReadArr[colToReadIn] <- 'numeric'
  colReadArr[!colToReadIn] <- 'NULL'
  
  if(verbose) cat('reading in ', length(colReadArr), ' samples\n')
  data.df <- read.csv(fileIn, colClasses=colReadArr)
  if(verbose) print(head(data.df))
  
  if(verbose) cat('making samples long table\n')
  data.df <- melt(data.df, id.vars=massHeader[massHeader %in% names(data.df)], variable.name='sample', value.name='intensity')
  
  if(verbose) cat('removing ', sum(data.df$intensity > 0), ' with 0 intensity\n')
  #removing instensities at or below zero
  #..this is where the 0 mass's are removed for each sample for final count
  data.df<-data.df[data.df$intensity > 0,]
  
  if(loadMassCharacteristics){
    if(verbose) cat('merging with mass characteristics\n')
    #combined tables to list m/z for specific samples
    data.df <- merge(massCharacteristic, data.df)
  }
  return(data.df)
}