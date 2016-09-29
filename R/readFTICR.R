#' Read in chuncks of an FT-ICR-MS file
#'
#' @param fileIn a string identifying the csv file to read our FT-ICR-MS data from
#' @param massHeader an array of strings or single string naming the mass to charge header
#' @param sampleRegStr a regular expression matching the sample headers (and ONLY the sample headers)
#' @param samplesToRead an array of indecies or regular expression or array of names of the samples to be read in
#' @param verbose boolean flag for verbose outputs
#'
#' @return a data frame with the long table format with the  mass ID, sample ID, and intensity
#' @examples
#' ##read samples by index (note this is the Nth sample not the Nth column)
#' #readFTICR(fileIn='testfile.csv', samplesToRead=1:2)
#' ##read samples by regular expression match of sample name
#' #readFTICR(fileIn='testFile.csv', samplesToRead='out2\\d')
#' ##read samples by matched name
#' #readFTICR(fileIn='testFile.csv', samplesToRead=c('X.out20', 'X.out21'))
#' @import reshape2 assertthat
#'
#' @export
readFTICR <-  function(fileIn, massHeader, sampleRegStr=NULL, samplesToRead = NULL, elementKey=NULL,
                       verbose=FALSE){
  
  assert_that(file.exists(fileIn))
  
  if(verbose) cat('reading in from file [', fileIn, ']\n')
  header.df <- read.csv(fileIn, nrows=1)
  if(verbose) cat('file has [', dim(header.df)[2], '] columns\n')
  
  #Flag the mass header to be read in
  if(is.null(elementKey)){
    elementHeader <- c()
  }else{
    elementHeader <- names(elementKey)
  }
  colToReadIn <- names(header.df) %in% c(massHeader, elementHeader)
  massColCount <- sum(colToReadIn)
  
  ##Read in samples if specified
  readInCount <- 0
  if(!is.null(sampleRegStr)){
    if(is.numeric(samplesToRead)){
      sampleCols <- grepl(sampleRegStr, names(header.df))
      readInCount <- sum(sampleCols)
      if(any(samplesToRead > readInCount)){
        warning('Requesting sample index past maximum sample count. Truncating sample indecies.')
        samplesToRead <- samplesToRead[samplesToRead <= readInCount]
      }
      
      sampleCols[sampleCols][-1*samplesToRead] <- FALSE
      #Flag the N columns to read in 
      colToReadIn[sampleCols] <- TRUE
    }else if(is.character(samplesToRead)){
      if(length(samplesToRead) == 1){
        colToReadIn <- colToReadIn | grepl(samplesToRead, names(header.df))    
        readInCount <- sum(grepl(samplesToRead, names(header.df)))
      }else{
        colToReadIn <- colToReadIn | names(header.df) %in% samplesToRead
        readInCount <- sum(names(header.df) %in% samplesToRead)
      }
    }else{
      stop('Badly defined samples: ', samplesToRead)
    }
  }
  
  colReadArr <- as.character(colToReadIn)
  colReadArr[colToReadIn] <- 'numeric'
  colReadArr[!colToReadIn] <- 'NULL'
  
  if(verbose) cat('reading in ', massColCount, ' id columns and ', readInCount,'samples [total:', sum(colToReadIn),']\n')
  data.df <- read.csv(fileIn, colClasses=colReadArr)
  if(verbose){cat('orginal column names\n'); print(head(data.df))}
  data.df <- rename(data.df, elementKey)

  if(verbose){cat('renaming colums according to element key\n'); print(head(data.df))}
  
  if(!is.null(sampleRegStr)){
    if(verbose) cat('making samples long table\n')
    data.df <- melt(data.df, id.vars=c(massHeader[massHeader %in% names(data.df)], elementKey), variable.name='sample', value.name='intensity')
    
    if(verbose) cat('removing ', sum(data.df$intensity <= 0), ' with no intensity signal\n')
    #removing instensities at or below zero
    #..this is where the 0 mass's are removed for each sample for final count
    data.df<-data.df[data.df$intensity > 0,]
    sampleIDs <- c(massHeader[massHeader %in% names(data.df)], 
                   'sample', 'intensity')
  }else{
    sampleIDs <- massHeader[massHeader %in% names(data.df)]
  }
  
  
  if(!is.null(elementKey)){
    if(verbose)cat('make ratios and indexes for compounds\n')
    data.df <- normalizeMassNotation(data.df, sampleIDs=sampleIDs, elementKey=elementKey, verbose=verbose)
  }
  
  return(data.df)
}