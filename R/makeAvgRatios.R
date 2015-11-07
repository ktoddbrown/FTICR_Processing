
library(reshape2)
library(plyr)
library(assertthat)

makeAvgRatios <- function(fileIn='FTICRInputFile.csv',
                        massHeader = c('Mass','m.z'),
                        ratioHeaders = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13'),
                        sampleRegStr = '(X.out)|(^X\\d+$)|(std)|(IntCal_)',
                        weighByIntensity=TRUE,
                        verbose=TRUE){
  data.df <- readFTICR(fileIn, massHeader, ratioHeaders, sampleRegStr, samplesToRead=samplesToRead, verbose)
  
  header.df <- read.csv(fileIn, nrows=1)
  
  ###Create mass to charge ratio table for all masses
  massCols <- as.character(names(header.df) %in% c(massHeader, ratioHeaders))
  massCols[names(header.df) %in% c(massHeader, ratioHeaders)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, ratioHeaders)] <- 'NULL'
  massCharacteristic <- read.csv(fileIn, colClasses=massCols)
  massCharacteristic$OtoC <- massCharacteristic$O/massCharacteristic$C
  massCharacteristic$HtoC <- massCharacteristic$H/massCharacteristic$C
  massCharacteristic$AImod <- (1 + massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$H*0.5)/(massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
  
  ans.df <- data.frame()
  sampleCols <- grepl(sampleRegStr, names(header.df))
  splitCol.ls <- split(1:sum(sampleCols), ceiling(seq_along(1:sum(sampleCols))/maxColReads))
  
  for(sampleIndex in splitCol.ls){
    if(verbose) cat(sprintf('processing sample: %d of %d ...\n', max(sampleIndex), sum(sampleCols)))
    
    #Flag the mass header to be read in
    colToReadIn <- names(header.df) %in% massHeader
    #Flag the N columns to read in
    colToReadIn[sampleCols][sampleIndex] <- TRUE
    
    colReadArr <- as.character(colToReadIn)
    colReadArr[colToReadIn] <- 'numeric'
    colReadArr[!colToReadIn] <- 'NULL'
    
    data.df <- read.csv(fileIn, colClasses=colReadArr)
    data.df <- melt(data.df, id.vars=massHeader[massHeader %in% names(data.df)], variable.name='sample', value.name='intensity')
    
    #removing instensities at or below zero
    #..this is where the 0 mass's are removed for each sample for final count
    data.df<-data.df[data.df$intensity > 0,]
    
    #combined tables to list m/z for specific samples
    data.df <- merge(massCharacteristic, data.df)
    
    if(weighByIntensity){
      temp.df <- ddply(data.df, c('sample'), summarize, OtoCMean=sum(OtoC*intensity, na.rm=TRUE)/sum(is.finite(OtoC)*intensity, na.rm=TRUE), HtoC=sum(HtoC*intensity, na.rm=TRUE)/sum(is.finite(HtoC)*intensity, na.rm=TRUE))
    }else{
      temp.df <- ddply(data.df, c('sample'), summarize, OtoCMean=mean(OtoC, na.rm=TRUE), HtoC=mean(HtoC, na.rm=TRUE))
    }
    ans.df <- rbind.fill(ans.df, temp.df)
  }
  return(ans.df)
}
  
