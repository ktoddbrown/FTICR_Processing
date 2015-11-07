library(ggplot2)

makeVanKrevelen <- function(fileIn='FTICRInputFile.csv',
                            massHeader = c('Mass','m.z'),
                               ratioHeaders = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13'),
                               sampleRegStr = '(X.out)|(^X\\d+$)|(std)|(IntCal_)',
                               samplesToPlot=1:9, colorBy=c('SN', 'intensity')[2],
                               verbose=TRUE){
  
  header.df <- read.csv(fileIn, nrows=1)
  
  ###Create mass to charge ratio table for all masses
  massCols <- as.character(names(header.df) %in% c(massHeader, ratioHeaders))
  massCols[names(header.df) %in% c(massHeader, ratioHeaders)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, ratioHeaders)] <- 'NULL'
  massCharacteristic <- read.csv(fileIn, colClasses=massCols)
  massCharacteristic$OtoC <- massCharacteristic$O/massCharacteristic$C
  massCharacteristic$HtoC <- massCharacteristic$H/massCharacteristic$C
  massCharacteristic$AImod <- (1 + massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$H*0.5)/(massCharacteristic$C - massCharacteristic$O*0.5 - massCharacteristic$S - massCharacteristic$N - massCharacteristic$P)
  
  ###Create the mass counts for each sample
  compounds.df <- data.frame()
  sampleCols <- grepl(sampleRegStr, names(header.df))
  #Flag the mass header to be read in
  colToReadIn <- names(header.df) %in% massHeader
  #Flag the N columns to read in
  colToReadIn[sampleCols][samplesToPlot] <- TRUE
  
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
  
  if(colorBy %in% 'intensity'){
  ans <- ggplot(data.df) + geom_point(aes(x=HtoC, y=OtoC, color=log(intensity))) + facet_wrap(~sample)
  }
  if(colorBy %in% 'SN'){
    ans <- ggplot(data.df) + geom_point(aes(x=HtoC, y=OtoC, color=(S>0 | N > 0)), alpha=0.5)+ facet_wrap(~sample)
  }
  
  return(ans)
}