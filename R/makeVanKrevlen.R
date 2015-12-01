#' Create a Van Krevlen plot
#'
#' @param fileIn csv file with the FT-ICR-MS data
#' @param massHeader array of strings with the names of the mass to charge ratio column
#' @param ratioHeaders array of strings with the names of the elemental count columns
#' @param sampleRegStr regular expression defining the sample columns
#' @param samplesToPlot number of samples to plot
#' @param colorBy string to color plot by intensity or by inclusion of S or N in element assignment
#' @param verbose boolean flag for verbose outputs
#'
#' @import ggplot2 reshape2
#' @return a ggplot object with the graph
#' @export
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
  ans <- ggplot(data.df) + geom_point(aes(x=OtoC, y=HtoC, color=log(intensity))) + facet_wrap(~sample)
  }
  if(colorBy %in% 'SN'){
    ans <- ggplot(data.df) + geom_point(aes(x=OtoC, y=HtoC, color=(S>0 | N > 0)), alpha=0.5)+ facet_wrap(~sample)
  }
  
  return(ans)
}