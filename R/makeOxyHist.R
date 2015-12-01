
makeOxyHist <- function(fileIn='FTICRInputFile.csv',
                        massHeader = c('Mass','m.z'),
                        ratioHeaders = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13'),
                        sampleRegStr = '(X.out)|(^X\\d+$)|(std)|(IntCal_)',
                        samplesToPlot=1:9, density=TRUE,
                        verbose=TRUE){
  data.df <- readFTICR(fileIn, massHeader, ratioHeaders, sampleRegStr, samplesToRead=samplesToPlot, verbose)
  
  if(density){
  ggplot(data.df) + geom_histogram(aes(x=O, y=..density..)) + xlim(1, max(data.df$O, na.rm=TRUE))+ facet_wrap(~sample)
  }else{
    ggplot(data.df) + geom_histogram(aes(x=O)) + xlim(1, max(data.df$O, na.rm=TRUE))+ facet_wrap(~sample)
  }
}