## ----librarys------------------------------------------------------------
library(FTICRProcessing)
library(ggplot2)
library(knitr)

## ----echo=FALSE----------------------------------------------------------
inputFile <- '../private/data/FTICRoutput.csv'

## ----eval=FALSE----------------------------------------------------------
#  inputFile <- 'data/ProcessedOutput.csv'

## ----loadcsv-------------------------------------------------------------
fileFormat <- read.csv(inputFile)[,c('m.z',  'X.out1.',  'X.out10.',  'X.out11.', 'X.out12.', 'C', 'H',  'O', 'N', 'X13C', 'S', 'P')]
kable(head(fileFormat[fileFormat[, 'X.out1.'] > 0,]))

## ----vanKrevlen, warning=FALSE, fig.width=4, fig.height=4----------------
mc <- readMass(inputFile, massHeader = 'm.z')

for(setNum in 1){
  msData <- readFTICR(inputFile, massHeader='m.z', sampleRegStr='X.out', samplesToRead=1:9+setNum*9)
  print(ggplot(merge(msData, mc)) + 
          geom_point(aes(x=OtoC, y=HtoC, color=log(intensity)), alpha=0.5) + 
          facet_wrap(~sample))
  print(ggplot(merge(msData, mc)) + 
          geom_point(aes(x=OtoC, y=HtoC, color=(S>0 | N > 0)), alpha=0.5) + 
          facet_wrap(~sample))
}

## ----oxyHist, warning=FALSE, fig.width=4, fig.height=4-------------------
mc <- readMass(inputFile, massHeader = 'm.z')

for(setNum in 1){
  msData <- readFTICR(inputFile, massHeader='m.z', sampleRegStr='X.out', samplesToRead=1:9+setNum*9)
  data.df <- merge(msData, mc)
  print(ggplot(data.df) + 
          geom_histogram(aes(x=O, y=..density..)) + xlim(1, max(data.df$O, na.rm=TRUE)) +
          facet_wrap(~sample))
  print(ggplot(data.df) + 
          geom_histogram(aes(x=O)) + xlim(1, max(data.df$O, na.rm=TRUE)) + 
          facet_wrap(~sample))
}

## ----Kendrick, warning=FALSE, fig.width=4, fig.height=4------------------
CO2mass <- 44.01
CH4mass <- 16.04
for(setNum in 1){
  msData <- readFTICR(inputFile, massHeader='m.z', sampleRegStr='X.out', samplesToRead=1:9+setNum*9)
  
  print(ggplot(msData) +
          geom_point(aes(x=m.z, y=m.z/CO2mass)))
}

