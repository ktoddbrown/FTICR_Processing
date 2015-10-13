#' Convert FT-ICR data to counts of compound classes
#' 
#' @decription
#' To convert FT-ICR intensity to compound classes, this function 1) discards masses that are of 0 intensity
#' 2) counts the number of masses for each sample which fall between specified C:H:O ratios which define the following
#' compound classes: (lipids, unsaturated hydrocarbons, proteins, lignins, carbohydrates, amino sugars, tannins,
#' condensed hydrocarbons, total number of masses which fell into those compounds, total number of masses 
#' which did not) as well as more basic grouping: (CHO, CHOS, CHOP, CHONS, CHOSP, CHONSP) and 
#' average O:C and H:C ratio.
#' 
#' @param fileIn A string specifiying the input csv file with header names corresponding to massHeader, and ratioHeaders.
#' @param fileOut A string specifying the output csv file.
#' @param massHeader A string or array of strings specifying the name of the colum with the mass value (normally first colum). If
#' this is an array then only one of the names should be valid for a given file (ie there can not be two columns named 'Mass' and 'm.z' if the defaults are not changed).
#' @param ratioHeader A string or array of strings specifying the columns with the elemental counts associated with
#' each mass. There can be extra strings in this array which do not match column headers, this will be ignored.
#' @param sampleRegStr A string that specifies the regular expression used to define sample columns. Note
#' that R will frequently stick an 'X' infront of headers that are strictly digits.
#' @param maxColReads A integer specifying the maximum number of columns to be read in at one time.
#' @param verbose A boolean flag to print out useful flags during processing.
#' 
#' @output  a data.frame with the count table (this is also saved as a csv to the 
#' specified output file)
#' 
#' @date October 2015
#' @author K Todd-Brown \email{ktoddbrown@gmail.com}
#' @author A P Smith
#' @author M Tfaily
#' 
#' @import reshape2 plyr
#' 

library(reshape2)
library(plyr)

countCompoundTypes <- function(fileIn='countCompoundTypesInput.csv', fileOut='countCompoundTypesOutput.csv', 
                               massHeader = c('Mass','m.z'),
                               ratioHeaders = c('C', 'H','O', 'N', 'X13C', 'S', 'P', 'C13', 'Na'),
                               sampleRegStr = '(X.out)|(^X\\d+$)|(std)',
                               maxColReads = 10, verbose=TRUE){
  
  header.df <- read.csv(fileIn, nrows=1)
  
  ###Create mass to charge ratio table for all masses
  massCols <- as.character(names(header.df) %in% c(massHeader, ratioHeaders))
  massCols[names(header.df) %in% c(massHeader, ratioHeaders)] <- 'numeric'
  massCols[!names(header.df) %in% c(massHeader, ratioHeaders)] <- 'NULL'
  massChargeRatio <- read.csv(fileIn, colClasses=massCols)
  massChargeRatio$OtoC <- massChargeRatio$O/massChargeRatio$C
  massChargeRatio$HtoC <- massChargeRatio$H/massChargeRatio$C
  
  ###Create the mass counts for each sample
  compounds.df <- data.frame()
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
    data.df <- merge(massChargeRatio, data.df)
    
    #assign counts to each sample based on the ratios associated with each mass
    temp.df <- ddply(data.df, c('sample'), function(xx){
      data.frame (
        Lipids = sum(xx$HtoC >= 1.55 & xx$OtoC <= 0.3, na.rm = TRUE),
        UnSaturated_Hydrocarbons = sum(xx$OtoC > 0.05 & xx$OtoC < 0.15 & 
                                         xx$HtoC < 1.5 & xx$HtoC > 0.7, na.rm = TRUE) ,
        Proteins = sum(xx$OtoC > 0.3 & xx$OtoC < 0.55 & 
                         xx$HtoC > 1.45 & xx$HtoC < 2, na.rm = TRUE) ,
        Lignin = sum(xx$OtoC > 0.28 & xx$OtoC < 0.65 & 
                       xx$HtoC < 1.45 & xx$HtoC > 0.81, na.rm = TRUE) , 
        Carbohydrates = sum(xx$OtoC < 1 & xx$OtoC > 0.68 & 
                              xx$HtoC > 1.48 & xx$HtoC < 2.15, na.rm = TRUE) , 
        Amino_Sugars = sum(xx$OtoC < .71 & xx$OtoC > 0.54 & 
                             xx$HtoC > 1.34 & xx$HtoC < 1.8, na.rm = TRUE) , 
        Tannins = sum(xx$OtoC < 1.05 & xx$OtoC > 0.65 & 
                        xx$HtoC > 0.7 & xx$HtoC < 1.3, na.rm = TRUE) , 
        Condensed_Hydrocarbons = sum(xx$OtoC < .7 & xx$OtoC > 0.12 & 
                                       xx$HtoC < .81 & xx$HtoC > .3, na.rm = TRUE) , 
        Count = sum(xx$C >= 0, na.rm = TRUE),
        CHO = sum(xx$C > 0 & xx$N == 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE ), 
        CHON = sum(xx$C > 0 & xx$N > 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE ),
        CHOS = sum(xx$C > 0 & xx$N == 0 & xx$S > 0 & xx$P == 0, na.rm = TRUE ),
        CHOP = sum(xx$C > 0 & xx$N == 0 & xx$S == 0 & xx$P > 0, na.rm = TRUE ),
        CHONS = sum(xx$C > 0 & xx$N > 0 & xx$S > 0 & xx$P == 0, na.rm = TRUE ),
        CHONP = sum(xx$C > 0 & xx$N > 0 & xx$S == 0 & xx$P > 0, na.rm = TRUE ),
        CHOSP = sum(xx$C > 0 & xx$N == 0 & xx$S > 0 & xx$P > 0, na.rm = TRUE ),
        CHONSP = sum(xx$C > 0 & xx$N > 0 & xx$S > 0 & xx$P > 0, na.rm = TRUE ),
        Not_Assigned =sum(xx$C == 0 & xx$N == 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE ),
        aveOtoC = mean(xx$OtoC, na.rm = TRUE),
        aveHtoC = mean(xx$HtoC, na.rm = TRUE)
      )
    })
    #New column called "Other" other is the Count - Sum of all compounds
    temp.df$Other <- temp.df$Count-rowSums(temp.df[,c('Lipids', 'UnSaturated_Hydrocarbons', 'Proteins', 'Lignin', 'Carbohydrates', 'Amino_Sugars', 'Tannins', 'Condensed_Hydrocarbons')])
    
    compounds.df <- rbind.fill(compounds.df, temp.df)
  }
  
  compounds.df <- compounds.df[,c(1:10, 22, 11:21)]
  write.csv(file=fileOut, compounds.df)
  #return(list(counts=compounds.df, massInfo=massChargeRatio)) #massChargeRatio is 3.2 Mb
  return(compounds.df)
}

