#' Convert FT-ICR data to counts of compound classes
#' 
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
#' @param sampleRegStr A string that specifies the regular expression used to define sample columns. Note
#' that R will frequently stick an 'X' infront of headers that are strictly digits.
#' @param maxColReads A integer specifying the maximum number of columns to be read in at one time.
#' @param verbose A boolean flag to print out useful flags during processing.
#' 
#' @return  a data.frame with the count table (this is also saved as a csv to the 
#' specified output file)
#' 
#' 
#' @import reshape2 plyr assertthat
#' @export
countCompoundTypes <- function(fileIn='countCompoundTypesInput.csv', fileOut='countCompoundTypesOutput.csv', 
                               massHeader = c('Mass','m.z'),
                               carbonHeader = c('C', 'X13C', 'C13')[1],
                               sampleRegStr = '(X.out)|(^X\\d+$)|(std)|(IntCal_)',
                               maxColReads = 10, 
                               calculateClass=c('compound', 'molecular', 'aromaticity'), 
                               verbose=FALSE){
  #@date October 2015
  #' @author K Todd-Brown \email{ktoddbrown@gmail.com}
  #' @author A P Smith
  #' @author M Tfaily
  header.df <- read.csv(fileIn, nrows=1)
  
  massCharacteristic <- readMass(fileIn=fileIn, massHeader=massHeader, carbonHeader=carbonHeader, verbose=verbose)
 
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
    data.df <- merge(massCharacteristic, data.df)
    
    categoryCuts <- function(xx, ranges=list(OtoC=c(0,1)), inclusive=list(OtoC=c(TRUE, TRUE)), 
                             na.rm=TRUE){
      assert_that(identical(names(ranges), names(inclusive) ))
      ans = 1
      for(nameStr in names(ranges)){
        ans <- ans *((ranges[[nameStr]][1] < xx[[nameStr]]) + 
                       inclusive[[nameStr]][1]*(ranges[[nameStr]][1]==xx[[nameStr]]))*
          ((ranges[[nameStr]][2] > xx[[nameStr]]) + 
             inclusive[[nameStr]][2]*(ranges[[nameStr]][2]==xx[[nameStr]]))
      }
      
      return(sum(ans > 0, na.rm=TRUE))
    }
    
    #assign counts to each sample based on the ratios associated with each mass
    temp.df <- ddply(data.df, c('sample'), function(xx){
      ans2 <- data.frame(total_peaks = sum(xx$C >= 0, na.rm = TRUE),
                         aveOtoC = mean(xx$OtoC, na.rm = TRUE),
                         aveHtoC = mean(xx$HtoC, na.rm = TRUE))
      
      if('molecular' %in% calculateClass){
        ans2 <- cbind(ans2, data.frame(
          CHO = sum(xx$C > 0 & xx$N == 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE ), 
          CHON = sum(xx$C > 0 & xx$N > 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE ),
          CHOS = sum(xx$C > 0 & xx$N == 0 & xx$S > 0 & xx$P == 0, na.rm = TRUE ),
          CHOP = sum(xx$C > 0 & xx$N == 0 & xx$S == 0 & xx$P > 0, na.rm = TRUE ),
          CHONS = sum(xx$C > 0 & xx$N > 0 & xx$S > 0 & xx$P == 0, na.rm = TRUE ),
          CHONP = sum(xx$C > 0 & xx$N > 0 & xx$S == 0 & xx$P > 0, na.rm = TRUE ),
          CHOSP = sum(xx$C > 0 & xx$N == 0 & xx$S > 0 & xx$P > 0, na.rm = TRUE ),
          CHONSP = sum(xx$C > 0 & xx$N > 0 & xx$S > 0 & xx$P > 0, na.rm = TRUE ),
          Molecular_NA =sum(xx$C == 0 & xx$N == 0 & xx$S == 0 & xx$P == 0, na.rm = TRUE )))
      }
      if('compound' %in% calculateClass){
        #November 2015 BOUNDARIES that we have been using: (M Tfaily)
        #class   O:C(low) O:C(high) H:C(low) H:C(high)
        #lipid     >0     0.3       1.5        2.5
        #unsatHC    0     0.125     0.8       <1.5
        #condHC     0     0.95      0.2       <0.8
        #protein   >0.3   0.55      1.5        2.3
        #aminosugar>0.55  0.7       1.5        2.2
        #carb      >0.7   1.5       1.5        2.5
        #lignin    >0.125 0.65      0.8       <1.5
        #tannin    >0.65  1.1       0.8       <1.5
        ans2 <- cbind(ans2, data.frame(
          Lipids = categoryCuts(xx, ranges=list(OtoC=c(0, 0.3), HtoC=c(1.5, 2.5)), 
                                inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE,TRUE))),
          UnSaturated_Hydrocarbons = categoryCuts(xx, ranges=list(OtoC=c(0, 0.125), HtoC=c(0.8, 1.5)), 
                                                  inclusive=list(OtoC=c(TRUE, TRUE), HtoC=c(TRUE,FALSE))),
          Condensed_Hydrocarbons   = categoryCuts(xx, ranges=list(OtoC=c(0, 0.95), HtoC=c(0.2, 0.8)), 
                                                  inclusive=list(OtoC=c(TRUE, TRUE), HtoC=c(TRUE,FALSE))), 
          Proteins = categoryCuts(xx, ranges=list(OtoC=c(0.3, 0.55), HtoC=c(1.5, 2.3)), 
                                  inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE,TRUE))) ,
          Amino_Sugars = categoryCuts(xx, ranges=list(OtoC=c(0.55, 0.7), HtoC=c(1.5,2.2)), 
                                      inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE, TRUE))) , 
          Carbohydrates = categoryCuts(xx, ranges=list(OtoC=c(0.7, 1.5), HtoC=c(1.5, 2.5)), 
                                       inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE,TRUE))) , 
          Lignin = categoryCuts(xx, ranges=list(OtoC=c(0.125, 0.65), HtoC=c(0.8, 1.5)), 
                                inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE, FALSE))) , 
          Tannins = categoryCuts(xx, ranges=list(OtoC=c(0.65, 1.1), HtoC=c(0.8, 1.5)), 
                                 inclusive=list(OtoC=c(FALSE, TRUE), HtoC=c(TRUE, FALSE)))
        ))
        ans2$Compounds_NA <- ans2$total_peaks-rowSums(ans2[,c('Lipids', 'UnSaturated_Hydrocarbons', 'Proteins', 'Lignin', 'Carbohydrates', 'Amino_Sugars', 'Tannins', 'Condensed_Hydrocarbons')])
      }
      if('aromaticity' %in% calculateClass){
        ans2 <- cbind(ans2, data.frame(
          Aliphatics = sum(xx$HtoC >= 1.5 & xx$HtoC < 2.0 & xx$N ==0, na.rm = TRUE),
          AliphaticsN = sum(xx$HtoC >= 1.5 & xx$HtoC < 2.0 & xx$N > 0, na.rm = TRUE),
          Saturated = sum(xx$HtoC >= 2.0 | xx$OtoC >= 0.9, na.rm = TRUE), 
          Condensed_Aromatics = sum(xx$AImod > 0.66, na.rm = TRUE),
          Aromatic = sum(xx$AImod <= 0.66 & xx$AImod > 0.5, na.rm = TRUE),
          LigninPhenolics = sum(xx$AImod <= 0.5 & xx$HtoC < 2.0, na.rm = TRUE)
        ))
        ans2$Aromaticity_NA <- ans2$total_peaks-rowSums(ans2[, c('Aliphatics', 'AliphaticsN', 'Saturated', 'Condensed_Aromatics', 'Aromatic', 'LigninPhenolics')])
        }
      
      return(ans2)
    })
    
    compounds.df <- rbind.fill(compounds.df, temp.df)
  }
  
  
  write.csv(file=fileOut, compounds.df)
  #return(list(total_peakss=compounds.df, massInfo=massCharacteristic)) #massCharacteristic is 3.2 Mb
  return(compounds.df)
}

