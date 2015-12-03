# Testing code for the RCMIP5 scripts in 'readFTICR.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("readFTICR")

test_that('readFTICR runs on test file with defaults for sample only read', {
  temp <- readFTICR(fileIn='testdata/smallTest.csv', massHeader = 'm.z', 
                    sampleRegStr='X.out', samplesToRead=1:60)
  expect_equal(names(temp), c('m.z', 'sample', 'intensity'))
  expect_equal(dim(temp), c(108, 3))
})

test_that('readFTICR errors on bad filename for sample only read',{
  expect_error(readFTICR(fileIn='badFileName.csv'))
})

test_that('readFTICR works with index, regexp, and arrays for sample only read',{
  regExpRead <- readFTICR(fileIn='testdata/smallTest.csv', massHeader = 'm.z', sampleRegStr='X.out',
                    samplesToRead='X.out2.')
  namedRead <- readFTICR(fileIn='testdata/smallTest.csv', massHeader = 'm.z', sampleRegStr='X.out',
                    samplesToRead=c('X.out2.', 'X.out20.', 'X.out21.', 'X.out22.', 'X.out23.', 'X.out24.', 'X.out25.', 'X.out26.', 'X.out27.', 'X.out28.', 'X.out29.'))
  indexRead <- readFTICR(fileIn='testdata/smallTest.csv', massHeader = 'm.z', sampleRegStr='X.out',
                     samplesToRead=12:22)
  
  expect_equal(indexRead, namedRead)
  expect_equal(indexRead, regExpRead)
})

####Test mass only reads
test_that('readFTICR pulls the correct columns for mass only read', {
  temp <- readFTICR(fileIn='testdata/smallTest.csv', massHeader='m.z',
                   elementKey = list(C=c('C', 'X13C'), H='H', O='O', N='N', S='S', P='P'))
  expect_equal(names(temp), c('m.z', 'C', 'H', 'O', 'N', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod'))
  expect_equal(dim(temp), c(49, 12))
})

test_that('readFTICR handles different element naming conventions for mass only read', {
  temp <- readFTICR(fileIn='testdata/smallTest.csv',
                   massHeader='m.z',
                   elementKey = list(C=c('C', 'X13C'), 
                                     H='H', O='O', N='N', S='S', P='P'))
  temp2 <-  readFTICR(fileIn='testdata/smallTestAltElm.csv',
                     massHeader='m.z',
                     elementKey = list(C=c('numC', 'num13C'), 
                                       H='numH', O='numO', N='numN', S='numS', P='numP'))
  temp <- temp[,c('m.z', 'C', 'H', 'O', 'N', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod')]
  temp2 <- temp2[,c('m.z', 'C', 'H', 'O', 'N', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod')]
  
  expect_identical(temp, temp2)
})
test_that('readFTICR warns if elements not present for mass only read', {
  expect_warning(readFTICR(fileIn='testdata/smallTest.csv',
                          massHeader='m.z',
                          elementKey = list(C='C')), '^O, H, and C not found')
  expect_warning(readFTICR(fileIn='testdata/smallTest.csv',
                          massHeader='m.z',
                          elementKey = list(C='C', H='H', O='O')), '^N, S, and P not found')
})

test_that('readFTICR errors if there is a bad element name passed for mass only read',{
  expect_error(readFTICR(fileIn='testdata/smallTest.csv',
                        massHeader='m.z',
                        elementKey = list(C='C', H='H', O='numO'), 'column name in elementKey does not match'))
})
