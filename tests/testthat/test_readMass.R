# Testing code for the RCMIP5 scripts in 'readMass.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("readMass")

test_that('readMass pulls the correct columns', {
  temp <- readMass(fileIn='testdata/smallTest.csv', massHeader='m.z',
                   elementKey = list(C=c('C', 'X13C'), 
                                     H='H', O='O', N='N', S='S', P='P'))
  expect_equal(names(temp), c('m.z', 'C', 'H', 'O', 'N', 'X13C', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod'))
  expect_equal(dim(temp), c(49, 13))
})

test_that('readMass handles different element naming conventions', {
  temp <- readMass(fileIn='testdata/smallTest.csv',
                   massHeader='m.z',
                   elementKey = list(C=c('C', 'X13C'), 
                                     H='H', O='O', N='N', S='S', P='P'))
  temp2 <-  readMass(fileIn='testdata/smallTestAltElm.csv',
                     massHeader='m.z',
                     elementKey = list(C=c('numC', 'num13C'), 
                             H='numH', O='numO', N='numN', S='numS', P='numP'))
  temp <- temp[,c('m.z', 'C', 'H', 'O', 'N', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod')]
  temp2 <- temp2[,c('m.z', 'C', 'H', 'O', 'N', 'S', 'P', 'OtoC', 'HtoC', 'DBE', 'AI', 'AImod')]
  
  expect_identical(temp, temp2)
})
test_that('readMass warns if elements not present', {
  expect_warning(readMass(fileIn='testdata/smallTest.csv',
                   massHeader='m.z',
                   elementKey = list(C='C')), '^O, H, and C not found')
  expect_warning(readMass(fileIn='testdata/smallTest.csv',
                          massHeader='m.z',
                          elementKey = list(C='C', H='H', O='O')), '^N, S, and P not found')
})

test_that('readMass errors if there is a bad element name passed',{
  expect_error(readMass(fileIn='testdata/smallTest.csv',
                        massHeader='m.z',
                        elementKey = list(C='C', H='H', O='numO'), 'column name in elementKey does not match'))
})
