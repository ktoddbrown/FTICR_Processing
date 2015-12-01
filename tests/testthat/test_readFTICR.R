# Testing code for the RCMIP5 scripts in 'readFTICR.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("readFTICR")

test_that('readFTICR runs on test file with defaults', {
  temp <- readFTICR(fileIn='testdata/smallTest.csv', samplesToRead=1:60)
  expect_equal(names(temp), c('m.z', 'sample', 'intensity'))
  expect_equal(dim(temp), c(108, 3))
})

test_that('readFTICR errors on bad filename',{
  expect_error(readFTICR(fileIn='badFileName.csv'))
})

test_that('readFTICR works with index, regexp, and arrays',{
  regExpRead <- readFTICR(fileIn='testdata/smallTest.csv',
                    samplesToRead='X.out2.')
  namedRead <- readFTICR(fileIn='testdata/smallTest.csv',
                    samplesToRead=c('X.out2.', 'X.out20.', 'X.out21.', 'X.out22.', 'X.out23.', 'X.out24.', 'X.out25.', 'X.out26.', 'X.out27.', 'X.out28.', 'X.out29.'))
  indexRead <- readFTICR(fileIn='testdata/smallTest.csv',
                     samplesToRead=12:22)
  
  expect_equal(indexRead, namedRead)
  expect_equal(indexRead, regExpRead)
})