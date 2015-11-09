# Testing code for the RCMIP5 scripts in 'readFTICR.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("readFTICR")

test_that('readFTICR pulls the correct columns', {
  temp <- readFTICR(fileIn='testdata/smallTest.csv')
  expect_equal(names(temp), c('m.z', 'C', 'H', 'O', 'N', 'X13C', 'S', 'P', 'OtoC', 'HtoC', 'AImod'))
  expect_equal(dim(temp), c(49, 11))
})