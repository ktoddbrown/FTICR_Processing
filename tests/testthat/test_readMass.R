# Testing code for the RCMIP5 scripts in 'readMass.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("readMass")

test_that('readMass pulls the correct columns', {
  temp <- readMass(fileIn='testdata/smallTest.csv')
  expect_equal(names(temp), c('m.z', 'C', 'H', 'O', 'N', 'X13C', 'S', 'P', 'OtoC', 'HtoC', 'AImod'))
  expect_equal(dim(temp), c(49, 11))
})