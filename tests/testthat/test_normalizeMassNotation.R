# Testing code for the RCMIP5 scripts in 'normalizeMassNotation.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("normalizeMassNotation")

test_that('normalizeMassNotation runs', {
  data.df <- data.frame(m.z=c(123, 132), C=c(33, 32), H=c(40, 18), O=c(38, 37), N=c(0, 2), X13C=c(1, 0), S=c(1, 0), P=c(0, 1))
  expect_silent(temp <- normalizeMassNotation(data.df, sampleIDs=c('m.z'), 
                        elementKey = list(C='C', H='H', O='O', N='N', S='S', P='P')))
})