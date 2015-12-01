# Testing code for the RCMIP5 scripts in 'countCompoundTypes.R'

# Uses the testthat package
# See http://journal.r-project.org/archive/2011-1/RJournal_2011-1_Wickham.pdf

context("countCompoundTypes")

test_that('countCompoundTypes executes with defaults', {
  expect_is(countCompoundTypes(fileIn='testdata/smallTest.csv'), 'data.frame')
})

test_that('countCompoundTypes returns the correct headers', {

  expect_named(countCompoundTypes(fileIn='testdata/smallTest.csv', calculateClass=c('compound', 'molecular', 'aromaticity')), c('sample', 'total_peaks', 
                       c("aveOtoC", "aveHtoC"), 
                       c("CHO", "CHON", "CHOS", "CHOP", "CHONS", "CHONP", "CHOSP", "CHONSP", 'Molecular_NA'),
                       c('Lipids', 'UnSaturated_Hydrocarbons', 'Condensed_Hydrocarbons', 'Proteins', 'Amino_Sugars', 'Carbohydrates', 'Lignin', 'Tannins', 'Compounds_NA'),
                       c('Aliphatics', 'AliphaticsN', 'Saturated', 'Condensed_Aromatics', 'Aromatic', 'LigninPhenolics', 'Aromaticity_NA')))
  
  expect_named(countCompoundTypes(fileIn='testdata/smallTest.csv', calculateClass=c('compound')), c('sample', 'total_peaks', 
                       c("aveOtoC", "aveHtoC"), 
                       c('Lipids', 'UnSaturated_Hydrocarbons', 'Condensed_Hydrocarbons', 'Proteins', 'Amino_Sugars', 'Carbohydrates', 'Lignin', 'Tannins', 'Compounds_NA')))
  
  expect_named(countCompoundTypes(fileIn='testdata/smallTest.csv', calculateClass=c('molecular')), c('sample', 'total_peaks', c("aveOtoC", "aveHtoC"), c("CHO", "CHON", "CHOS", "CHOP", "CHONS", "CHONP", "CHOSP", "CHONSP", 'Molecular_NA')))
  
  expect_named(countCompoundTypes(fileIn='testdata/smallTest.csv', calculateClass=c('aromaticity')), c('sample', 'total_peaks', 
                                                                                                                            c("aveOtoC", "aveHtoC"), c('Aliphatics', 'AliphaticsN', 'Saturated', 'Condensed_Aromatics', 'Aromatic', 'LigninPhenolics', 'Aromaticity_NA')))
})