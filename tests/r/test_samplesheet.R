testthat::test_that('samplesheet test fixture declares required columns', {
  x <- utils::read.csv('tests/test_data/samplesheet.csv', check.names = FALSE)
  testthat::expect_true(all(c('sample_id', 'idat_red', 'idat_green') %in% names(x)))
})
