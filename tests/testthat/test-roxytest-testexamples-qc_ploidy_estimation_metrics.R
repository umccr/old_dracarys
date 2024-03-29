# Generated by roxytest: Do not edit by hand!

context("File R/qc_ploidy_estimation_metrics.R: @testexamples")

test_that("Function read_ploidy_estimation_metrics() @ L23", {
  
  x <- system.file("extdata/COLO829.ploidy_estimation_metrics.csv.gz", package = "dracarys")
  (pm <- read_ploidy_estimation_metrics(x))
  
  expect_equal(colnames(pm), c("label", "var", "value"))
})

