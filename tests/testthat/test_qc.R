context("qc-metrics")

mm <- read_mapping_metrics(system.file("extdata/COLO829.mapping_metrics.csv.gz", package = "dracarys"))
fl <- read_fragment_length_hist(system.file("extdata/COLO829.fragment_length_hist.csv.gz", package = "dracarys"))


test_that("column names are correct", {
  expect_equal(colnames(mm), c("category", "Phenotype", "RG", "var", "var_abbrev", "count", "pct"))
  expect_equal(colnames(fl), c("fragmentLength", "count", "sample"))
})
