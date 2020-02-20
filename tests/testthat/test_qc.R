context("qc-metrics")

mm <- read_mapping_metrics(system.file("extdata/COLO829.mapping_metrics.csv.gz", package = "dracarys"))
fl <- read_fragment_length_hist(system.file("extdata/COLO829.fragment_length_hist.csv.gz", package = "dracarys"))
cm <- read_wgs_coverage_metrics(system.file("extdata/COLO829.wgs_coverage_metrics_tumor.csv.gz", package = "dracarys"),
                                label = "tumor")
tm <- read_time_metrics(system.file("extdata/COLO829.time_metrics.csv.gz", package = "dracarys"))
vm <- read_varcalling_metrics(system.file("extdata/COLO829.vc_metrics.csv.gz", package = "dracarys"))
replay <- read_replay(system.file("extdata/COLO829-replay.json.gz", package = "dracarys"))
ploidy_metrics <- read_ploidy_estimation_metrics(system.file("extdata/COLO829.ploidy_estimation_metrics.csv.gz", package = "dracarys"))
cov_contig <- read_wgs_contig_coverage(system.file("extdata/COLO829.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys"), label = "tumor")


test_that("column names are correct", {
  expect_equal(colnames(mm), c("category", "Phenotype", "RG", "var", "var_abbrev", "count", "pct"))
  expect_equal(colnames(fl), c("fragmentLength", "count", "sample"))
  expect_equal(colnames(cm), c("label", "var", "var_abbrev", "pct", "count"))
  expect_equal(colnames(tm), c("Step", "Time"))
  expect_equal(colnames(vm), c("category", "sample", "var", "count", "pct"))
  expect_equal(names(replay), c("command_line", "dragen_config", "inputs", "system"))
  expect_equal(colnames(ploidy_metrics), c("var", "value"))
  expect_equal(colnames(cov_contig), c("label", "chrom", "n_bases", "coverage"))
})
