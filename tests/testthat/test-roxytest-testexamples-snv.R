# Generated by roxytest: Do not edit by hand!

context("File R/snv.R: @testexamples")

test_that("Function snv_table() @ L25", {
  
  env_path <- paste(Sys.getenv("DRACARYS_ENV"),
                    system.file("src", package = "dracarys"), "$PATH", sep = ":")
  x <- system.file("extdata/COLO829.hard-filtered.vcf.gz", package = "dracarys")
  (snvs <- snv_table(x, prefix = "COLO829", outdir = tempdir(), env_path = env_path))
  
  expect_equal(length(snvs), 2)
  expect_equal(names(snvs), c("variants", "samples"))
  expect_true(inherits(snvs$variants, "tbl_df"))
})

