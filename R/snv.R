#' Generate Table with Single Nucleotide/INDEL Variants
#'
#' @param x Path to `hard-filtered.vcf.gz` file.
#' @param prefix Prefix for output variant and sample files.
#' @param outdir Output directory for results.
#' @param env_path Path to the directory containing the python environment
#' for snv_table.
#'
#' @return A list with two elements:
#'         - A tibble with SNV/INDEL variants.
#'         - A tibble with sample names corresponding to the VCF columns.
#'
#' @examples
#' env_path <- paste(Sys.getenv("DRACARYS_ENV"),
#'                   system.file("src", package = "dracarys"), "$PATH", sep = ":")
#' x <- system.file("extdata/COLO829.hard-filtered.vcf.gz", package = "dracarys")
#' (snvs <- snv_table(x, prefix = "COLO829", outdir = tempdir(), env_path = env_path))
#'
#' @testexamples
#' expect_equal(length(snvs), 2)
#' expect_equal(names(snvs), c("variants", "samples"))
#' expect_true(inherits(snvs$variants, "tbl_df"))
#'
#' @export
snv_table <- function(x, prefix, outdir, env_path) {

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = TRUE)
  }
  outdir <- normalizePath(outdir)
  .snv_table_run <- function() {
    assertthat::assert_that(file.exists(x), length(prefix) == 1, nchar(prefix) > 0)

    export_path <- glue::glue("export PATH={env_path}:$PATH")
    snv_table_py <- shQuote(system.file("src/snv_table", package = "dracarys"))
    cmd <- glue::glue("{export_path} && {snv_table_py} -i {x} -p {prefix} -o {outdir}")
    system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
  }

  .snv_table_run()

  ctypes <- c(
    CHROM = "c", START = "i", VARTYPE = "c", FILTER = "c", DP = "d", MQ = "d", TLOD = "d", NLOD = "d",
    GT1 = "c", GT2 = "c", AF1 = "c", AF2 = "c", AD1 = "c", AD2 = "c", DP1 = "d", DP2 = "d")

  test1 <- readr::read_lines(glue::glue("{outdir}/snv_table/{prefix}_variants.tsv"), n_max = 1)
  test2 <- readr::read_lines(glue::glue("{outdir}/snv_table/{prefix}_samples.tsv"), n_max = 1)
  assertthat::assert_that(test1 == paste(names(ctypes), collapse = "\t"),
                          test2 == paste(c("col1", "col2"), collapse = "\t"))
  variants <- readr::read_tsv(glue::glue("{outdir}/snv_table/{prefix}_variants.tsv"), col_types = paste(ctypes, collapse = ""))
  samples <- readr::read_tsv(glue::glue("{outdir}/snv_table/{prefix}_samples.tsv"), col_types = "cc")

  list(variants = variants,
       samples = samples)
}
