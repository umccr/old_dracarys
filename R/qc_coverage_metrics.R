#' Read WGS Coverage Metrics File
#'
#' Reads the `wgs_coverage_metrics_<phenotype>.csv` file, which contains
#' read depth of coverage metrics.
#'
#' @param x Path to `wgs_coverage_metrics_<phenotype>.csv` file.
#' @param label Label for the file e.g. 'tumor' or 'normal'.
#' @return tibble with following columns:
#'   - label: file label.
#'   - var: variable name.
#'   - var_abbrev: variable abbreviation.
#'   - pct: percentage value.
#'   - count: count value.
#'
#' @examples
#' x <- system.file("extdata/COLO829.wgs_coverage_metrics_normal.csv.gz", package = "dracarys")
#' y <- system.file("extdata/COLO829.wgs_coverage_metrics_tumor.csv.gz", package = "dracarys")
#'
#' read_wgs_coverage_metrics(x, label = "normal")
#' read_wgs_coverage_metrics(y, label = "tumor")
#' @export
read_wgs_coverage_metrics <- function(x, label) {

  abbrev_nm <- c(
    "Aligned bases" = "Aln bases",
    "Aligned bases in genome" = "Aln Bases Genome",
    "Average alignment coverage over genome" = "Avg Cov Genome",
    "Average chr X coverage over genome" = "Avg Cov chrX",
    "Average chr Y coverage over genome" = "Avg Cov chrY",
    "Average mitochondrial coverage over genome" = "Avg Cov chrM",
    "Average autosomal coverage over genome" = "Avg Cov Autos",
    "Median autosomal coverage over genome" = "Med Cov Autos",
    "Mean/Median autosomal coverage ratio over genome" = "Cov Ratio",
    "Aligned reads" = "Aln Reads",
    "Aligned reads in genome" = "Aln Reads Genome")

  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("COVERAGE SUMMARY", d[1]))

  d <- d %>%
    tibble::enframe(name = "name", value = "value") %>%
    tidyr::separate(.data$value, into = c("category", "dummy1", "extra"), sep = ",", extra = "merge") %>%
    tidyr::separate(.data$extra, into = c("var", "value"), sep = ",", extra = "merge") %>%
    dplyr::mutate(label = label)

  pct <- d %>%
    dplyr::filter(grepl("PCT", .data$var)) %>%
    dplyr::mutate(
      value = as.numeric(.data$value),
      var_abbrev = dplyr::case_when(
        grepl("PCT of genome", .data$var) ~ sub("PCT of genome with coverage", "%genome", var),
        grepl("Uniformity", .data$var) ~ "uniformity (% > 0.2*mean)",
        TRUE ~ "FOO")) %>%
    dplyr::select(.data$label, .data$var, .data$var_abbrev, pct = .data$value)

  cnt <- d %>%
    dplyr::filter(!grepl("PCT", .data$var)) %>%
    tidyr::separate(.data$value, into = c("count", "pct"), sep = ",", fill = "right", convert = TRUE) %>%
    dplyr::mutate(var_abbrev = dplyr::recode(.data$var, !!!abbrev_nm)) %>%
    dplyr::select(.data$label, .data$var, .data$var_abbrev, .data$count, .data$pct)

  dplyr::bind_rows(pct, cnt)
}
