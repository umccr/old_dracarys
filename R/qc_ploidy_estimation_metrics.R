#' Read Ploidy Estimation Metrics File
#'
#' Reads the `ploidy_esimation_metrics.csv` file
#'
#' @param x Path to `ploidy_estimation_metrics.csv` file.
#' @return tibble with the following columns:
#'     - category: summary or read group
#'     - Phenotype: e.g. tumor, normal
#'     - RG: read group
#'     - var: metric variable
#'     - var_abbrev: metric variable abbreviation
#'     - count: count of reads
#'     - pct: percentage of reads
#'
#' @examples
#' x <- system.file("extdata/COLO829.ploidy_estimation_metrics.csv.gz", package = "dracarys")
#' read_ploidy_estimation_metrics(x)
#' @export
read_ploidy_estimation_metrics <- function(x) {

  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("PLOIDY ESTIMATION", d[1]))

  d %>%
    tibble::enframe(name = "name", value = "value") %>%
    tidyr::separate(.data$value, into = c("dummy1", "dummy2", "var", "value"), sep = ",", convert = FALSE) %>%
    dplyr::select(.data$var, .data$value)

}
