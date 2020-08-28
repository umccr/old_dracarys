#' Read ROH Metrics File
#'
#' @param x Path to `roh_metrics.csv` file.
#'
#' @return A tibble with two columns:
#'         - `pct_snv_in_3mb_roh`: percent of SNVs in large (>= 3Mb) ROH regions.
#'         - `n_3mb_roh`: number of large (>= 3Mb) ROH regions.
#'
#' @examples
#' x <- system.file("extdata/COLO829_N.roh_metrics.csv", package = "dracarys")
#' (roh <- read_roh_metrics(x))
#'
#' @testexamples
#' expect_equal(nrow(roh), 1)
#'
#' @export
read_roh_metrics <- function(x) {
  readr::read_csv(x, col_names = c("dummy1", "dummy2", "key", "value"),
                  col_types = "cccd") %>%
    dplyr::select(.data$key, .data$value) %>%
    dplyr::mutate(key = ifelse(.data$key == "Percent SNVs in large ROH ( >= 3000000)", "pct_snv_in_3mb_roh",
                               ifelse(.data$key == "Number of large ROH ( >= 3000000)", "n_3mb_roh", "OTHER"))) %>%
    tidyr::pivot_wider(names_from = .data$key, values_from = .data$value)

}

#' Read ROH BED File
#'
#' @param x Path to `roh.bed` file.
#'
#' @return A tibble with 4 columns:
#'         - `chrom`, `start`, `end`: coordinates for one region of homozygosity (0-based, half-open interval).
#'           The `chrom` is returned with a `hs` prefix instead of `chr`.
#'         - `score`: function of the number of homozygous and heterozygous variants,
#'           where each homozygous variant increases the score by a pre-defined value,
#'           and each heterozygous variant reduces the score by (1 - pre-defined value),
#'           where the pre-defined value is in range (0, 1).
#'
#' @examples
#' x <- system.file("extdata/COLO829_N.roh.bed", package = "dracarys")
#' (roh <- rohbed2circos(x))
#'
#' @testexamples
#' expect_equal(ncol(roh), 4)
#'
#' @export
rohbed2circos <- function(x) {
  readr::read_tsv(x, col_names = c("chrom", "start", "end", "value", "hom", "het"),
                  col_types = "cddddd") %>%
    dplyr::mutate(`#chromosome` = sub("chr", "hs", .data$chrom)) %>%
    dplyr::select(.data$`#chromosome`, .data$start, .data$end, .data$value)

}
