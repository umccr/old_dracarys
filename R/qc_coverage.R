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

#' Read WGS Contig Mean Coverage File
#'
#' Reads the `wgs_contig_mean_cov_<phenotype>.csv` file, which contains
#' the estimated coverage for all contigs, and an autosomal estimated coverage.
#' Columns are:
#'   - contig name.
#'   - number of bases aligned to contig
#'     (excludes bases from duplicate marked reads, reads with MAPQ=0, and clipped bases).
#'   - coverage estimate = col2 / contig length
#'
#' @param x Path to `wgs_contig_mean_cov_<phenotype>.csv` file.
#' @param phenotype Phenotype of file e.g. 'tumor', 'normal'.
#' @param keep_alt Keep the ALT + Mito chromosomes?
#' @return tibble with following columns:
#'   - phenotype
#'   - chrom
#'   - n_bases
#'   - coverage
#'
#' @examples
#' x <- system.file("extdata/COLO829.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' y <- system.file("extdata/COLO829.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#'
#' read_wgs_contig_coverage(x, phenotype = "normal")
#' read_wgs_contig_coverage(y, phenotype = "tumor")
#' @export
read_wgs_contig_coverage <- function(x, phenotype, keep_alt = FALSE) {
  readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") %>%
    dplyr::filter(
      if (!keep_alt) {
        !grepl("chrM|MT|_|Autosomal", .data$chrom)
      } else {
        TRUE
      }) %>%
    dplyr::mutate(phenotype = phenotype) %>%
    dplyr::select(.data$phenotype, .data$chrom, .data$n_bases, .data$coverage)
}

#' Plot WGS Contig Mean Coverage Files
#'
#' Plots the `wgs_contig_mean_cov_<phenotype>.csv` files for tumor and normal.
#'
#' @param tumor Path to `wgs_contig_mean_cov_tumor.csv` file.
#' @param normal Path to `wgs_contig_mean_cov_normal.csv` file.
#' @param colours Colours for normal and tumor sample, in that order.
#'
#' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
#'
#' @examples
#' normal <- system.file("extdata/COLO829.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' tumor <- system.file("extdata/COLO829.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#'
#' plot_wgs_contig_coverage(tumor = tumor, normal = normal)
#' @export
plot_wgs_contig_coverage <- function(tumor, normal, colours = c("#009E73", "#D55E00")) {
  assertthat::assert_that(length(colours) == 2)
  cov_contig_normal <- dracarys::read_wgs_contig_coverage(normal, phenotype = "normal")
  cov_contig_tumor <- dracarys::read_wgs_contig_coverage(tumor, phenotype = "tumor")

  cov_contig <- dplyr::bind_rows(cov_contig_normal, cov_contig_tumor)
  chrom_order <- gtools::mixedsort(unique(cov_contig$chrom))

  cov_contig %>%
    dplyr::mutate(chrom = factor(.data$chrom, levels = chrom_order)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$chrom, y = .data$coverage,
                                 colour = .data$phenotype, group = .data$phenotype)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(colour = "Phenotype") +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Coverage") +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))
}

#' Read WGS Fine Hist File
#'
#' Reads the `wgs_fine_hist_<phenotype>.csv` file, which contains
#' two columns: Depth and Overall.
#' The value in the Depth column ranges from 0 to 1000+ and the Overall
#' column indicates the number of loci covered at the corresponding depth.
#'
#' @param x Path to `wgs_fine_hist_<phenotype>.csv` file.
#' @param label Label for the file e.g. 'tumor' or 'normal'.
#' @return tibble with following columns:
#'
#' @examples
#' x <- system.file("extdata/COLO829.wgs_coverage_metrics_normal.csv.gz", package = "dracarys")
#' y <- system.file("extdata/COLO829.wgs_coverage_metrics_tumor.csv.gz", package = "dracarys")
#'
#' read_wgs_coverage_metrics(x, label = "normal")
#' read_wgs_coverage_metrics(y, label = "tumor")
#' @export
read_wgs_fine_hist <- function(x, label) {
  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("COVERAGE SUMMARY", d[1]))
}
