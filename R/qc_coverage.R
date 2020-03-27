#' Read WGS Coverage Metrics File
#'
#' Reads the `wgs_coverage_metrics_<phenotype>.csv` file, which contains
#' read depth of coverage metrics.
#'
#' @param x Path to `wgs_coverage_metrics_<phenotype>.csv` file.
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
#' (cm_x <- read_wgs_coverage_metrics(x))
#' (cm_y <- read_wgs_coverage_metrics(y))
#'
#' @testexamples
#' expect_equal(colnames(cm_x), c("label", "var", "var_abbrev", "pct", "count"))
#' expect_equal(colnames(cm_y), c("label", "var", "var_abbrev", "pct", "count"))
#'
#' @export
read_wgs_coverage_metrics <- function(x) {

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

  b <- basename(x)
  suffix <- ifelse(grepl("_normal\\.csv", b), "_N",
                   ifelse(grepl("_tumor\\.csv", b), "_T",
                          ""))
  nm <- sub("(.*)\\.wgs_coverage_metrics.*", "\\1", b)
  label <- paste0(nm, suffix)

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
#' @param keep_alt Keep the ALT + Mito chromosomes?
#' @return tibble with following columns:
#'   - label
#'   - chrom
#'   - n_bases
#'   - coverage
#'
#' @examples
#' x <- system.file("extdata/COLO829.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' y <- system.file("extdata/COLO829.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#'
#' (cc_x <- read_wgs_contig_coverage(x))
#' (cc_y <- read_wgs_contig_coverage(y))
#'
#' @testexamples
#' expect_equal(colnames(cc_x), c("label", "chrom", "n_bases", "coverage"))
#' expect_equal(colnames(cc_y), c("label", "chrom", "n_bases", "coverage"))
#'
#' @export
read_wgs_contig_coverage <- function(x, keep_alt = FALSE) {

  b <- basename(x)
  suffix <- ifelse(grepl("_normal\\.csv", b), "_N",
                   ifelse(grepl("_tumor\\.csv", b), "_T",
                          ""))
  nm <- sub("(.*)\\.wgs_contig_mean_cov.*", "\\1", b)
  label <- paste0(nm, suffix)

  readr::read_csv(x, col_names = c("chrom", "n_bases", "coverage"), col_types = "cdd") %>%
    dplyr::filter(
      if (!keep_alt) {
        !grepl("chrM|MT|_|Autosomal|HLA-", .data$chrom)
      } else {
        TRUE
      }) %>%
    dplyr::mutate(label = label) %>%
    dplyr::select(.data$label, .data$chrom, .data$n_bases, .data$coverage)
}

#' Plot WGS Contig Mean Coverage Files
#'
#' Plots the `wgs_contig_mean_cov_<phenotype>.csv` files.
#'
#' @param xs Character vector containing paths to `wgs_contig_mean_cov_<phenotype>.csv` files.
#' @param top_alt_n Number of top covered alt contigs to plot per phenotype.
#'
#' @return A ggplot2 object with chromosomes on X axis, and coverage on Y axis.
#'
#' @examples
#' normal <- system.file("extdata/COLO829.wgs_contig_mean_cov_normal.csv.gz", package = "dracarys")
#' tumor <- system.file("extdata/COLO829.wgs_contig_mean_cov_tumor.csv.gz", package = "dracarys")
#' plot_wgs_contig_coverage(xs = c(tumor, normal))
#'
#' @export
plot_wgs_contig_coverage <- function(xs, top_alt_n = 15) {
  assertthat::assert_that(length(top_alt_n) == 1, top_alt_n >= 0, is.numeric(top_alt_n))

  cov_contig <-
    purrr::map(.x = xs, read_wgs_contig_coverage, keep_alt = TRUE) %>%
    dplyr::bind_rows()

  # Display chr1-22, X, Y at top (M goes to bottom).
  # Display top 20 of the rest, plus rest as 'other', at bottom
  main_chrom1 <- c(1:22, "X", "Y")
  main_chrom2 <- c(paste0("chr", main_chrom1))
  main_chrom <- c(main_chrom1, main_chrom2, "Autosomal regions")

  cov_contig <- cov_contig %>%
    dplyr::mutate(panel = ifelse(.data$chrom %in% main_chrom, "main", "alt"),
                  panel = factor(.data$panel, levels = c("main", "alt")))

  main_panel <- cov_contig %>%
    dplyr::filter(.data$panel == "main") %>%
    dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)
  alt_panel <- cov_contig %>%
    dplyr::filter(.data$panel == "alt") %>%
    dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)

  top_alt <- alt_panel %>%
    dplyr::group_by(.data$label) %>%
    dplyr::top_n(top_alt_n, wt = .data$coverage) %>%
    dplyr::arrange(dplyr::desc(.data$coverage)) %>%
    dplyr::pull(.data$chrom) %>%
    unique()

  alt_panel2 <- alt_panel %>%
    dplyr::mutate(alt_group = ifelse(.data$chrom %in% top_alt, "top", "bottom"))

  alt_panel_final <- alt_panel2 %>%
    dplyr::group_by(.data$alt_group, .data$label) %>%
    dplyr::summarise(mean_cov = mean(.data$coverage)) %>%
    dplyr::inner_join(alt_panel2, by = c("alt_group", "label")) %>%
    dplyr::mutate(
      chrom = ifelse(.data$alt_group == "bottom", "OTHER", .data$chrom),
      coverage = ifelse(.data$alt_group == "bottom", .data$mean_cov, .data$coverage)) %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::select(.data$label, .data$chrom, .data$coverage, .data$panel)

  chrom_fac_levels <- c(main_chrom, "chrM", "MT", top_alt[!top_alt %in% c("chrM", "MT")], "OTHER")
  d <- dplyr::bind_rows(main_panel, alt_panel_final) %>%
    dplyr::mutate(chrom = factor(.data$chrom, levels = chrom_fac_levels))

  d %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data$chrom, y = .data$coverage,
                   colour = .data$label, group = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(
      limits = c(0, NA), expand = c(0, 0), labels = scales::comma,
      breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Mean Coverage Per Chromosome", colour = "Label") +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Coverage") +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 6),
      plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"),
      panel.spacing = ggplot2::unit(2, "lines")) +
    ggplot2::facet_wrap(ggplot2::vars(.data$panel), nrow = 2, scales = "free")
}

#' Read WGS Fine Hist File
#'
#' Reads the `wgs_fine_hist_<phenotype>.csv` file, which contains
#' two columns: Depth and Overall.
#' The value in the Depth column ranges from 0 to 1000+ and the Overall
#' column indicates the number of loci covered at the corresponding depth.
#'
#' @param x Path to `wgs_fine_hist_<phenotype>.csv` file.
#' @return tibble with three columns:
#'   - label
#'   - depth
#'   - number of loci with given depth
#'
#' @examples
#' x <- system.file("extdata/COLO829.wgs_fine_hist_normal.csv.gz", package = "dracarys")
#' y <- system.file("extdata/COLO829.wgs_fine_hist_tumor.csv.gz", package = "dracarys")
#'
#' (fh_x <- read_wgs_fine_hist(x))
#' (fh_y <- read_wgs_fine_hist(y))
#'
#' @testexamples
#' expect_equal(colnames(fh_x), c("label", "depth", "n_loci"))
#' expect_equal(colnames(fh_y), c("label", "depth", "n_loci"))
#'
#' @export
read_wgs_fine_hist <- function(x) {
  d <- readr::read_csv(x, col_types = "cd")
  assertthat::assert_that(all(colnames(d) == c("Depth", "Overall")))

  b <- basename(x)
  suffix <- ifelse(grepl("_normal\\.csv", b), "_N",
                   ifelse(grepl("_tumor\\.csv", b), "_T",
                          ""))
  nm <- sub("(.*)\\.wgs_fine_hist.*", "\\1", b)
  label <- paste0(nm, suffix)

  d %>%
    dplyr::mutate(label = label,
                  Depth = ifelse(grepl("1000+", .data$Depth), 1000, .data$Depth),
                  Depth = as.integer(.data$Depth)) %>%
    dplyr::select(.data$label, depth = .data$Depth, n_loci = .data$Overall)
}


#' Plot WGS Fine Hist File
#'
#' Plots the `wgs_fine_hist_<phenotype>.csv` files.
#'
#' @param xs Character vector containing paths to `wgs_fine_hist_<phenotype>.csv` files.
#' @param x_lim X axis range to plot.
#'
#' @return A ggplot2 object with depth of coverage on X axis, and number of loci with
#'   that depth on Y axis.
#'
#' @examples
#' normal <- system.file("extdata/COLO829.wgs_fine_hist_normal.csv.gz", package = "dracarys")
#' tumor <- system.file("extdata/COLO829.wgs_fine_hist_tumor.csv.gz", package = "dracarys")
#'
#' plot_wgs_fine_hist(xs = c(tumor, normal))
#' plot_wgs_fine_hist(xs = c(tumor, normal), x_lim = c(0, 500))
#' @export
plot_wgs_fine_hist <- function(xs, x_lim = c(0, 300)) {
  assertthat::assert_that(length(x_lim) == 2)

  cov <-
    purrr::map(xs, read_wgs_fine_hist) %>%
    dplyr::bind_rows()

  cov %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$depth, y = .data$n_loci,
                                 colour = .data$label, group = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::coord_cartesian(xlim = x_lim) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Coverage Distribution", colour = "Label") +
    ggplot2::xlab("Depth of Coverage") +
    ggplot2::ylab("Number of Loci with Given Coverage") +
    ggplot2::theme(
      legend.position = c(0.9, 0.9),
      legend.justification = c(1, 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      # axis.text.x = ggplot2::element_text(angle = 0, vjust = 1, hjust = 1),
      plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"))
}


#' Plot WGS Fine Hist File Cumsum
#'
#' Plots the `wgs_fine_hist_<phenotype>.csv` coverage cumsum.
#'
#' @param xs Character vector containing paths to `wgs_fine_hist_<phenotype>.csv` files.
#' @param y_lim Y axis lower limit.
#'
#' @return A ggplot2 object with depth of coverage on X axis, and fraction
#'   of loci with that depth on Y axis.
#'
#' @examples
#' normal <- system.file("extdata/COLO829.wgs_fine_hist_normal.csv.gz", package = "dracarys")
#' tumor <- system.file("extdata/COLO829.wgs_fine_hist_tumor.csv.gz", package = "dracarys")
#'
#' plot_wgs_fine_cumsum(xs = c(tumor, normal))
#' @export
plot_wgs_fine_cumsum <- function(xs, y_lim = 0.003) {
  assertthat::assert_that(length(y_lim) == 1)
  cov <-
    purrr::map(xs, read_wgs_fine_hist) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$label) %>%
    dplyr::mutate(freq = .data$n_loci / sum(.data$n_loci)) %>%
    dplyr::arrange(dplyr::desc(.data$depth)) %>%
    dplyr::mutate(cumsum = cumsum(.data$freq)) %>%
    dplyr::filter(.data$cumsum >= y_lim) %>%
    dplyr::ungroup()

  cov %>%
    ggplot2::ggplot(
      ggplot2::aes(x = .data$depth, y = .data$cumsum, colour = .data$label)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Cumulative Coverage", colour = "Label") +
    ggplot2::xlab("Depth of Coverage") +
    ggplot2::ylab("Fraction of Loci with >= Given Coverage") +
    ggplot2::theme(
      legend.position = c(0.9, 0.9),
      legend.justification = c(1, 1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"))
}
