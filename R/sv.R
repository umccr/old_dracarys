#' Plot SVs in Circos
#'
#' Generates a Circos plot of  the SVs found in the `sv.vcf.gz` VCF file output
#' by DRAGEN.
#'
#' @param manta Path to `sv.vcf.gz` file.
#' @param circos_path Path to the directory containing the circos binary (see examples).
#' @param outdir Directory to write the output circos file along with the configs used.
#' @param outfile Name of PNG file.
#' @param genome Genome assembly used to generate the SV calls.
#' @param filter_pass Logical. Should we keep only the PASSed variants?
#'
#' @return Path to circos PNG file output by the Perl circos program.
#'
#' @examples
#' \dontrun{
#' manta_sv <- system.file("extdata/COLO829.sv.vcf.gz", package = "dracarys")
#'
#' circos_png <- plot_circos(
#'    manta = manta_sv,
#'    circos_path = "path/to/bin_with_circos",
#'    outdir = file.path("circos"))
#' }
#' @export
plot_circos <- function(manta, circos_path, outdir, outfile = "circos.png", genome = "hg38", filter_pass = TRUE) {
  sv_vcf_circos <- rock::circos_prep(manta = manta, genome = genome, outdir = outdir, filter_pass = filter_pass)
  export_path <- glue::glue("export PATH={circos_path}:$PATH")
  cmd <- glue::glue("{export_path} && ",
                    "circos -nosvg -conf {outdir}/circos.conf ",
                    "-outputdir {outdir} -outputfile {outfile}")
  system(cmd, ignore.stdout = FALSE)
  return(file.path(outdir, outfile))
}


#' Plot SV Counts per Contig
#'
#' Plots number of SV breakpoints per contig for mate1 and mate2 of each
#' SV event.
#'
#' @param manta Path to `sv.vcf.gz` file or object of class 'sv' from the
#' `rock` package.
#'
#' @return A ggplot2 object plotting number of SV breakpoints per contig for
#' mate1 and mate2 of each SV event.
#'
#' @examples
#' manta <- system.file("extdata/COLO829.sv.vcf.gz", package = "dracarys")
#' plot_sv_contig_counts(manta)
#'
#' manta_rock <- rock::prep_manta_vcf(manta, filter_pass = TRUE)
#' plot_sv_contig_counts(manta_rock)
#'
#' @export
plot_sv_contig_counts <- function(manta) {
  if (inherits(manta, "sv")) {
    sv_vcf <- manta$sv
  } else {
    sv_vcf <- rock::prep_manta_vcf(manta, filter_pass = TRUE)$sv
  }

  noalt_chrom <- c(1:22, "X", "Y", "M")
  alt_label <- "ALT"

  d <- sv_vcf %>%
    dplyr::select(.data$chrom1, .data$chrom2) %>%
    tidyr::pivot_longer(cols = c(.data$chrom1, .data$chrom2), names_to = "mate", values_to = "chromosome") %>%
    dplyr::mutate(is_alt = !.data$chromosome %in% noalt_chrom,
                  chromosome = ifelse(.data$is_alt, alt_label, .data$chromosome)) %>%
    dplyr::count(.data$mate, .data$chromosome) %>%
    dplyr::mutate(chromosome = factor(.data$chromosome, levels = c(noalt_chrom, alt_label)))

  max_n <- max(d$n)

  ggplot2::ggplot(data = d, ggplot2::aes(x = .data$chromosome, y = .data$n, group = .data$mate, colour = .data$mate)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 10), limits = c(0, max_n + 1), expand = c(0, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Count of Breakpoints per Chromosome", colour = "Mate") +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Count") +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
      plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"),
      panel.spacing = ggplot2::unit(2, "lines"))
}
