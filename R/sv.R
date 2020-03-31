#' Plot SVs in Circos
#'
#' Generates a Circos plot of  the SVs found in the `sv.vcf.gz` VCF file output
#' by DRAGEN.
#'
#' @param manta Path to `sv.vcf.gz` file.
#' @param env_path Path to the directory containing the circos binary (see examples).
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
#'    env_path = "path/to/bin_with_circos",
#'    outdir = file.path("circos"))
#' }
#' @export
plot_circos <- function(manta, env_path, outdir, outfile = "circos.png", genome = "hg38", filter_pass = TRUE) {
  rock::circos_prep(manta = manta, genome = genome, outdir = outdir, filter_pass = filter_pass)
  export_path <- glue::glue("export PATH={env_path}:$PATH")
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
#' @param manta Path to `sv.vcf.gz` file or object of class 'sv'.
#'
#' @return A ggplot2 object plotting number of SV breakpoints per contig for
#' mate1 and mate2 of each SV event.
#'
#' @examples
#' manta <- system.file("extdata/COLO829.sv.vcf.gz", package = "dracarys")
#' plot_sv_contig_counts(manta)
#'
#' manta_table <- dracarys::sv_table(manta)
#' plot_sv_contig_counts(manta_table)
#'
#' @export
plot_sv_contig_counts <- function(manta) {
  if (inherits(manta, "sv")) {
    sv <- manta$sv
  } else {
    sv <- sv_table(manta)$sv
  }

  noalt_chrom <- c(1:22, "X", "Y", "M")
  alt_label <- "ALT"

  d <- sv %>%
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

#' Generate Table with Structural Variants
#'
#' @param manta Path to `sv.vcf.gz` file or object of class 'sv'.
#'
#' @return A tibble with structural variants
#'
#' @examples
#' manta <- system.file("extdata/COLO829.sv.vcf.gz", package = "dracarys")
#' sv_table(manta)
#'
#' @export
sv_table <- function(manta) {
  x <- bedr::read.vcf(manta, split.info = TRUE, verbose = FALSE)
  standard_vcf_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "FORMAT")
  format_fields_description <- structure(x$header$FORMAT[, "Description"], names = x$header$FORMAT[, "ID"])
  info_fields_description <- structure(x$header$INFO[, "Description"], names = x$header$INFO[, "ID"])
  req_format_fields <- c("PR", "SR")
  req_info_fields <- c("MATEID", "SVTYPE", "SOMATICSCORE", "IMPRECISE", "CIPOS",
                       "CIGAR", "BND_DEPTH", "MATE_BND_DEPTH")
  sample_nms <- setdiff(colnames(x$vcf), c(standard_vcf_cols, names(info_fields_description)))

  assertthat::assert_that(length(sample_nms) > 0)
  assertthat::assert_that(all(req_format_fields == names(format_fields_description)))
  assertthat::assert_that(all(req_info_fields %in% names(info_fields_description)))

  sample_cols <- tibble::as_tibble(x$vcf[sample_nms]) %>%
    dplyr::mutate(num = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = sample_nms, names_to = "sample") %>%
    tidyr::separate(col = .data$value, into = req_format_fields, sep = ":", fill = "right") %>%
    tidyr::pivot_longer(cols = req_format_fields, names_to = "field") %>%
    tidyr::pivot_wider(names_from = c(.data$sample, .data$field), values_from = .data$value, id_cols = .data$num) %>%
    dplyr::arrange(.data$num)

  # bring them all together
  DF <- tibble::tibble(chrom1 = as.character(sub("chr", "", x$vcf$CHROM)),
                       pos1 = as.integer(x$vcf$POS),
                       pos2 = as.integer(x$vcf$END),
                       filter = x$vcf$FILTER,
                       id = x$vcf$ID) %>%
    dplyr::bind_cols(x$vcf[req_info_fields]) %>%
    dplyr::bind_cols(sample_cols) %>%
    dplyr::select(.data$num, dplyr::everything())

  # BNDs
  df_bnd <- DF %>%
    dplyr::filter(.data$SVTYPE == "BND") %>%
    dplyr::bind_cols(., .[match(.$id, .$MATEID), c("chrom1", "pos1")]) %>%
    dplyr::rename(chrom2 = .data$chrom11) %>%
    dplyr::mutate(pos2 = ifelse(is.na(.data$pos2), .data$pos11, .data$pos2),
                  bndid = substring(.data$id, nchar(.data$id)))

  orphan_mates <- df_bnd %>%
    dplyr::filter(.data$chrom2 %in% NA) %>%
    dplyr::mutate(orphan = paste0(.data$chrom1, ":", .data$pos1)) %>%
    dplyr::pull(.data$orphan)

  df_bnd <- df_bnd %>%
    dplyr::filter(!is.na(.data$chrom2)) %>%
    dplyr::filter(.data$bndid == "1") %>%
    dplyr::select(-c(.data$bndid, .data$pos11)) %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$chrom2,
                  .data$pos2, .data$id, .data$MATEID, .data$SVTYPE, .data$filter, dplyr::everything())

  if (length(orphan_mates) > 0) {
    warning(glue::glue("The following {length(orphan_mates)} orphan BND ",
                       "mates are removed:\n",
                       paste(orphan_mates, collapse = "\n")))
  }

  assertthat::assert_that(rock:::.manta_proper_pairs(df_bnd$id, df_bnd$MATEID))

  # Non-BNDs
  df_other <- DF %>%
    dplyr::filter(.data$SVTYPE != "BND") %>%
    dplyr::mutate(chrom2 = .data$chrom1) %>%
    dplyr::select(.data$chrom1, .data$pos1, .data$chrom2, .data$pos2,
                  .data$id, .data$MATEID, .data$SVTYPE, .data$filter, dplyr::everything())

  # All together now
  sv <- df_other %>%
    dplyr::bind_rows(df_bnd)

  structure(list(sv = sv), class = "sv")

}
