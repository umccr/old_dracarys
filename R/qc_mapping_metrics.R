#' Read Mapping Metrics File
#'
#' Reads the `mapping_metrics.csv` file, which contains
#' mapping and aligning metrics, like the metrics computed by the Samtools
#' Flagstat command. These metrics are available on an aggregate level
#' (over all input data), and on a per read group level.
#' Unless explicitly stated, the metrics units are in reads (i.e., not in
#' terms of pairs or alignments).
#'
#' @param x Path to `mapping_metrics.csv` file.
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
#' x <- system.file("extdata/COLO829.mapping_metrics.csv.gz", package = "dracarys")
#' (mm <- read_mapping_metrics(x))
#'
#' @testexamples
#' expect_equal(colnames(mm), c("category", "Phenotype", "RG", "var", "var_abbrev", "count", "pct"))
#'
#'
#' @export
read_mapping_metrics <- function(x) {

  abbrev_nm <- c(
    "Total Reads per RG" = "Tot",
    "Number of duplicate marked reads" = "Dup",
    "Number of duplicate marked and mate reads removed" = "Dup Rem",
    "Number of unique reads (excl. duplicate marked reads)" = "Unique",
    "Reads with mate sequenced" = "Mated",
    "Reads without mate sequenced" = "noMated",
    "QC-failed reads" = "Failed",
    "Mapped reads" = "Mapped",
    "Mapped reads R1" = "R1map",
    "Mapped reads R2" = "R2map",
    "Number of unique & mapped reads (excl. duplicate marked reads)" = "UniqueMap",
    "Unmapped reads" = "Unmapped",
    "Singleton reads (itself mapped; mate unmapped)" = "Singleton",
    "Paired reads (itself & mate mapped)" = "PairedMap",
    "Properly paired reads" = "PairedProper",
    "Not properly paired reads (discordant)" = "PairedDisc",
    "Paired reads mapped to different chromosomes" = "DiffChrom",
    "Paired reads mapped to different chromosomes (MAPQ>=10)" = "DiffChrom MQ10",
    "Reads with MAPQ [40:inf)" = "MQ 40+",
    "Reads with MAPQ [30:40)" = "MQ 30-40",
    "Reads with MAPQ [20:30)" = "MQ 20-30",
    "Reads with MAPQ [10:20)" = "MQ 10-20",
    "Reads with MAPQ [ 0:10)" = "MQ 0-10",
    "Reads with MAPQ NA (Unmapped reads)" = "MQ NA",
    "Reads with indel R1" = "R1 indel",
    "Reads with indel R2" = "R2 indel",
    "Total bases" = "TotalBases",
    "Total bases R1" = "TotBasesR1",
    "Total bases R2" = "TotBasesR2",
    "Mapped bases R1" = "MappedBasesR1",
    "Mapped bases R2" = "MappedBasesR2",
    "Soft-clipped bases R1" = "SoftClipR1",
    "Soft-clipped bases R2" = "SoftClipR2",
    "Mismatched bases R1" = "MismatchR1",
    "Mismatched bases R2" = "MismatchR2",
    "Mismatched bases R1 (excl. indels)" = "MismatchR1 NI",
    "Mismatched bases R2 (excl. indels)" = "MismatchR2 NI",
    "Q30 bases" = "Q30",
    "Q30 bases R1" = "R1Q30",
    "Q30 bases R2" = "R2Q30",
    "Q30 bases (excl. dups & clipped bases)" = "Q30 nondup",
    "Total alignments" = "TotAlign",
    "Secondary alignments" = "SecAlign",
    "Supplementary (chimeric) alignments" = "ChimericAlign",
    "Estimated read length" = "Read Length",
    "Bases in reference genome" = "Genome Bases",
    "Bases in target bed [% of genome]" = "BedBases %Genome",
    "Average sequenced coverage over genome" = "Coverage Avg",
    "Insert length: mean" = "InsertLength Mean",
    "Insert length: median" = "InsertLength Median",
    "Insert length: standard deviation" = "InsertLength StdDev",
    "Provided sex chromosome ploidy" = "Ploidy SexChrom",
    "DRAGEN mapping rate [mil. reads/second]" = "Map Rate"
  )

  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("MAPPING/ALIGNING", d[1]))

  d %>%
    tibble::enframe(name = "name", value = "value") %>%
    tidyr::separate(.data$value, into = c("category", "RG", "extra"), sep = ",", extra = "merge") %>%
    tidyr::separate(.data$extra, into = c("var", "value"), sep = ",", extra = "merge") %>%
    tidyr::separate(.data$value, into = c("count", "pct"), sep = ",", fill = "right", convert = TRUE) %>%
    dplyr::mutate(
      Phenotype = dplyr::case_when(
        grepl("TUMOR", .data$category) ~ "tumor",
        grepl("NORMAL", .data$category) ~ "normal",
        TRUE ~ "unknown"),
      category = dplyr::case_when(
        grepl("ALIGNING SUMMARY", .data$category) ~ "summary",
        grepl("ALIGNING PER RG", .data$category) ~ "readgroup",
        TRUE ~ "unknown"),
      RG = ifelse(.data$RG == "", "TOTAL", .data$RG),
      var = ifelse(grepl("Total.*reads", .data$var), "Total Reads per RG", .data$var),
      var_abbrev = dplyr::recode(.data$var, !!!abbrev_nm)) %>%
    dplyr::select(.data$category, .data$Phenotype, .data$RG,
                  .data$var, .data$var_abbrev, .data$count, .data$pct)
}
