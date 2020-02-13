#' Plot SVs
#'
#' Plots the SVs found in the `sv.vcf.gz` VCF file output by DRAGEN.
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
#' circos_png <- plot_circos(
#'    manta = "path/to/sample.sv.vcf",
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
  system(cmd, ignore.stdout = TRUE)
  return(file.path(outdir, outfile))
}
