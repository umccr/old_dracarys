#' Set Up dracarys Variables
#'
#' Reads the `replay.json` file, which contains the DRAGEN command line,
#' parameters, version and inputs for the specific run. It then pre-creates
#' expected file names for DRAGEN output.
#'
#' @param x Path to `replay.json` file.
#'
#' @return A list with several elements pointing to DRAGEN output file names.
#'
#'
#' @examples
#' x <- system.file("extdata/COLO829-replay.json.gz", package = "dracarys")
#' setup_vars(x)
#'
#' @export
setup_vars <- function(x) {

  assertthat::assert_that(length(x) == 1, is.character(x))

  # - Get enable-sv, enable-variant-caller from replay file
  # - See if there is '--tumor-fastq' in the command line i.e. it's run in T/N mode.
  # - Output directory and file names for mapping/coverage metrics, SNV/SV VCFs etc.

  nm <- sub("-replay.json.*", "", basename(x))
  replay <- dracarys::read_replay(x)

  .replay_info <- function(replay) {

    d <- replay[["dragen_config"]]
    sv <- ("enable-sv" %in% d$name) && (d$value[d$name == "enable-sv"] == "true")
    vc <- ("enable-variant-caller" %in% d$name) && (d$value[d$name == "enable-variant-caller"] == "true")

    cl <- replay[["command_line"]]
    tn <- grepl("--tumor-fastq", cl)
    resdir <- dirname(x)

    tibble::tribble(
      ~tn, ~vc, ~sv, ~resdir,
      tn, vc, sv, resdir
    )
  }

  ri <- .replay_info(replay)

  suffix <- if (ri$tn) {
    c("_tumor", "_normal")
  } else {
    ""
  }

  list(
    replay_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}-replay.json")), lab = nm),
    fraglen_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.fragment_length_hist.csv")), lab = nm),
    map_met_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.mapping_metrics.csv")), lab = nm),
    ploidy_est_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.ploidy_estimation_metrics.csv")), lab = nm),
    time_met_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.time_metrics.csv")), lab = nm),

    cov_contig_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.wgs_contig_mean_cov{suffix}.csv")), lab = paste0(nm, suffix)),
    cov_met_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.wgs_coverage_metrics{suffix}.csv")), lab = paste0(nm, suffix)),
    cov_finehist_fn = list(fn = file.path(ri$resdir, glue::glue("{nm}.wgs_fine_hist{suffix}.csv")), lab = paste0(nm, suffix)),

    vc_met_fn = dplyr::if_else(ri$vc, file.path(ri$resdir, glue::glue("{nm}.vc_metrics.csv")), NA_character_),
    snv_fn = dplyr::if_else(ri$vc, file.path(ri$resdir, glue::glue("{nm}.hard-filtered.vcf.gz")), NA_character_),
    sv_fn = dplyr::if_else(ri$sv, file.path(ri$resdir, glue::glue("{nm}.sv.vcf.gz")), NA_character_),
    info = ri,
    label = nm
  )
}

