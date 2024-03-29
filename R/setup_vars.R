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
#' @examples
#' x <- system.file("extdata/COLO829-replay.json.gz", package = "dracarys")
#' (v <- setup_vars(x))
#'
#' @testexamples
#' expect_equal(class(v), "list")
#' expect_equal(length(v), 1)
#' expect_equal(length(v[[1]]), 14)
#' expect_equal(names(v[[1]]), c("name", "res_dir", "replay_fn", "fraglen_fn", "map_met_fn",
#' "ploidy_est_fn", "time_met_fn", "cov_contig_fn", "cov_met_fn",
#' "cov_finehist_fn", "vc_met_fn", "snv_fn", "sv_fn", "steps_run"))
#'
#' @export
setup_vars <- function(x) {

  assertthat::assert_that(length(x) == 1, is.character(x))

  # - Get enable-sv, enable-variant-caller from replay file
  # - See if there is '--tumor-fastq' in the command line i.e. it's run in T/N mode.
  # - Output directory and file names for mapping/coverage metrics, SNV/SV VCFs etc.

  resdir <- normalizePath(dirname(x))
  nm <- sub("-replay.json.*", "", basename(x))
  replay <- read_replay(x)

  .replay_info <- function(replay) {

    d <- replay[["dragen_config"]]
    sv <- ("enable-sv" %in% d$name) && (d$value[d$name == "enable-sv"] == "true")
    vc <- ("enable-variant-caller" %in% d$name) && (d$value[d$name == "enable-variant-caller"] == "true")

    cl <- replay[["command_line"]]
    tn <- grepl("--tumor-fastq", cl)

    list(
      tn = tn, vc = vc, sv = sv
    )
  }

  ri <- .replay_info(replay)

  suffix <- if (ri$tn) {
    c("_tumor", "_normal")
  } else {
    ""
  }

  l <- list(
    name = nm,
    res_dir = resdir,
    replay_fn = file.path(resdir, glue::glue("{nm}-replay.json")),
    fraglen_fn = file.path(resdir, glue::glue("{nm}.fragment_length_hist.csv")),
    map_met_fn = file.path(resdir, glue::glue("{nm}.mapping_metrics.csv")),
    ploidy_est_fn = file.path(resdir, glue::glue("{nm}.ploidy_estimation_metrics.csv")),
    time_met_fn = file.path(resdir, glue::glue("{nm}.time_metrics.csv")),

    cov_contig_fn = c(file.path(resdir, glue::glue("{nm}.wgs_contig_mean_cov{suffix}.csv"))),
    cov_met_fn = c(file.path(resdir, glue::glue("{nm}.wgs_coverage_metrics{suffix}.csv"))),
    cov_finehist_fn = c(file.path(resdir, glue::glue("{nm}.wgs_fine_hist{suffix}.csv"))),

    vc_met_fn = dplyr::if_else(ri$vc, file.path(resdir, glue::glue("{nm}.vc_metrics.csv")), NA_character_),
    snv_fn = dplyr::if_else(ri$vc, file.path(resdir, glue::glue("{nm}.hard-filtered.vcf.gz")), NA_character_),
    sv_fn = dplyr::if_else(ri$sv, file.path(resdir, glue::glue("{nm}.sv.vcf.gz")), NA_character_),

    steps_run = ri
  )

  res <- purrr::set_names(list(l), nm)
  res
}
