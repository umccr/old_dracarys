#' Read Replay File
#'
#' Reads the `replay.json` file, which contains the DRAGEN command line,
#' parameters, version and inputs for the specific run.
#'
#' @param x Path to `replay.json` file.
#'
#' @return A list with the following elements:
#'   - `command_line`: character. DRAGEN command line used.
#'   - `dragen_config`: tibble. Parameters used for the DRAGEN run.
#'   - `inputs`: tibble. Germline input FASTQ files (if run in Tumor/Normal mode).
#'   - `system`: list. DRAGEN version, node name, and kernel release.
#'
#'
#' @examples
#' x <- system.file("extdata/COLO829-replay.json.gz", package = "dracarys")
#' read_replay(x)
#'
#' @export
read_replay <- function(x) {

  res <- x %>%
    jsonlite::read_json(simplifyVector = TRUE) %>%
    purrr::map_if(is.data.frame, tibble::as_tibble)

  assertthat::assert_that(
    all(names(res) %in% c("command_line", "dragen_config", "inputs", "system")))

  res$dragen_config <- res$dragen_config %>%
    dplyr::mutate(value = ifelse(grepl("https", .data$value), "pre-signed URL maybe?", .data$value),
                  value = ifelse(.data$name == "vc-decoy-contigs", "lots of contigs", .data$value)) %>%
    dplyr::filter(!grepl("credentials", .data$name)) %>%
    dplyr::arrange(.data$name)

  res$system <- res$system %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(dplyr::everything())

  res$label <- sub("-replay.json.*", "", basename(x))

  res
}
