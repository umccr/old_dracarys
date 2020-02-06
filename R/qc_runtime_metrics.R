#' Read Time Metrics File
#'
#' Reads the `time_metrics.csv` file, which contains
#' a breakdown of the run duration for each DRAGEN process.
#'
#' @param x Path to `time_metrics.csv` file.
#' @return tibble with the following columns:
#'   - Step: DRAGEN step
#'   - Time: time in HH:MM:SS
#'
#' @examples
#' x <- system.file("extdata/COLO829.time_metrics.csv.gz", package = "dracarys")
#' read_time_metrics(x)
#' @export
read_time_metrics <- function(x) {
  cn <- c("dummy1", "dummy2", "Step", "time_hrs", "time_sec")
  ct <- readr::cols(.default = "c", time_sec = "d")
  d <- readr::read_csv(x, col_names = cn, col_types = ct)
  assertthat::assert_that(d$dummy1[1] == "RUN TIME", is.na(d$dummy2[1]))


  d %>%
    dplyr::mutate(Step = tools::toTitleCase(sub("Time ", "", .data$Step)),
                  Time = sub("\\..*", "", .data$time_hrs)) %>%
    dplyr::select(.data$Step, .data$Time)
}