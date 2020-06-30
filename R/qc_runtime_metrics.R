#' Read Time Metrics File
#'
#' Reads the `time_metrics.csv` file, which contains
#' a breakdown of the run duration for each DRAGEN process.
#'
#' @param x Path to `time_metrics.csv` file.
#' @return tibble with the following columns:
#'   - Step: DRAGEN step
#'   - Time: time in HH:MM
#'
#' @examples
#' x <- system.file("extdata/COLO829.time_metrics.csv.gz", package = "dracarys")
#' (tm <- read_time_metrics(x))
#'
#' @testexamples
#' expect_equal(colnames(tm), c("Label", "Step", "Time"))
#'
#' @export
read_time_metrics <- function(x) {
  cn <- c("dummy1", "dummy2", "Step", "time_hrs", "time_sec")
  ct <- readr::cols(.default = "c", time_hrs = readr::col_time(format = "%T"), time_sec = "d")
  d <- readr::read_csv(x, col_names = cn, col_types = ct)
  assertthat::assert_that(d$dummy1[1] == "RUN TIME", is.na(d$dummy2[1]))
  assertthat::assert_that(inherits(d$time_hrs, "hms"))

  label <- sub(".time_metrics.csv.*", "", basename(x))

  d %>%
    dplyr::mutate(Step = tools::toTitleCase(sub("Time ", "", .data$Step)),
                  Label = label,
                  Time = substr(.data$time_hrs, 1, 5)) %>%
    dplyr::select(.data$Label, .data$Step, .data$Time)
}

#' Process Time Metrics File
#'
#' Processes the `time_metrics.csv` file.
#'
#' @param x Path to `time_metrics.csv` file. It can take multiple inputs (as a character vector).
#' @param id ID for each input (used to disambiguate files run on same samples). Default: index from 1 to length of `x`.
#' @return tibble with the following columns:
#'   - Step: DRAGEN step
#'   - Time: time in HH:MM
#'
#' @examples
#' x <- system.file("extdata/COLO829.time_metrics.csv.gz", package = "dracarys")
#' (tm <- time_metrics_process(c(x, x), id = c("run1", "run2")))
#'
#' @testexamples
#' expect_equal(nrow(tm), 2)
#'
#' @export
time_metrics_process <- function(x, id = seq_len(length(x))) {
  x %>%
    purrr::map(read_time_metrics) %>%
    purrr::set_names(id) %>%
    dplyr::bind_rows(.id = "ID") %>%
    tidyr::pivot_wider(id_cols = c("ID", "Label"), names_from = "Step", values_from = "Time") %>%
    dplyr::relocate(.data$`Total Runtime`, .after = .data$Label)
}
