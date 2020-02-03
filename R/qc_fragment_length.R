#' Read fragment length hist file
#'
#' Reads the `fragment_length_hist.csv` file, which contains the
#' insert length distribution for each sample.
#'
#' @param x Path to `fragment_length_hist.csv` file.
#' @return A tibble with the following columns:
#'   - fragmentLength
#'
#' @examples
#' x <- system.file("extdata/Colo829.fragment_length_hist.csv.gz", package = "dracarys")
#' read_fragment_length_hist(x)
#' @export
read_fragment_length_hist <- function(x) {
  stopifnot(file.exists(x))
  readr::read_lines(x) %>%
    tibble::enframe() %>%
    dplyr::mutate(
      sample = ifelse(
        grepl("#Sample", .data$value),
        sub("#Sample: (.*)", "\\1", .data$value), NA_character_)) %>%
    tidyr::fill(.data$sample, .direction = "down") %>%
    dplyr::filter(!grepl("#Sample: |FragmentLength,Count", .data$value)) %>%
    tidyr::separate(.data$value, c("fragmentLength", "count"), convert = TRUE) %>%
    dplyr::select(-.data$name)
}

#' Plot insert length distribution
#'
#' Plots the insert length distributions as written in the `fragment_length_hist.csv` file.
#'
#' @param x Path to `fragment_length_hist.csv` file.
#' @param min_count Minimum read count to be plotted (Default: 10).
#' @return A ggplot2 plot containing insert lengths on X axis and read counts
#'   on Y axis for each sample.
#'
#' @examples
#' x <- system.file("extdata/Colo829.fragment_length_hist.csv.gz", package = "dracarys")
#' p <- plot_fragment_length_hist(x, min_count = 10)
#' @export
plot_fragment_length_hist <- function(x, min_count = 10) {
  stopifnot(is.numeric(min_count), min_count >= 0)
  read_fragment_length_hist(x) %>%
    dplyr::filter(.data$count >= min_count) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$fragmentLength, y = .data$count, colour = sample)) +
    ggplot2::geom_line() +
    ggplot2::xlab("Insert Length (bp)") +
    ggplot2::ylab(glue::glue("Read Count (min: {min_count})")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = c(0.9, 0.9),
                   legend.justification = c(1, 1),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank())
}
