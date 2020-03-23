#' Read Fragment Length Hist File
#'
#' Reads the `fragment_length_hist.csv` file, which contains the
#' insert length distribution for each sample.
#'
#' @param x Path to `fragment_length_hist.csv` file.
#' @return A tibble with the following columns:
#'   - fragmentLength: estimated fragment length
#'   - count: number of reads with estimated fragment length
#'   - sample: name of sample
#'
#' @examples
#' x <- system.file("extdata/Colo829.fragment_length_hist.csv.gz", package = "dracarys")
#' (fl <- read_fragment_length_hist(x))
#'
#' @testexamples
#' expect_equal(colnames(fl), c("fragmentLength", "count", "sample"))
#'
#' @export
read_fragment_length_hist <- function(x) {

  d <- readr::read_lines(x)
  assertthat::assert_that(grepl("#Sample", d[1]))

  d %>%
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

#' Plot Fragment Length Distribution
#'
#' Plots the fragment length distributions as written in the `fragment_length_hist.csv` files.
#'
#' @param xs Character vector containing paths to `fragment_length_hist.csv` files.
#' @param min_count Minimum read count to be plotted (Default: 10).
#' @return A ggplot2 plot containing insert lengths on X axis and read counts
#'   on Y axis for each sample.
#'
#' @examples
#' x <- system.file("extdata/Colo829.fragment_length_hist.csv.gz", package = "dracarys")
#'
#' plot_fragment_length_hist(xs = x, min_count = 10)
#' @export
plot_fragment_length_hist <- function(xs, min_count = 10) {
  assertthat::assert_that(is.numeric(min_count), min_count >= 0, length(xs) > 0)
  d <- purrr::map(xs, read_fragment_length_hist) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(.data$count >= min_count)

  d %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$fragmentLength, y = .data$count, colour = sample)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Fragment Length Distribution") +
    ggplot2::xlab("Fragment Length (bp)") +
    ggplot2::ylab(glue::glue("Read Count (min: {min_count})")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = c(0.9, 0.9),
      legend.justification = c(1, 1),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(colour = "#2c3e50", size = 14, face = "bold"))
}
