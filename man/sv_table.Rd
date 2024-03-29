% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{sv_table}
\alias{sv_table}
\title{Generate Table with Structural Variants}
\usage{
sv_table(manta)
}
\arguments{
\item{manta}{Path to \code{sv.vcf.gz} file or object of class 'sv'.}
}
\value{
A tibble with structural variants
}
\description{
Generate Table with Structural Variants
}
\examples{
manta <- system.file("extdata/COLO829.sv.vcf.gz", package = "dracarys")
sv_table(manta)

}
