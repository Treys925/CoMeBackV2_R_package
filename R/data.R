#' Documentation for objects needed when using build hs1
#'
#' @format ## `Anno450k`
#' A data frame with 450k rows and 2 columns:
#' \describe{
#'   \item{rownames}{CpG names}
#'   \item{CHR}{chromosome}
#'   \item{MAPINFO}{CpG position}
#'   ...
#' }
#' @source <https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/>
"Anno450k"
#'
#' @format ## `AnnoEPIC`
#' A data frame with 850k rows and 2 columns:
#' \describe{
#'   \item{rownames}{CpG names}
#'   \item{CHR}{chromosome}
#'   \item{MAPINFO}{CpG position}
#'   ...
#' }
#' @source <https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/>
"AnnoEPIC"
#`
#' @format ## `AnnoEPICv2`
#' A data frame with 900k rows and 2 columns:
#' \describe{
#'   \item{rownames}{CpG names}
#'   \item{CHR}{chromosome}
#'   \item{MAPINFO}{CpG position}
#'   ...
#' }
#' @source <https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/>
"AnnoEPICv2"
#`
#' @format ## `backgroundProbes`
#' A list of 24:
#' \describe{
#'   \item{Chr1}{chromosome 1 CpG Positions}
#'   \item{Chr2}{chromosome 2 CpG Positions}
#'   \item{Chr3}{chromosome 3 CpG Positions}
#'   ...
#' }
#' @source <https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hs1.html>
"backgroundProbes"
#`
