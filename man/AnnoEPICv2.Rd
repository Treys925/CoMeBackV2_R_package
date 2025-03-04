\name{AnnoEPICv2}
\alias{AnnoEPICv2}
\docType{data}
\title{
Illumina EPICv2 Annotation
}
\description{
This data set gives the CpG positions for in Illumina EPICv2
array in build hs1
}
\usage{data("AnnoEPICv2")}
\format{
A data frame with 900k rows and 2 columns:
  \describe{
    \item{\code{CHR}}{chromosome}
    \item{\code{MAPINFO}}{CpG position}
  }
}
\details{
Created using chain files from UCSC
}
\source{
<https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/>
}
\references{
<https://www.bioconductor.org/packages/release/data/annotation/html/EPICv2manifest.html>
}
\examples{
data("AnnoEPICv2")
## view(AnnoEPICv2) ...
}
\keyword{datasets}
