\name{AnnoEPIC}
\alias{AnnoEPIC}
\docType{data}
\title{
Illumina EPIC Annotation
}
\description{
This data set gives the CpG positions for in Illumina EPIC
array in build hs1
}
\usage{data("AnnoEPIC")}
\format{
A data frame with 850k rows and 2 columns:
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
<https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylationEPICmanifest.html>
}
\examples{
data("AnnoEPIC")
## view(AnnoEPIC) ...
}
\keyword{datasets}
