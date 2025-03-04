\name{Anno450k}
\alias{Anno450k}
\docType{data}
\title{
Illumina 450k Annotation
}
\description{
This data set gives the CpG positions for in Illumina 450k
array in build hs1
}
\usage{data("Anno450k")}
\format{
A data frame with 450k rows and 2 columns:
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
<https://bioconductor.org/packages/release/data/annotation/html/IlluminaHumanMethylation450kmanifest.html>
}
\examples{
data("Anno450k")
## view(Anno450k) ...
}
\keyword{datasets}
