\name{backgroundProbes}
\alias{backgroundProbes}
\docType{data}
\title{
Chromosome CpGs
}
\description{
List of CpGs by chromosome based on build hs1
}
\usage{data("backgroundProbes")}
\format{
  A list of 24
  \describe{
    \item{\code{Ch1}}{a list of CpG positions}
    \item{\code{Ch2}}{a list of CpG positions}
    ...
  }
}
\details{

}
\source{
<https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hs1.html>
}
\references{
Sergey Nurk \emph{et al} ., The complete sequence of a human genome.
\emph{Science} 376, 44-53(2022). DOI:10.1126/science.abj6987
}
\examples{
data(backgroundProbes)
## maybe str(backgroundProbes) ;  ...
}
\keyword{datasets}
