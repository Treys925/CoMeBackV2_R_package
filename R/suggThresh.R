#' @title Calculate suggestive P-Value threshold
#'
#' @description Function \code{suggThresh} calculates a suggestive P-Value threshold using SVD (Singular Value Decomposition),
#' which estimates the number of linearly independent combinations of probes
#'
#' @param Data matrix of betas, preferably cell-type corrected, or M-values, where rows are observations and columns are array probes
#' @param Columns list of column i.e. probe names
#
#' @return numeric value representing a suggestive P-Value threshold
#' @examples
#' \dontrun{
#' suggThresh("Mdata_ProbesAreColumns_matrix")
#' }
#' @author Trey Smith, \email{treysm@@umich.edu}
#'
#' @export suggThresh
#' 
#' 

suggThresh <- function(Data,Columns){
  
  MeanCenteredAutoScaledData <- Data #center and scale data
  
    for(i in 1:length(Columns)){ #get means and SDs
      Entry <- Columns[i]
      Mean <- mean(Data[,Entry])
      Std <- sd(Data[,Entry])
      MeanCenteredAutoScaledData[,Entry] <- (MeanCenteredAutoScaledData[,Entry] - Mean)/Std
    }
  
  SVD <- svd(x=as.matrix(MeanCenteredAutoScaledData[,Columns])) #run singular value decomposition
  
  EigenValues <- SVD$d^2 / dim(SVD$u)[1] #get eigenvalues
  M <- length(EigenValues)
  L <- M - 1
  Var <- var(EigenValues)
  MEff <- M*(1.0-(L*Var/M^2))
  IntEVals <- as.numeric(EigenValues>=1.0)
  NonIntEVals <- EigenValues - floor(EigenValues)
  MEffLi <- ceiling(sum(IntEVals) + sum(NonIntEVals))
  
  thresh <- 0.05/MEffLi #calculate threshold using bonferroni correction
  
  return(thresh)
  
}