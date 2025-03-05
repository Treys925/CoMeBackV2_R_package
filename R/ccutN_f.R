#' @title Correlation cut-off as function of sample size
#'
#' @description Empirical function based on simulation
#'
#' @param N is the sample size
#' @return the correlation cut-off value
#' @author Trey Smith, \email{treysm@@umich.edu}
ccutN_f=function(N){
  if (N<25 & N>0){res=0.5}
  else if(N <=0){res=NA}
  else {res=0.05+0.45*10^(-((N-25)^0.46)/21)} #max(0.1,0.38*10^(-((N-50)^0.48)/36))
  return(res)
}
