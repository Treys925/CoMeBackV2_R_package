#' @title Composite Methylation for Co-methylated Regions.
#'
#' @description Function \code{cmr_comp} estimates a composite mehtylation measure for co-methylated regions.
#'
#' @param cmrs=NULL is the list of CMRs construced with a call to the cmr function.  If not provided, CMRs are estimated with the methylation data.
#' @param Mdata matrix of DNA methylation values (Beta or M), rows are observations, columns are probes
#' @param cmethod the aggregarion method for the composite measure. One of "pca" (default), "mean" or "median"
#' @param retweights whether to return the weights derived from the 1st PC, default is FALSE.  If set to TRUE, the output is a list of 2 elemnts: the composite methylation matrix, and a list of per-cmr weights
#' @param verbose report progress, default is TRUE
#' @return matrix of CMR composite methylation values, one column per CMR.
#' @examples
#' \dontrun{
#' cmr_comp(CMRs,Mdata_ProbesAreColumns_matrix)
#' }
#' @author Trey Smith, \email{treysm@@umich.edu}
#' @importFrom stats median prcomp
#' @export cmr_comp
cmr_comp=function(cmrs,Mdata, cmethod=c("pca","mean","median"), retweights=F, verbose=T){
  cmethod <- match.arg(cmethod)
  nrMm=nrow(Mdata)
  ncMm=ncol(Mdata)
  nrTst=min(1000,nrMm)
  ncTst=min(1000,ncMm)
  if((nrMm>ncMm)&&(length(grep("cg",colnames(Mdata)[1:ncTst]))<1)&&(length(grep("cg",rownames(Mdata)[1:nrTst]))>0)) {
    if (verbose) print("Note: Probes should be in columns, so slowing down now to gently transpose the data for you.\nPlease be patient and if this operation is not over by the time you finish reading this, then it would appear that your computer is rather slow, so you may want to go for a stroll, or just chill and comeback later.")
    Mdata=t(Mdata) # for the chilled doof who do not RTFM
  }

  if(is.null(cmrs)){
    if (verbose) print("CMRs not provided, estimating now. This will take awhile..")
    cmrs=unlist(cmr(Mdata,verbose=verbose),recursive = F)
  } else {
    if(length(cmrs)<24) if (length(cmrs[[1]][[1]])!=1) cmrs=unlist(cmrs,recursive = F) # flat list, if passing chromosomes
  }

  ncmr=length(cmrs)

  cmr_frstPrb_nam=sapply(cmrs,function(x){x[[1]]}) # names for the CMRs - the first probe

  if (ncmr>1000) messs=" CMRs, this may take awhile."
  else messs=" CMRs."
  if (verbose) {
    print(paste0("Found ",ncmr,messs))
    print("Starting CMRs composite estimation.")}

  if(cmethod=="pca"){
    if(!retweights){
      res=as.matrix(do.call(cbind,lapply(cmrs,function(x){
      pc1=prcomp(Mdata[,x],retx = FALSE, center = TRUE, scale. = FALSE, rank. = 1)
      pcld=pc1$rotation
      cwght=pcld/sum(pcld)
      return(abs(Mdata[,x] %*% cwght))
        })))
    } else{ # also return the weights from 1st PC
      w1pc=lapply(cmrs,function(x){
        pc1=prcomp(Mdata[,x],retx = FALSE, center = TRUE, scale. = FALSE, rank. = 1)
        pcld=pc1$rotation
        cwght=pcld/sum(pcld)
        whichmaxw=which.max(abs(cwght))
        if(cwght[[whichmaxw]]<0){cwght= -cwght}
        return(cwght)
      })
      names(w1pc)=cmr_frstPrb_nam
      res=as.matrix(do.call(cbind,lapply(1:ncmr,function(x){
        return(Mdata[,cmrs[[x]]] %*% w1pc[[x]])
      })))
      }

    } else if(cmethod=="mean"){
      res=as.matrix(do.call(cbind,lapply(cmrs,function(x){
        return(apply(Mdata[,x],1,mean))
      })))
    } else{
      res=as.matrix(do.call(cbind,lapply(cmrs,function(x){
        return(apply(Mdata[,x],1,median))
      })))
    }

  rownames(res)=rownames(Mdata)
  colnames(res)=cmr_frstPrb_nam
  if((cmethod=="pca")&retweights){
    res=list(compMe=res,w1pc=w1pc)}
  if (verbose) print("Finished.")
  return(res)
}
