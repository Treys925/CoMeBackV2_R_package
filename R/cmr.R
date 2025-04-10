#' @title Construct Co-methylated Regions.
#'
#' @description Function \code{cmr} constructs co-methylated regions based on pair-wise positive correlations above a correlation threshold.
#' The threshold can be optionally set to vary with genomic CpG density, in which case it increases linearly from a low value for low-density CpG background to a high value for high density background.
#' Pearson Correlation is used by default, appropriate for cell-type corrected methylation betas.
#' Creating internal and external data for CoMeBackV2 Package
#'
#'
#' @param Mdata matrix of betas, preferably cell-type corrected, or M-values, where rows are observations and colums are array probes.
#' @param meds medians of probe betas, used if constraining adjacent probe maximum difference
#' @param Iarray which Illumina array platform, one of "450K" (default), "EPIC", or "EPICv2".
#' @param Build which Genome build, one of "hg38" (default) or "hs1". NOTE: for hs1 all files from /data must be loaded into the global environment.
#' @param cormethod a character string indicating which correlation is to be computed. One of "pearson" (default), "kendall", or "spearman".
#' @param corlo low threshold for correlation, when max distance between adjacent genomic CpGs is corlodst; optional - default is NULL
#' @param corhi high threshold for correlation, when max distance between adjacent genomic CpGs is corhidst; optional - default is NULL
#' @param corlodst the maximum distance between CpGs considered, defines low CpG density, default is 400bp
#' @param corhidst the minimum distance between CpGs considered, defines high CpG density, default is 1bp
#' @param maxprbdst Maximum bp distance between adjacent probes, default is 2000bp
#' @param maxdlvl the maximum difference between the medians of the adjacent probes
#' @param verbose whether to report progress for chromosomes, default is TRUE
#' @return List of chromosome-lists of co-methylated probe set character vectors.
#' @examples
#' \dontrun{
#' cmr("Mdata_ProbesAreColumns_matrix")
#' }
#' @author Trey Smith, \email{treysm@@umich.edu}
#' @importFrom stats cor
#' @export cmr
#' 
utils::globalVariables(c("Anno450k", "AnnoEPIC", "AnnoEPICv2", "backgroundProbes"))

cmr=function(Mdata,meds=NULL,Iarray=c("450K", "EPIC", "EPICv2"),Build = c("hg38", "hs1"),cormethod=c("pearson", "kendall", "spearman"), 
             corlo=NULL,corhi=corlo,corlodst=400,corhidst=1, maxprbdst=2000, maxdlvl=Inf, verbose=TRUE){ #
  # scan probes sequentially, get the number of CpGs between them, and based on that estimate correlation or not,
  # add next probe to current CMR if  correlation passes variable cutoff that depends on max CpG gap

  tryCatch(
            {

  cormethod <- match.arg(cormethod)
  Iarray <- match.arg(Iarray)
  Build <- match.arg(Build)

  if(is.null(corlo) && is.null(corhi)){
    corlo=corhi=ccutN_f(nrow(Mdata))
    print(paste0("No correlation cut-off specifified, using ad-hoc sample-size based cut-off ",corlo))
  }


  if (verbose) print("Getting CpG info, about to start estimation any day now.")
  
  if (Build == "hs1"){
    chr_seq_GpCpos=backgroundProbes
    if(Iarray=="450K") {EPIC_Manifest=Anno450k}
    else if (Iarray == "EPIC"){EPIC_Manifest=AnnoEPIC}
    else {EPIC_Manifest=AnnoEPICv2}
  }
    
    else{

    chr_seq_GpCpos=init_data$chr_seq_GpCpos
    if(Iarray=="450K") {EPIC_Manifest=init_data$I450K_Manifest}
    else if (Iarray == "EPIC"){EPIC_Manifest=init_data$EPIC_Manifest}
    else {EPIC_Manifest=init_data$EPICv2_Manifest}
    
  }

    EPIC_OK_prb=rownames(EPIC_Manifest)


    # take the given probes and sort them into chromosomes
    prb_MMat=colnames(Mdata)

    wkprbs=intersect(EPIC_OK_prb,prb_MMat)
    numprbs=length(wkprbs)

    EPIC_Manifest=EPIC_Manifest[wkprbs,]
    EPIC_Manifest$Name=wkprbs
    EPIC_Manifest=EPIC_Manifest[order(EPIC_Manifest$CHR,EPIC_Manifest$MAPINFO),]

    wkprbs_crd=EPIC_Manifest$MAPINFO
    names(wkprbs_crd)=rownames(EPIC_Manifest)

    chrom_nams=unique(EPIC_Manifest$CHR)
    chrom_nams=na.omit(chrom_nams)
    chrom_nams=chrom_nams[order(chrom_nams)]
    chroms=length(chrom_nams) # number of chromosomes present in the Mdata

  cmr_ac=vector(mode = "list", length = chroms)
  names(cmr_ac)=paste0("chr",chrom_nams)

  if (verbose) print(paste("Found ",numprbs," probes from ",chroms," out of 24 chromosomes",collapse = "" ))

  if((maxdlvl==1)|(is.null(meds))){

    for (i in 1:chroms){
      cni=chrom_nams[[i]]
      manchr=EPIC_Manifest[EPIC_Manifest$CHR==cni,c("Name","MAPINFO")]
      manchr=manchr[order(manchr$"MAPINFO"),]

      x=manchr$Name
      x_prbcrd=manchr$MAPINFO

      cni=paste0("chr",chrom_nams[[i]])
      cpgi=match(x_prbcrd,chr_seq_GpCpos[[cni]])
      cpgi_tf=!is.na(cpgi) # use only the probes that have genomic CpG match
      cpgi=cpgi[cpgi_tf]
      x=x[cpgi_tf]
      x_prbcrd=x_prbcrd[cpgi_tf]

      # for each pair of consecutive probes in this chromosome, are they correlated above the variable cut-off
      ccut_tf=sapply(2:length(cpgi), function(y){ # for each pair of consecutive probes in this cluster
        res=FALSE # the correlation is above cutoff
        # for this chromosomes' consecutive probes, all gaps between consecutive CpGs
        ym1=y-1
        dcpg=chr_seq_GpCpos[[cni]][(cpgi[[ym1]]+1):cpgi[[y]]]-chr_seq_GpCpos[[cni]][cpgi[[ym1]]:(cpgi[[y]]-1)]
        maxdcpg=max(dcpg)
        # if the maximum CpG gap is below the window, evaluate the correlation and check if cor is below threshold
        if((maxdcpg <= corlodst)&(abs(wkprbs_crd[[x[[(ym1)]]]]-wkprbs_crd[[x[[(y)]]]])<=maxprbdst)){
          corcut=corhi-(maxdcpg-corhidst)*(corhi-corlo)/(corlodst-corhidst) # the correlation cut-off for the pair of probes
          res=(cor(x=Mdata[,x[[(ym1)]]],y=Mdata[,x[[y]]],use = "pairwise.complete.obs", method = cormethod)>corcut)
        }
        return(res)
      })

      ccut_i=which(ccut_tf) # which pairs of adjacent probes are correlated above cutoff
      lpc=length(ccut_i)

      cmr1dcr=NULL
      if(lpc>0){
        cc=c(ccut_i[[1]]) # current cluster starts with first pair above cutoff
        cp=1 # pair count
        res=list()
        cct=0 # count of clusters
        if(lpc>1) {
          while(cp<lpc){
            cp1=cp+1
            if((ccut_i[[cp1]]-ccut_i[[cp]])==1){
              cc=c(cc,ccut_i[[cp1]]) # add next pair to the cluster
            } else { # next pair starts the next cluster
              cct=cct+1
              res[[cct]]=cc # add completed last cluster to the result
              cc=c(ccut_i[[cp1]]) # start new current cluster with the last pair
              if(cp1==lpc) res[[cct+1]]=cc # this is the last pair alone in the last cluster
            }
            cp=cp+1
          }
          # only one pair total
        } else {res[[1]]=cc}

        cmr1dcr=lapply(res,function(z){ # go from vector of pairs to vector of probes
          lz=length(z)
          if(lz==0) print(res)
          if(lz==1) {
            cl=x[z:(z+1)]
          } else{
            cl=x[z[[1]]:(z[[lz]]+1)]
          }
          return(cl)
        })
      }

      cmr_ac[[cni]]=cmr1dcr
      cmr_ac[[cni]]=cmr_ac[[cni]][sapply(cmr_ac[[cni]],length)>0]

      if (verbose) print(paste("Done chr",as.character(cni),"with",length(cmr_ac[[i]]),"cmrs"))
    } # chromosome loop
   } else {
    print(sprintf("Applying adjacent probe filter: max level difference %f",maxdlvl))
    names(meds)=colnames(Mdata)
    for (i in 1:chroms){
      cni=chrom_nams[[i]]
      manchr=EPIC_Manifest[EPIC_Manifest$CHR==cni,c("Name","MAPINFO")]
      manchr=manchr[order(manchr$"MAPINFO"),]

      x=manchr$Name
      x_prbcrd=manchr$MAPINFO

      cni=paste0("chr",chrom_nams[[i]])
      cpgi=match(x_prbcrd,chr_seq_GpCpos[[cni]])
      cpgi_tf=!is.na(cpgi) # use only the probes that have genomic CpG match
      cpgi=cpgi[cpgi_tf]
      x=x[cpgi_tf]
      x_prbcrd=x_prbcrd[cpgi_tf]

      # for each pair of consecutive probes in this chromosome, are they correlated above the variable cut-off
      ccut_tf=sapply(2:length(cpgi), function(y){ # for each pair of consecutive probes in this cluster
        res=FALSE # the correlation is above cutoff
        # for this chromosomes' consecutive probes, all gaps between consecutive CpGs
        ym1=y-1
        dcpg=chr_seq_GpCpos[[cni]][(cpgi[[ym1]]+1):cpgi[[y]]]-chr_seq_GpCpos[[cni]][cpgi[[ym1]]:(cpgi[[y]]-1)]
        maxdcpg=max(dcpg)
        # if the maximum CpG gap is below the window, evaluate the correlation and threshold, check if cor is below threshold
        if((maxdcpg <= corlodst)&(abs(wkprbs_crd[[x[[(ym1)]]]]-wkprbs_crd[[x[[(y)]]]])<=maxprbdst)){
          corcut=corhi-(maxdcpg-corhidst)*(corhi-corlo)/(corlodst-corhidst) # the correlation cut-off for the pair of probes
          res=((cor(x=Mdata[,x[[(ym1)]]],y=Mdata[,x[[y]]],use = "pairwise.complete.obs")>corcut)&(abs(meds[[x[[(ym1)]]]]-meds[[x[[y]]]])<maxdlvl))
        }
        return(res)
      })

      ccut_i=which(ccut_tf) # which pairs of adjacent probes are correlated above cutoff
      lpc=length(ccut_i)

      cmr1dcr=NULL
      if(lpc>0){
        cc=c(ccut_i[[1]]) # current cluster starts with first pair above cutoff
        cp=1 # pair count
        res=list()
        cct=0 # count of clusters
        if(lpc>1) {
          while(cp<lpc){
            cp1=cp+1
            if((ccut_i[[cp1]]-ccut_i[[cp]])==1){
              cc=c(cc,ccut_i[[cp1]]) # add next pair to the cluster
            } else { # next pair starts the next cluster
              cct=cct+1
              res[[cct]]=cc # add completed last cluster to the result
              cc=c(ccut_i[[cp1]]) # start new current cluster with the last pair
              if(cp1==lpc) res[[cct+1]]=cc # this is the last pair alone in the last cluster
            }
            cp=cp+1
          }
          # only one pair total
        } else {res[[1]]=cc}

        cmr1dcr=lapply(res,function(z){ # go from vector of pairs to vector of probes
          lz=length(z)
          if(lz==0) print(res)
          if(lz==1) {
            cl=x[z:(z+1)]
          } else{
            cl=x[z[[1]]:(z[[lz]]+1)]
          }
          return(cl)
        })
      }

      cmr_ac[[cni]]=cmr1dcr
      cmr_ac[[cni]]=cmr_ac[[cni]][sapply(cmr_ac[[cni]],length)>0]

      if (verbose) print(paste("Done",as.character(cni),"with",length(cmr_ac[[cni]]),"cmrs"))
    } # chromosome loop
  }

  return(cmr_ac)

  },
      error = function(cond){
            message("Not using a full set of CpGs")
            message("Here's the original error message:")
            message(conditionMessage(cond))
            # Choose a return value in case of error
            return(cmr_ac)
        }
  )

}
