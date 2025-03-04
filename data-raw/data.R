#' Creating internal and external data for CoMeBackV2 Package
#'
#' Use files from UCSC Genome Browser and Illumina Manifests
#' 
#'
#' @format ## `Anno450k`
#' A data frame with 450k rows and 2 columns:
#' \describe{
#'   \item{rownames}{CpG names}
#'   \item{CHR}{chromosome}
#'   \item{MAPINFO}{CpG position}
#'   ...
#' }
"Anno450k"
#' @source <https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/>
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
### 
#
## Code for installing necessary libraries is commented out
#
#
## First we need to create our list of background probes and annotation files for hg38

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#
#library(BSgenome.Hsapiens.UCSC.hg38)

# Go chromosome by chromosome and create lists of CpGs

chrs <- names(Hsapiens)[!grepl("alt|random|M|Un|fix", 
                               names(Hsapiens))] #filter out the "alternate", "random", and mitochondrial chromosomes

  CGs <- lapply(chrs, function(x) start(matchPattern("CG",
                                                   Hsapiens[[x]])))
    names(CGs) <- chrs


## Get and format Illumina Manifest files

#BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
chrEPIC2 <- sub('chr', '', IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations@listData$chr)
mapEPIC2 <- IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations@listData$pos
EPICv2_Manifest <- cbind(chrEPIC2, mapEPIC2)
EPICv2_Manifest <- data.frame(EPICv2_Manifest)
colnames(EPICv2_Manifest) <- c("CHR", "MAPINFO")
EPICv2_Manifest$probeID <- substr(IlluminaHumanMethylationEPICv2anno.20a1.hg38::Locations@rownames,1,10)
EPICv2_Manifest <- EPICv2_Manifest[!duplicated(EPICv2_Manifest$probeID), ]
row.names(EPICv2_Manifest) <- EPICv2_Manifest$probeID
EPICv2_Manifest <- EPICv2_Manifest[,-3]
EPICv2_Manifest$MAPINFO <- as.integer(EPICv2_Manifest$MAPINFO)

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
chrEPIC <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations@listData$chr
mapEPIC <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations@listData$pos
EPIC_Manifest <- cbind(chrEPIC, mapEPIC)
EPIC_Manifest <- data.frame(EPIC_Manifest)
colnames(EPIC_Manifest) <- c("CHR", "MAPINFO")
row.names(EPIC_Manifest) <- IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations@rownames
EPIC_Manifest$MAPINFO <- as.integer(EPIC_Manifest$MAPINFO)

#BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
chr450k <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations@listData$chr
map450k <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations@listData$pos
I450K_Manifest <- cbind(chr450k, map450k)
I450K_Manifest <- data.frame(I450K_Manifest)
colnames(I450K_Manifest) <- c("CHR", "MAPINFO")
row.names(I450K_Manifest) <- IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations@rownames
I450K_Manifest$MAPINFO <- as.integer(I450K_Manifest$MAPINFO)


### Chain hg19 to hg38 for EPIC and 450K
#
## Make GRanges objects for chaining
#
# library(GenomicRanges)

GR1 <- makeGRangesFromDataFrame(EPIC_Manifest,seqnames.field = "CHR", 
                                start.field = "MAPINFO", end.field = "MAPINFO", 
                                 na.rm = TRUE)

GR1@ranges@NAMES <- paste0(row.names(EPIC_Manifest), "_", EPIC_Manifest$CHR)
GR1@elementMetadata@rownames <- EPIC_Manifest$CHR

GR2 <- makeGRangesFromDataFrame(I450K_Manifest,seqnames.field = "CHR", 
                                start.field = "MAPINFO", end.field = "MAPINFO", 
                                na.rm = TRUE)

GR2@ranges@NAMES <- paste0(row.names(I450K_Manifest), "_", I450K_Manifest$CHR)
GR2@elementMetadata@rownames <- I450K_Manifest$CHR


#Get liftOver chain file from https://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/
# 
# BiocManager::install("liftOver")
# 
# library(liftOver)

ch <- import.chain("./hg19ToHg38.over.chain")
EPIC_hg38 <- unlist(liftOver(GR1, ch))

length(GR1) - length(EPIC_hg38) #Check number of probes lost
# 247

I450k_hg38 <- unlist(liftOver(GR2, ch))

length(GR2) - length(I450k_hg38)
# 168

# library(stringr)

AnnoEPIC = do.call("rbind", strsplit(EPIC_hg38@ranges@NAMES, split = "_"))
AnnoEPIC = data.frame(AnnoEPIC)
AnnoEPIC$MAPINFO <- EPIC_hg38@ranges@start
AnnoEPIC <- AnnoEPIC[!(duplicated(AnnoEPIC$X1)),]
AnnoEPIC <- tibble::column_to_rownames(AnnoEPIC, var = "X1")
colnames(AnnoEPIC)[1] <- "CHR"
AnnoEPIC$CHR <- str_replace(AnnoEPIC$CHR, "chr", "")


Anno450k = do.call("rbind", strsplit(I450k_hg38@ranges@NAMES, split = "_"))
Anno450k = data.frame(Anno450k)
Anno450k$MAPINFO <- I450k_hg38@ranges@start
Anno450k <- Anno450k[!(duplicated(Anno450k$X1)),]
Anno450k <- tibble::column_to_rownames(Anno450k, var = "X1")
colnames(Anno450k)[1] <- "CHR"
Anno450k$CHR <- str_replace(Anno450k$CHR, "chr", "")


init_data <- list(CGs, Anno450k, AnnoEPIC, EPICv2_Manifest)
names(init_data) <- c("chr_seq_GpCpos", "I450K_Manifest", "EPIC_Manifest",  "EPICv2_Manifest")

usethis::use_data(init_data, internal = TRUE)


### Create files for T2T Genome Build
# 
# BiocManager::install("BSgenome.Hsapiens.UCSC.hs1", force = TRUE)
#
# library(BSgenome.Hsapiens.UCSC.hs1)

chrs <- names(BSgenome.Hsapiens.UCSC.hs1::Hsapiens)[!grepl("alt|random|M|Un|fix", 
    names(BSgenome.Hsapiens.UCSC.hs1::Hsapiens))] #filter out the "alternate", "random", and mitochondrial chromosomes
  
  backgroundProbes <- lapply(chrs, function(x) start(matchPattern("CG", #get CpGs
                        BSgenome.Hsapiens.UCSC.hs1::Hsapiens[[x]])))
    
    names(backgroundProbes) <- chrs

    
### Chain hg38 to hs1 for all three arrays

AnnoEPIC$CHR <- paste0("chr" ,AnnoEPIC$CHR)

GR1 <- makeGRangesFromDataFrame(AnnoEPIC,seqnames.field = "CHR", 
                                start.field = "MAPINFO", end.field = "MAPINFO",
                                 na.rm = TRUE)

  GR1@ranges@NAMES <- paste0(row.names(AnnoEPIC), "_", AnnoEPIC$CHR)

    GR1@elementMetadata@rownames <- AnnoEPIC$CpG_chrm

    
Anno450k$CHR <- paste0("chr" ,Anno450k$CHR)


GR2 <- makeGRangesFromDataFrame(Anno450k,seqnames.field = "CHR", 
                                start.field = "MAPINFO", end.field = "MAPINFO",
                                na.rm = TRUE)

  GR2@ranges@NAMES <- paste0(row.names(Anno450k), "_", Anno450k$CHR)

    GR2@elementMetadata@rownames <- Anno450k$CHR


        
EPICv2_Manifest$CHR <- paste0("chr" ,EPICv2_Manifest$CHR)

  GR3 <- makeGRangesFromDataFrame(EPICv2_Manifest,seqnames.field = "CHR", 
                                start.field = "MAPINFO", end.field = "MAPINFO",
                                na.rm = TRUE)

    GR3@ranges@NAMES <- paste0(row.names(EPICv2_Manifest), "_", EPICv2_Manifest$CHR)

      GR3@elementMetadata@rownames <- EPICv2_Manifest$CHR


#Get liftOver chain file from https://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/

ch <- import.chain("./hg38ToHs1.over.chain")
EPIC_T2T <- unlist(liftOver(GR1, ch))

length(GR1) - length(EPIC_T2T) #Check number of probes lost
# 1048

I450k_T2T <- unlist(liftOver(GR2, ch))

length(GR2) - length(I450k_T2T)
# 733

EPICv2_T2T <- unlist(liftOver(GR3, ch))

length(GR3) - length(EPICv2_T2T)
# 1038


## create Annotations from GRanges objects

AnnoEPIC = do.call("rbind", strsplit(EPIC_T2T@ranges@NAMES, split = "_"))
AnnoEPIC = data.frame(AnnoEPIC)
AnnoEPIC$MAPINFO <- EPIC_T2T@ranges@start
AnnoEPIC <- AnnoEPIC[!(duplicated(AnnoEPIC$X1)),]
AnnoEPIC <- tibble::column_to_rownames(AnnoEPIC, var = "X1")
colnames(AnnoEPIC)[1] <- "CHR"
AnnoEPIC$CHR <- str_replace(AnnoEPIC$CHR, "chr", "")


Anno450k = do.call("rbind", strsplit(I450k_T2T@ranges@NAMES, split = "_"))
Anno450k = data.frame(Anno450k)
Anno450k$MAPINFO <- I450k_T2T@ranges@start
Anno450k <- Anno450k[!(duplicated(Anno450k$X1)),]
Anno450k <- tibble::column_to_rownames(Anno450k, var = "X1")
colnames(Anno450k)[1] <- "CHR"
Anno450k$CHR <- str_replace(Anno450k$CHR, "chr", "")


AnnoEPICv2 <- do.call("rbind", strsplit(EPICv2_T2T@ranges@NAMES, split = "_"))
AnnoEPICv2 <- data.frame(AnnoEPICv2)
AnnoEPICv2[3] <- EPICv2_T2T@ranges@start
AnnoEPICv2 <- AnnoEPICv2[!(duplicated(AnnoEPICv2$X1)),]
AnnoEPICv2 <- tibble::column_to_rownames(AnnoEPICv2, var = "X1")
colnames(AnnoEPICv2)[1:2] <- c("CHR", "MAPINFO")
AnnoEPICv2$CHR <- str_replace(AnnoEPICv2$CHR, "chr", "")


## Save CpG lists and annotations

save(backgroundProbes, file = "./HS1_CpGs.rda")

save(Anno450k, file = "./HS1_450k.rda")

save(AnnoEPIC, file = "./HS1_EPIC.rda")

save(AnnoEPICv2, file = "./HS1_EPICv2.rda")

