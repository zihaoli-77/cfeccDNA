rm(list = ls())

library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(Rsamtools)
library(devtools)
#library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)

###
setwd("/your/path/to/bamfile") 
output_dir <- "/your/path/to/outfile"
###


eccDNA_EDM <- function(bamfile){
  
  param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                           hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                           isFirstMateRead = NA, isSecondMateRead = NA,
                                           isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                           isDuplicate = NA, isSupplementaryAlignment = NA),
                        what = "seq")
  primary <- readGAlignments(bamfile, param = param,use.names = T) 
  
  a<-mcols(primary)
  a.seqnames <- seqnames(primary)
  a.start<-start(primary)
  a.end<-end(primary)
  a.strand<-strand(primary)
  a.cigar <- cigar(primary)
  df1<-(a)
  a.id<-rownames(df1)
  df1$p.chr <- a.seqnames
  df1$p.start<-a.start
  df1$p.end<-a.end
  df1$p.strand<-a.strand
  df1$p.cigar <- a.cigar
  
  hg.38 <- BSgenome.Hsapiens.UCSC.hg38

  non_grch38_chr <- setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY"))
 
  
  df1 <- df1[!(df1$p.chr %in% non_grch38_chr),]
  
  start1<-df1$p.start
  
  end1 <- df1$p.end
  stand1 <- df1$p.strand
  stand1 <- as.factor(stand1)
  pA.chr <- df1$p.chr[stand1 == '-']
  pB.chr <- df1$p.chr[stand1 == '+']
  
  start <- df1$p.start[stand1=='+']
  end <- end1[stand1=='-']
  
  gr1 <- GRanges(seqnames = pA.chr,ranges = IRanges(start = end-3,end = end))
  gr2 <- GRanges(seqnames = pB.chr,ranges = IRanges(start = start,end = start+3))
  pa.seq<-getSeq(hg.38,gr1)
  pb.seq<-getSeq(hg.38,gr2)
  pa.seq.df <- as.character(pa.seq)
  pb.seq.df <- as.character(pb.seq)
  
  fu<-as.data.frame(pa.seq.df)
  zh<-as.data.frame(pb.seq.df)
  names(fu)<-'motif'
  names(zh)<-'motif'  
  com <- rbind(fu,zh)
  com <- as.data.frame(com[!grepl("N", com$motif), ])
  edm<- as.data.frame(prop.table(table(com)))
  names(edm) <- c('motif','frequency')
  number_motif<-"4"
  base_pairs <- c("A", "T", "C", "G")
  all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
  all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
  empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
  names(empty.df)<-'motif'
  EDM<-merge(empty.df,edm,by='motif',all.x=TRUE)
  EDM[is.na(EDM)] <- 0
  
  return(EDM)
  
}  

file_list_pri <- list.files(pattern = ".*quality\\.sorted\\.bam$")

input_files <- lapply(file_list_pri, function(file) {

  sample_name <- sub("_split\\.freq2\\.quality\\.sorted\\.bam$", "", file)
  
  list(input_file = file,
       output_file = paste0(sample_name, "_EDM", ".txt"))
})

for (i in seq_along(input_files)) {
  result <- eccDNA_EDM(input_files[[i]]$input_file)
  output_file <- file.path(output_dir, input_files[[i]]$output_file)
  write.table(result, output_file, sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
}