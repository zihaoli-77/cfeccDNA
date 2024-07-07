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
library(mgsub)
setwd("/public/home/hpc216085/eccDNA/sample/split.sorted")
rm(list = ls())


eccDNA_OLM <- function(bamfile){
param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                         hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                         isFirstMateRead = NA, isSecondMateRead = NA,
                                         isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                         isDuplicate = NA, isSupplementaryAlignment = NA),
                      what = 'seq',tag = 'SA')#设置读入bam文件的参数
differ <- readGAlignments( bamfile, param = param,use.names = T) #读bam文件
sa <- as.data.frame(differ@elementMetadata$SA)
#将SA的内容按照“，”分列
names(sa) <- 'SA'
sa <- separate(sa,SA,into = c('chr','pos','strand','cigar2','mapq','nm'),sep = ",")

a<-mcols(differ)
a.seqnames <- seqnames(differ)
a.start<-start(differ)
a.end<-end(differ)
a.strand<-strand(differ)
a.cigar <- cigar(differ)
df1<-(a)
a.id<-rownames(df1)
df1$p.chr <- a.seqnames
df1$p.start<-a.start
df1$p.end<-a.end
df1$p.strand<-a.strand
df1$p.cigar <- a.cigar
df1$p.signal <- NA
df1$p.signal <- ifelse(grepl("^[0-9]+M.*[0-9]+[SH]$|^[0-9][SH][0-9]+M", df1$p.cigar), "A","B")
df1$id<-a.id
df1$cigar2 <- sa$cigar2
df1$signal2 <- ifelse(grepl("^[0-9]+M.*[0-9]+[SH]$|^[0-9][SH][0-9]+M", df1$cigar2), "A","B")
df2<-df1[df1[, "p.cigar"] %>% str_detect("^[0-9]+[SHM][0-9]+[MSH]$"), ]
df2<-df2[df2[, "cigar2"] %>% str_detect("^[0-9]+[SHM][0-9]+[MSH]$"), ]

hg.38 <- BSgenome.Hsapiens.UCSC.hg38
# 获取所有非“chr1-22，X，Y”的seqnames
non_grch38_chr <- setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY"))
# 将seqnames列中非“chr1-22，X，Y”的行删除

#正链
df3 <- df2[df2[, "p.strand"] == "+", ] #正链


Mlong <- unlist(lapply(str_extract_all(df3$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))
Slong <- unlist(lapply(str_extract_all(df3$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))
range <- Mlong-Slong
range <- which(range > 3)#控制重叠片段的碱基数
df3 <- df3[range,]
df3$len<- unlist(lapply(str_extract_all(df3$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))-unlist(lapply(str_extract_all(df3$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))

qname1<- df3$id
start1<-df3$p.start
end1 <- df3$p.end
df3$end_1<-end1-df3$len
df3$start_1<-start1+df3$len
p.signal <- df3$p.signal
##
pA.qname <- qname1[p.signal == 'A']
df3$p.chr <- as.factor(df3$p.chr)
pA.chr <- df3$p.chr[p.signal == 'A']

pA.end <- end1[p.signal == 'A']
pA.end_1 <- df3$end_1[p.signal == 'A']
pB.qname <- qname1[p.signal == 'B']
pB.chr <- df3$p.chr[p.signal == 'B']
pB.start <- start1[p.signal == 'B']
pB.start_1 <- df3$start_1[p.signal == 'B']
gr1 <- GRanges(seqnames = pA.chr,ranges = IRanges(start = pA.end_1-3,end = pA.end_1))
gr2 <- GRanges(seqnames = pB.chr,ranges = IRanges(start = pB.start_1,end = pB.start_1+3))
pa.seq<-getSeq(hg.38,gr1)
pb.seq<-getSeq(hg.38,gr2)
pa.seq.df <- as.character(pa.seq)
pb.seq.df <- as.character(pb.seq)

pa.df<- as.data.frame(pA.qname)
pa.df$chr.1 <- pA.chr
pa.df$motif.1 <- pa.seq.df
names(pa.df) <- c("Qname","Chr.1","Motif.1")

pb.df <- as.data.frame(pB.qname)
pb.df$chr.1 <- pB.chr
pb.df$motif.1 <- pb.seq.df
names(pb.df) <- c("Qname","Chr.1","Motif.1")

one <- rbind(pa.df,pb.df)

one <- one[!grepl("N", one$Motif.1), ]
one<-as.data.frame(one$Motif.1)
names(one)<-'motif'
#负链
df4 <- df2[df2[, "p.strand"] == "-", ] #正链


Mlong2 <- unlist(lapply(str_extract_all(df4$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))
Slong2 <- unlist(lapply(str_extract_all(df4$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))
range2 <- Mlong2-Slong2
range2 <- which(range2 > 3)#控制重叠片段的碱基数
df4 <- df4[range2,]
df4$len<- unlist(lapply(str_extract_all(df4$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))-unlist(lapply(str_extract_all(df4$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))

qname2<- df4$id
start2<-df4$p.start
end2 <- df4$p.end
df4$end_2<-end2-df4$len
df4$start_2<-start2+df4$len
p.signal2 <- df4$p.signal
##
pA.qname2 <- qname2[p.signal2 == 'A']
df4$p.chr <- as.factor(df4$p.chr)
pA.chr2 <- df4$p.chr[p.signal2 == 'A']

pA.end2 <- end2[p.signal2 == 'A']
pA.end_2 <- df4$end_2[p.signal2 == 'A']
pB.qname2 <- qname2[p.signal2 == 'B']
pB.chr2 <- df4$p.chr[p.signal2 == 'B']
pB.start2 <- start2[p.signal2 == 'B']
pB.start_2 <- df4$start_2[p.signal2 == 'B']
gr3 <- GRanges(seqnames = pA.chr2,ranges = IRanges(start = pA.end_2-3,end = pA.end_2))
gr4 <- GRanges(seqnames = pB.chr2,ranges = IRanges(start = pB.start_2,end = pB.start_2+3))
pa.seq2<-getSeq(hg.38,gr3)
pb.seq2<-getSeq(hg.38,gr4)
pa.seq.df2 <- as.character(pa.seq2)
pb.seq.df2 <- as.character(pb.seq2)

pa.df2<- as.data.frame(pA.qname2)
pa.df2$chr.1 <- pA.chr2
pa.df2$motif.1 <- pa.seq.df2
names(pa.df2) <- c("Qname","Chr.1","Motif.1")

pb.df2 <- as.data.frame(pB.qname2)
pb.df2$chr.1 <- pB.chr2
pb.df2$motif.1 <- pb.seq.df2
names(pb.df2) <- c("Qname","Chr.1","Motif.1")

two <- rbind(pa.df2,pb.df2)

two <- two[!grepl("N", two$Motif.1), ]

two<-as.data.frame(mgsub(two$Motif.1,c('A','T','C','G'),c('T','A','G','C')))
names(two)<-'motif'

olm<- as.data.frame(prop.table(table(rbind(one,two))))

names(olm) <- c('motif','frequency')
number_motif<-"4"
base_pairs <- c("A", "T", "C", "G")
all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
names(empty.df)<-'motif'
OLM<-merge(empty.df,olm,by='motif',all.x=TRUE)
OLM[is.na(OLM)] <- 0

return(OLM)
}



# 获取所有符合特定文件名模式的文件列表
file_list_pri <- list.files(pattern = ".*quality\\.sorted\\.bam$")

# 生成输入文件列表
input_files <- lapply(file_list_pri, function(file) {
  # 提取文件名中的样本名称和重复次数
  sample_name <- sub("_split\\.freq2\\.quality\\.sorted\\.bam$", "", file)
  
  # 返回一个列表，包含两个输入文件和输出文件的路径
  list(input_file = file,
       output_file = paste0(sample_name, "_OLM", ".txt"))
})

for (i in seq_along(input_files)) {
  result <- eccDNA_OLM(input_files[[i]]$input_file)
  output_file <- file.path("/public/home/hpc216085/eccDNA/sample/split.sorted/4merOLM", input_files[[i]]$output_file)
  write.table(result, output_file, sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)
}