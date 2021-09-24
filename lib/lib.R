library(data.table)
library(dtplyr)
library(plyr)
library(dplyr)

library(zoo)
library(mgcv)
library(visreg)
library(fitdistrplus)
#install.packages("itsadug")
library(itsadug)

##library(mgcViz)
#devtools::install_github('gavinsimpson/gratia')
#library(gratia)

library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(ggpubr)

library(ggExtra)

library(gtable)
library(gridExtra)
library(grid)
library(ggsci)
library(forcats)
library(RColorBrewer)
#install.packages("tidyverse")
library(forcats)
#install.packages("corrplot")
library(corrplot)
#install.packages("reshape")
library(reshape)
library(pheatmap)

library(devtools)



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
#
# BiocManager::install(c("GenomicRanges", "SomaticSignatures", "GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene", "BSgenome.Hsapiens.UCSC.hg19", "GRanges"))
library(biovizBase)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("karyoploteR")
library(karyoploteR)

get_BP_DT <- function(ISV_DT){
  #split SV in BP
  l <- list(ISV_DT[,c("chrom1", "start1", "end1", "strand1", "svclass", "donor_id", "cancer_id", "sv_id", "tumor_grade", "tumor_stage")],
            ISV_DT[,c("chrom2", "start2", "end2", "strand2", "svclass", "donor_id", "cancer_id", "sv_id", "tumor_grade", "tumor_stage")] )
  setattr(l, 'names', c("BP1", "BP2"))
  l$BP1[,BP_order := 1]
  l$BP2[,BP_order := 2]

  #Data Table of BP
  BP_DT <- rbindlist(l, use.names=FALSE, fill=FALSE, idcol=FALSE)
  colnames(BP_DT) <-c("chrom", "start", "end", "strand", "sv_class", "donor_id", "cancer_type", "sv_id", "tumor_grade", "tumor_stage", "BP_order")
  keycols <- c("chrom", "start", "end")
  setkeyv(BP_DT, cols=keycols , physical = TRUE) # haskey(BP_DT)

  return(BP_DT)
}

get_signal_BP_DT <- function(BP_DT){
  # Generating a regression function of BPRscore per BP getting a fit.RDscore value for each BP location
  # Compute the peaks and marking the center.peak as TRUE
  # browser()
  chr <- unique(as.character(BP_DT[, chrom]))
  List_BP_DT <- lapply(chr, get_sv_signal, target_DT= BP_DT, span = 0.25)
  signal_BP_DT <- rbindlist(List_BP_DT, use.names=FALSE, fill=FALSE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value
  signal_BP_DT [, BP.ps := (nlogRDf + abs(-log10(max(chr_arm_lengths[chr=="chr1",]$p_end))))/abs(-log10(max(chr_arm_lengths[chr=="chr1",]$p_end)))]

  chr <- unique(as.character(signal_BP_DT[chrom != "chrY" , chrom]))
  List_BP_DT <- lapply(chr, get_sv_signal2, target_DT=signal_BP_DT, window = 50, span = 0.25)
  signal_BP_DT <- rbindlist(List_BP_DT, use.names=FALSE, fill=FALSE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value

  #adding the chr arm info
  signal_BP_DT <- signal_BP_DT[ end <= p_end, arm := "p"]
  signal_BP_DT <- signal_BP_DT[ start >= p_end, arm := "q"]
  signal_BP_DT <- signal_BP_DT[ start <= p_end & end >= p_end , arm := "b"]

  signal_BP_DT$chrom <- factor(signal_BP_DT$chrom, levels = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'))
  keycols <- c("chrom", "start", "end")
  setkeyv(signal_BP_DT, cols=keycols , physical = TRUE)
  #key(signal_BP_DT)

  ## reordering the global_index and adding space between chromosome arms for graphical aims
  signal_BP_DT[, global_indx := as.numeric(row.names(signal_BP_DT))]
  signal_BP_DT[, global_indx := suppressWarnings(ifelse(substring(chrom, 4,nchar(as.character(chrom))) == "X", global_indx+22*5,
                                                        ifelse(substring(chrom, 4,nchar(as.character(chrom))) == "Y", global_indx+23*5,
                                                               global_indx+(as.numeric(substring(chrom, 4,nchar(as.character(chrom))))-1)*5)))]

  signal_BP_DT <- signal_BP_DT[(is.na(strand) | strand != "/") & chrom != "chrY" & start <= q_end, ]

  return(signal_BP_DT)
}

get_GC_covariate <- function( target_DT) {
  GC <-  as.data.table( apply(target_DT[,c("chrom", "start", "end")], 1, function(row) GCcontent(Hsapiens, GRanges(as.character(row[1]), IRanges(as.numeric(row[2])-50, as.numeric(row[3])+50)))) )
  colnames(GC) <- c("GC")
  BP_GC <- cbind( target_DT, GC )
  return(BP_GC)
}

get_LAD_covariate <- function(LAD, target_DT) {
  BP_LAD <- foverlaps(target_DT, LAD, type="any", mult = "last", nomatch=NA)
  BP_LAD <- BP_LAD[,-c(2,3)]
  BP_LAD <- BP_LAD[, c(1,3:dim(BP_LAD)[2], 2) , with=FALSE] #
  colnames(BP_LAD)[c(2,3)] <- c("start", "end")
  BP_LAD <-  BP_LAD[is.na(get('LAD')) , LAD := 0.65 ]
  BP_LAD[ , LAD.factor:= 0]
  BP_LAD[ LAD > 0.65, LAD.factor:= 1]
  BP_LAD$LAD.factor <- factor(BP_LAD$LAD.factor, levels = c(0, 1))
  BP_LAD$LAD.factor <- as.factor(BP_LAD$LAD.factor)
  return(BP_LAD)
}

get_ChromMark_covariate <- function(ChromHMM_Mark, target_DT) {
  BP_ChromMark <- foverlaps(target_DT, ChromHMM_Mark, type="any", mult = "last", nomatch=NA)
  BP_ChromMark <- BP_ChromMark[,-c(2,3)]
  BP_ChromMark <- BP_ChromMark[, c(1,3:dim(BP_ChromMark)[2], 2) , with=FALSE] #
  colnames(BP_ChromMark)[c(2,3)] <- c("start", "end")
  BP_ChromMark <-  BP_ChromMark[is.na(get('ChromMark')) , ChromMark := "-" ]
  BP_ChromMark$ChromMark <- factor(BP_ChromMark$ChromMark, levels = c("-", '15_Quies', '14_ReprPCWk', '13_ReprPC', '12_EnhBiv', '11_BivFlnk', '10_TssBiv', '9_Het', '8_ZNF/Rpts', '7_Enh', '6_EnhG', '5_TxWk', '4_Tx', '3_TxFlnk', '2_TssAFlnk', '1_TssA'))
  BP_ChromMark$ChromMark <- as.factor(BP_ChromMark$ChromMark)
  return(BP_ChromMark)
}

get_RepClass_covariate <- function(Rep.Class, target_DT) {
  BP_Rep.Class <- foverlaps(target_DT, Rep.Class[,c(6,7,8,12)], type="any", mult = "last", nomatch=NA)
  BP_Rep.Class <- BP_Rep.Class[,-c(2,3)]
  BP_Rep.Class <- BP_Rep.Class[, c(1,3:dim(BP_Rep.Class)[2], 2) , with=FALSE] #
  colnames(BP_Rep.Class)[c(2,3)] <- c("start", "end")
  BP_Rep.Class <-  BP_Rep.Class[is.na(get('repClass')) , repClass := "-" ]
  BP_Rep.Class$repClass <- factor(BP_Rep.Class$repClass, levels = c("-", as.vector(unique(BP_Rep.Class$repClass))[as.vector(unique(BP_Rep.Class$repClass))!= "-"])) #c("-", "Simple_repeat", "SINE", "LTR", "LINE", "DNA", "Low_complexity", "RC", "Other", "tRNA", "Satellite", "rRNA", "Unknown", "srpRNA", "snRNA", "RNA")
  BP_Rep.Class$repClass <- as.factor(BP_Rep.Class$repClass)
  return(BP_Rep.Class)
}

get_TADsegment_covariate <- function(TAD_segment, target_DT) {
  BP_TAD.segment <- foverlaps(target_DT, TAD_segment, type="any", mult = "last", nomatch=NA)
  colnames(BP_TAD.segment)[c(2,3,4,5,6)] <- c("TAD_segment.start", "TAD_segment.end", "TAD_segment.class", "start", "end")
  BP_TAD.segment <- BP_TAD.segment[, c(1,5:dim(BP_TAD.segment)[2], 2:4) , with=FALSE] #
  return(BP_TAD.segment)
}

get_RT_covariate <- function(RT, target_DT) {
  BP_RT <- foverlaps(target_DT, RT[,c(1,2,3,4,7)], type="any", mult = "last", nomatch=NA)
  BP_RT <- BP_RT[,-c(2,3)]
  BP_RT <- BP_RT[, c(1,4:dim(BP_RT)[2], 2,3) , with=FALSE] #
  colnames(BP_RT)[c(2,3)] <- c("start", "end")
  BP_RT <-  BP_RT[is.na(get('RT.mean')), c("RT.mean", "RTodd"):= .(0,1) ]
  BP_RT <-  BP_RT[is.na(get('RT.mean')), RTodd := 1 ]
  return(BP_RT)
}


get_Gene_density_covariate <- function(Gene_density, target_DT) {
  BP_Gene_density <- foverlaps(target_DT, Gene_density[,c(1,2,3,6)], type="any", mult = "last", nomatch=NA)
  BP_Gene_density <- BP_Gene_density[,-c(2,3)]
  BP_Gene_density <- BP_Gene_density[, c(1,3:dim(BP_Gene_density)[2], 2) , with=FALSE] #
  colnames(BP_Gene_density)[c(2,3)] <- c("start", "end")
  return(BP_Gene_density)
}


get_TAD_recurr_covariate <- function(TAD_recurr, target_DT) {
  BP_TAD_recurr <- foverlaps(target_DT, TAD_recurr, type="any", mult = "all", nomatch=NA)
  colnames(BP_TAD_recurr)[c(2,3,4,5,6,7)] <- c("TAD.start", "TAD.end", "TAD.data_set","TAD.recurrence", "start", "end")
  BP_TAD_recurr <- BP_TAD_recurr[, c(1,6:dim(BP_TAD_recurr)[2], 2:5) , with=FALSE] #
  BP_TAD_recurr[,TAD.recurrence:= TAD.recurrence/36]
  BP_TAD_recurr[, BP.TAD.recurrence := mean(TAD.recurrence), by = .(chrom, start, end, strand, BP_order, sv_class, donor_id, cancer_type, sv_id)]
  BP_TAD_recurr <- BP_TAD_recurr[ !duplicated(BP_TAD_recurr[ , .(chrom, start, end, strand, BP_order, sv_class, donor_id, cancer_type, sv_id, tumor_grade, tumor_stage, BP.TAD.recurrence)]), ]
  return(BP_TAD_recurr)
}

get_TAD.B_recurr_covariate <- function(TAD.B_recurr, target_DT) {
  BP_TADB_recurr <- foverlaps(target_DT, TAD.B_recurr, type="any", mult = "all", nomatch=NA) #[Data_set=="Boundaries_T470_raw-rep1_TADs.txt"]
  colnames(BP_TADB_recurr)[c(2,3,4,5,6,7)] <- c("TAD_boundary.start", "TAD_boundary.end", "TAD_boundary.data_set", "TAD_boundary.recurrence", "start", "end")
  BP_TADB_recurr <- BP_TADB_recurr[, c(1,6:dim(BP_TADB_recurr)[2], 2:5) , with=FALSE] #
  BP_TADB_recurr[, TAD_boundary.recurrence := TAD_boundary.recurrence/36]
  BP_TADB_recurr[, BP.TAD_boundary.recurrence := mean(TAD_boundary.recurrence), by = .(chrom, start, end, strand, BP_order, sv_class, donor_id, cancer_type, sv_id)]
  BP_TADB_recurr <- BP_TADB_recurr[ !duplicated(BP_TADB_recurr[ , .(chrom, start, end, strand, BP_order, sv_class, donor_id, cancer_type, sv_id, tumor_grade, tumor_stage, BP.TAD_boundary.recurrence)]), ]
  BP_TADB_recurr[is.na(get('BP.TAD_boundary.recurrence')) , BP.TAD_boundary.recurrence := 0 ]
  return(BP_TADB_recurr)
}

get_FS_covariate <- function(FS, target_DT) {
  BP_FS <- foverlaps(target_DT, FS, type="any", mult = "last", nomatch=NA)
  BP_FS <- BP_FS[,-c(2,3)]
  BP_FS <- BP_FS[, c(1,3:dim(BP_FS)[2], 2) , with=FALSE] #
  colnames(BP_FS)[c(2,3)] <- c("start", "end")
  BP_FS[ !is.na(get('FS')), FS := "1"]
  BP_FS[ is.na(get('FS')), FS := "0"]
  BP_FS$FS <- factor(BP_FS$FS, levels = c("0", "1"))
  BP_FS$FS <- as.factor(BP_FS$FS)
  return(BP_FS)
}



Gene.Density <- function( G.window =1e6){
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- genes(txdb)
kp <- plotKaryotype(genome="hg19")
kp <- kpPlotDensity(kp, all.genes, window.size= G.window)
Gene_density <- cbind(as.data.table(kp$latest.plot$computed.values$windows ), as.data.table( kp$latest.plot$computed.values$density) )
colnames(Gene_density)[c(1,6)] <- c("chrom", "Gene.density")
setkey(Gene_density, chrom, start, end, physical = TRUE)
return(Gene_density)
}

LAD <- function(){
  LAD <- fread( file= "./Input/gennome_annotation/LAD.bed", sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, col.names = c("chrom", "start", "end", "LAD"), key = c("chrom", "start", "end") )
return(LAD)
}


Rep.Time <- function(){
RT <- as.data.table(read.delim("./Input/gennome_annotation/RT.bedGraph", header=FALSE))
colnames(RT) <- c("chrom", "start", "end", "RT.mean", "early", "late")
RT[,RTodd := early/late]
setkey(RT, chrom, start, end, physical = TRUE)
return(RT)
}

TAD.Segment <- function(cohort_name){
 # TAD_segment <- as.data.table(read.delim("./Input/gennome_annotation/TAD_segment_ChromMark_class.txt", sep = "", quote = "\"'", header=TRUE))
  TAD_segment <- as.data.table(read.delim(paste0("./Input/cancer_organ/",cohort_name,"/TAD_segment_",cohort_name,"_ChromMark_class.txt") , sep = "", quote = "\"'", header=TRUE))
  TAD_segment <- unique(TAD_segment[ tissue == cohort_name, c("chr", "start", "end", "TAD.class")]) # == name_map[tisue==f_name]$chromHMM
  colnames(TAD_segment) <- c("chrom", "start", "end", "TAD_segment.class")
  TAD_segment$end <- as.integer(TAD_segment$end)
  TAD_segment$start <- as.integer(TAD_segment$start)
  TAD_segment$chrom <- as.factor(TAD_segment$chrom)
  setkey(TAD_segment, chrom, start, end, physical = TRUE)
return(TAD_segment)
}

TAD.Recurrence <- function(){
  TAD_recurr <- fread( file="./Input/gennome_annotation/TADS_count_70overlap.bed", sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, col.names = c("chrom", "start", "end", "Data_set", "recurr"), key = c("chrom", "start", "end") )
return(TAD_recurr)
}

TADB.Recurrence <- function(){
  TAD.B_recurr <- fread( file="./Input/gennome_annotation/TADSBoundaries_count_overlap.bed", sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, col.names = c("chrom", "start", "end", "Data_set", "recurr"), key = c("chrom", "start", "end") )
return(TAD.B_recurr)
}

ChromHMM.Mark <- function(cohort_name){
  HMM_ChromMark <- as.data.table( read.delim(paste0("./Input/gennome_annotation/ChromHMM/",cohort_name,".bed"), header=FALSE))
  #HMM_ChromMark <- as.data.table( read.delim(paste0("./Input/cancer_organ/",cohort_name,"/ChromHMM.bed"), header=FALSE))
  colnames(HMM_ChromMark) <- c("chrom", "start", "end", "ChromMark")
  setkey(HMM_ChromMark, chrom, start, end, physical = TRUE)
return(HMM_ChromMark)
}

Repeat.Class <- function(){
  Rep.Class <- as.data.table( read.delim("./Input/gennome_annotation/Repeats.bed"))
  colnames(Rep.Class)[c(6,7,8)] <- c("chrom", "start", "end")
  setkey(Rep.Class, chrom, start, end, physical = TRUE)
  Rep.Class[,repClass := gsub("\\?", '', repClass)]
return(Rep.Class)
}

Fragile.Site <- function(){
  FS <- fread( file="./Input/gennome_annotation/FS.bed", sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "FS"), colClasses=c("character", "integer", "integer", "factor"), key = c("chrom", "start", "end") )
return(FS)
}

get_GAM <- function( target_DT) {
  gam_model <- gam(curve.o ~ FS +  LAD.factor + repClass + ChromMark +
                    s( RT.mean,  bs = "cr", by = LAD.factor) +
                      s( GC,  bs = "cr", by = LAD.factor) +
                      s( BP.TAD.recurrence, bs = "cr", by = LAD.factor) +
                      s( BP.TAD_boundary.recurrence, bs = "cr", by = LAD.factor) +
                      s( Gene.density, bs = "cr", by = LAD.factor) +
                      te( Gene.density, RT.mean, bs = c("cr", "cr"), by = TAD_segment.class) +
                      te( Gene.density, BP.TAD.recurrence, bs = c("cr", "cr"), by = TAD_segment.class) +
                      te( RT.mean, BP.TAD.recurrence, bs = c("cr", "cr"), by = TAD_segment.class)
                    , data = target_DT, na.action=na.exclude,
                    family = Gamma(link = "log"), method="REML" )
  return(gam_model)
}

chr_arm_lengths <- as.data.table(readRDS("./Input/gennome_annotation/chr_arm_lengths.rds"))

#============================================
## funtion for Generating a regression function of RDscore per BP getting a fit.RDscore value for each BP locat
get_sv_signal <- function(ch, target_DT, span = 0.2) {
  chr_dt <- target_DT[chrom == ch] #BP_DT
  # chr_dt[,chrom := paste("chr",chrom,sep ="")]
  setkey(chr_dt, cols=chrom , physical = TRUE)
  setkey(chr_arm_lengths, cols=chr , physical = TRUE)
  chr_dt <- chr_dt[chr_arm_lengths, nomatch=0] #, on =.(chrom == chr, start <= p_end),
  l <- list ( p_arm = chr_dt[start <= p_end],
              cen = list(chr_dt[1,]$chrom, chr_dt[1,]$p_end, chr_dt[1,]$p_end, "/", "CENTROMER","CENTROMER", "CENTROMER", "CENTROMER", "CENTROMER", "CENTROMER", "/", chr_dt[1,]$p_start, chr_dt[1,]$p_end, chr_dt[1,]$q_start, chr_dt[1,]$q_end), #
              p_arm = chr_dt[start > p_end],
              tel = list(chr_dt[1,]$chrom, chr_dt[1,]$q_end, chr_dt[1,]$q_end, "/", "TELOMER", "TELOMER", "TELOMER", "TELOMER", "TELOMER", "TELOMER", "/", chr_dt[1,]$p_start, chr_dt[1,]$p_end, chr_dt[1,]$q_start, chr_dt[1,]$q_end) #
  )
  chr_dt <- rbindlist(l, use.names=FALSE, fill=FALSE, idcol=FALSE)

  chr_dt[,indx := as.numeric(row.names(chr_dt))]

  chr_dt[, F_RD := start - data.table::shift(start)] #, by = list(arm_end )
  chr_dt[is.na(F_RD), F_RD := start - p_start] # 0]#

  setorder(chr_dt, -end, -indx)
  chr_dt[, R_RD := abs(end -  data.table::shift(end) )]
  chr_dt[is.na(R_RD), R_RD := abs(end - q_end)] # 0]#

  setorder(chr_dt, start, indx)
  chr_dt[, BP_NR := (F_RD + R_RD)/2]


  chr_dt[start <= p_end, nlogRDf := -log10(BP_NR+1)] # 1-BP_NR] #-log10((BP_NR+1)/(p_end-p_start+1))]
  chr_dt[start  > p_end, nlogRDf := -log10(BP_NR+1)] # 1-BP_NR] #-log10((BP_NR+1)/(q_end-q_start+1))]

  #  w = ifelse(dim(chr_dt)[1] < 50 , dim(chr_dt)[1], 50) # the w (number of BP) should be a setup variable
  #  signal <- test(w, span, chr_dt$indx, chr_dt$nlogRDf, chr_dt[sv_class=="CENTROMER"]$indx)
  #chr_dt[as.data.table(signal["indx"]), on=.(indx), center.peak:= TRUE ]
  #  chr_dt <- cbind(chr_dt, as.data.table(signal["curve.o"]))

  #  return(chr_dt)
  #tail(chr_dt)
}
#--------------------------------------------
get_sv_signal2 <- function(ch, target_DT, window = 25, span = 0.2) {
  chr_dt <- target_DT[chrom == ch] #BP_DT
  setkey(chr_dt, cols=chrom , physical = TRUE)
  w = ifelse(dim(chr_dt)[1] < window , dim(chr_dt)[1], window ) # the w (number of BP) should be a setup variable
  signal <- test(w, span, chr_dt$indx, chr_dt$BP.ps, chr_dt[sv_class=="CENTROMER"]$indx)
  #chr_dt[as.data.table(signal["indx"]), on=.(indx), center.peak:= TRUE ]
  chr_dt <- cbind(chr_dt, as.data.table(signal["curve.o"]))

}

# Functions
test <- function(w, span, x, y, centromer) {
  peaks <- argmax(x, y, w=w, span=span)

  list(indx = x[peaks$i], curve.o = peaks$y.hat) # peak no in use
}

#----PEAK---
argmax <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  # print(loess(y ~ x, ...))
  y.smooth <- loess(y ~ x, ...)$fitted # family = c("gaussian")
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}

#---------------------------------------------


get_peaks_signal <- function(ch, target_DT, window = 25, span = 0.2) {

  chr_dt <- target_DT[chrom == ch] #BP_DT
  w = ifelse(dim(chr_dt)[1] < window , dim(chr_dt)[1], window) # the w (number of BP) should be a setup variable
  #  span = 0.2                                          # or 0.2  the span should be a setup variable
  signal <- get_peaks(w, span, chr_dt) #, chr_dt[sv_class=="CENTROMER"]$indx

  chr_dt[as.data.table(signal["indx.max"]), on=.(global_indx == indx.max), peak.top:= TRUE ]
  chr_dt[as.data.table(signal["indx.min"]), on=.(global_indx == indx.min), valley.bottom:= TRUE ]
  chr_dt <- cbind(chr_dt, as.data.table(signal["curve.c"]))

  #  return(chr_dt)
  #tail(chr_dt)

}

#----
get_peaks <- function(w, span, chr_dt) { #, centromer
  peaks <- significant_peaks(chr_dt, w=w, span=span)
  x <- chr_dt$global_indx
  list(indx.min = x[peaks$i.min], indx.max = x[peaks$i.max], curve.c = peaks$curve.c )
}

## Area peaks------
#funtion to compute the peak based on the y_c
significant_peaks <- function(chr_dt, w=25, ...) {
  x <- chr_dt$global_indx
  y <- chr_dt$y_c
  require(zoo)
  n <- length(y)
  y_c.smooth <- loess(y ~ x, ...)$fitted # family = c("gaussian")
  y.max <- rollapply(zoo(y_c.smooth), 2*w+1, max, align="center")
  find.max <- y.max - y_c.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(find.max <= 0) + w

  y.min <- rollapply(zoo(-1*y_c.smooth), 2*w+1, max, align="center")
  find.min <- y.min + y_c.smooth[-c(1:w, n+1-1:w)]
  i.min <- which(find.min <= 0) + w

  list(curve.c=y_c.smooth, i.min=i.min, y.min= y[i.min], i.max=i.max, y.max= y[i.max]) # y values not in use

}


#------------------------------

get_peaks_range <- function(ch, target_DT, PF=0.75) {
  local.min <- rbind (target_DT[chrom==ch][valley.bottom == TRUE],
                      target_DT[chrom==ch][global_indx == min(global_indx) | global_indx == max(global_indx),])
  keycols <- c("chrom", "start", "end")
  setkeyv(local.min, cols=keycols , physical = TRUE) # haskey(BP_DT)
  peaks_limit <- as.data.table( rollapply(local.min$global_indx, width=2, FUN=function(x) {  x  }))

  #PF <- seq(0, 1, by = 0.1) #I should loop over this
  #peaks_range <- lapply(PF, get_peaks_range_by_cutoff, target_DT, peaks_limit) # here I have the peaks ranges 0.75 --> 25
  #peaks_range <- rbindlist(peaks_range, use.names=FALSE, fill=FALSE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value

  peaks_range =  peaks_limit[, get_range2(V1, V2, target_DT[chrom==ch], PF), by = seq_len(nrow(peaks_limit))] #1:NROW(peaks_limit)
  #  work_BP_table_temp[as.data.table( peaks_range), on=.(chrom, global_indx), peak_range := TRUE ] # annotating the external table with the BP within peaks
  return(peaks_range)
}


get_range2 <- function(L, R, target_DT, PF) {
 # browser()

  peak_range <- target_DT[ global_indx > L & global_indx < R]
  if (peak_range[,sum( !is.na(peak.top))] > 0) {
    baseline <- if (target_DT[global_indx == L]$curve.c < target_DT[global_indx == R]$curve.c) target_DT[global_indx == L]$curve.c else target_DT[global_indx == R]$curve.c
    top <- (max(peak_range$curve.c)-baseline)*PF
    peak_range <- peak_range[curve.c > baseline+top]
    #some Peaks top are not in the sumit of the regions brain chr16_2
    if (peak_range[,sum( !is.na(peak.top))] == 0) {peak_range <- NULL}
  }else{peak_range <- NULL}
  return(peak_range)
}


get_peaks_range_by_cutoff <- function(PF, target_DT){
  list_peaks_ranges <- lapply(chr, get_peaks_range, target_DT, PF) # here I have the peaks ranges 0.75 --> 25
  peaks_ranges_table <- rbindlist(list_peaks_ranges, use.names=TRUE, fill=TRUE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value
  peaks_ranges_table$height_cutoff <- PF
  return(peaks_ranges_table)
}

get_peaks_info <- function(n, target_DT){
#  browser()
  peaks_ranges_table <- rbindlist(target_DT[n], use.names=FALSE, fill=FALSE, idcol=FALSE)
  # peaks_ranges_table <- rbindlist(list_peaks_ranges_table[1], use.names=FALSE, fill=FALSE, idcol=FALSE)

  peaks_ranges_table[, peak.ID := do.call(paste, .(chrom, seq_len, sep = "_"))]
  peaks_ranges_table[, N_sample := uniqueN(donor_id) , by = peak.ID]
  peaks_ranges_table[, N_SV := uniqueN(sv_id) ,by = peak.ID]

  peaks_ranges_table[, N_sample_cancer_type := uniqueN(donor_id) ,by = .(peak.ID, cancer_type) ]
  peaks_ranges_table[, Fract_sample_cancer_type := N_sample_cancer_type/N_sample ]

  peaks_ranges_table[, N_BP := .N ,by = peak.ID]
  peaks_ranges_table[ , N_BP_sample := .N , by = .(peak.ID, donor_id)]
  peaks_ranges_table[ , outliers_N_BP_sample := sum(boxplot.stats(N_BP_sample)$out > quantile(N_BP_sample, .75)), by = .(peak.ID)]

  peaks_ranges_table[, c("range_s", "range_d") := .(min(start), max(end)), by = peak.ID]
  peaks_ranges_table[, range_length := range_d - range_s]
  peaks_ranges_table[, Nsample_N_BP.fract := N_sample/N_BP]
  peaks_ranges_table[, peak.area := sum(curve.c - min(curve.c)) ,by = peak.ID]
  peaks_ranges_table[, peak.heights := max(curve.c - min(curve.c)) ,by = peak.ID] #
  peaks_ranges_table[, Recu_Posit_Sele := sqrt((N_sample/N_SV)*(peak.area/range_length))]

  peaks_ranges_table <- peaks_ranges_table[Sample_count.info, on = .(cancer_type)] # join the infomation of the number of sample per cohort
  peaks_ranges_table <- peaks_ranges_table[!is.na(peak.ID),]
  if(dim(Sample_count.info)[1] > 1) {    #binom.test
    peaks_ranges_table[ , cancer_type_enrichment.pvalue := prop.test( x = N_sample_cancer_type, n = N_sample, p= cancer_type.sample_fract, alternative = c("greater"), conf.level = 0.95)$p.value, by = .(peak.ID)]
    peaks_ranges_table[ , cancer_type_enrichment.qvalue.BH := p.adjust(cancer_type_enrichment.pvalue, method = "BH")]
    uniqueN(peaks_ranges_table[,peak.ID])

    cancer_type_enrichment.peak_data <- unique(peaks_ranges_table[peaks_ranges_table[, .I[Fract_sample_cancer_type == max(Fract_sample_cancer_type)], by=peak.ID]$V1] [,.(peak.ID, cancer_type, N_sample, N_sample_cancer_type, Fract_sample_cancer_type, cancer_type.sample_fract, cancer_type_enrichment.pvalue, cancer_type_enrichment.qvalue.BH)])
    cancer_type_enrichment.peak_data <- unique(cancer_type_enrichment.peak_data[cancer_type_enrichment.peak_data[, .I[cancer_type_enrichment.pvalue == min(cancer_type_enrichment.pvalue)], by=peak.ID]$V1] [,.(peak.ID, cancer_type, N_sample, N_sample_cancer_type, Fract_sample_cancer_type, cancer_type.sample_fract, cancer_type_enrichment.pvalue, cancer_type_enrichment.qvalue.BH)])
    cancer_type_enrichment.peak_data <- unique(cancer_type_enrichment.peak_data[cancer_type_enrichment.peak_data[, .I[cancer_type.sample_fract == max(cancer_type.sample_fract)], by=peak.ID]$V1] [,.(peak.ID, cancer_type, N_sample, N_sample_cancer_type, Fract_sample_cancer_type, cancer_type.sample_fract, cancer_type_enrichment.pvalue, cancer_type_enrichment.qvalue.BH)])

  } else {
    peaks_ranges_table$cancer_type_enrichment.pvalue <- NA
    peaks_ranges_table$cancer_type_enrichment.qvalue.BH <- NA
    uniqueN(peaks_ranges_table[,peak.ID])

    cancer_type_enrichment.peak_data <- unique(peaks_ranges_table[,.(peak.ID, cancer_type, N_sample, N_sample_cancer_type, Fract_sample_cancer_type, cancer_type.sample_fract, cancer_type_enrichment.pvalue, cancer_type_enrichment.qvalue.BH)])
  }


  peak_data <- unique(peaks_ranges_table[!is.na(peak.top)][, .(chrom, range_s, range_d, peak.ID, global_indx, curve.c,
                                                               range_length, N_BP, Nsample_N_BP.fract, outliers_N_BP_sample, N_SV,
                                                               peak.area, peak.heights, Recu_Posit_Sele,
                                                               height_cutoff), by = .(peak.ID) ]) #, peak.pvalue
  peak_data <- unique(peak_data[peak_data[, .I[curve.c == max(curve.c)], by=peak.ID]$V1] ) # fix if there is two summit in the peak [check out the get_peaks_range for possible update]

  peak_data <- peak_data[cancer_type_enrichment.peak_data, on = .(peak.ID)] # join the infomation of the number of sample per cohort and significant enreahminet

  if(dim(peak_data)[1] > 1){
  fit_g <- fitdist(data.frame(Value = peak_data$peak.heights)[,c("Value")], "gamma", method = c("mge")) # method = c("mle", "mme", "qme", "mge")
  #plotdist(data.frame(Value =  abs(peak_data$peak.heights))[,c("Value")], histo = TRUE, demp = TRUE, distr = "gamma", para = list(shape=fit_g$estimate[1], rate = fit_g$estimate[2]))
  peak_data[, peak.pvalue := .(Value = 1- pgamma(abs(peak_data$peak.heights), shape=fit_g$estimate[1], rate = fit_g$estimate[2], log = FALSE) )]
  peak_data[ ,peak.qvalue.BH := p.adjust(peak.pvalue, method = "BH")]

  single_sample.gamma_dist <- fitdist(data.frame(Value = (peak_data$N_BP/peak_data$N_sample))[,c("Value")], "gamma", method = c("mge"), gof = "ADL") # method = c("mle", "mme", "qme", "mge")
  #plotdist(data.frame(Value =  abs( peak_data$N_BP/peak_data$N_sample))[,c("Value")], histo = TRUE, demp = TRUE, distr = "gamma", para = list(shape=single_sample.gamma_dist$estimate[1], rate = single_sample.gamma_dist$estimate[2]))
  peak_data[, chrom.pvalue := .(Value = 1- pgamma(abs(peak_data$N_BP/peak_data$N_sample), shape=single_sample.gamma_dist$estimate[1], rate = single_sample.gamma_dist$estimate[2], log = FALSE) )]
  peak_data[ ,chrom.qvalue.BH := p.adjust(chrom.pvalue, method = "BH")]

  recurrent.gamma_dist <- fitdist(data.frame(Value = peak_data$Recu_Posit_Sele )[,c("Value")], "gamma", method = c("mge"), gof = "ADL") # method = c("mle", "mme", "qme", "mge")
  #plotdist(data.frame(Value =  peak_data$Recu_Posit_Sele )[,c("Value")], histo = TRUE, demp = TRUE, distr = "gamma", para = list(shape = recurrent.gamma_dist$estimate[1], rate = recurrent.gamma_dist$estimate[2]))
  peak_data[, sample_recurrence.pvalue := .(Value = 1- pgamma(abs(peak_data$Recu_Posit_Sele), shape=recurrent.gamma_dist$estimate[1], rate = recurrent.gamma_dist$estimate[2], log = FALSE) )]
  peak_data[ ,sample_recurrence.qvalue.BH := p.adjust(sample_recurrence.pvalue, method = "BH")]

  return( list(peak_data, recurrent.gamma_dist, single_sample.gamma_dist, fit_g) )
  }else {return( list(peak_data, NA, NA, NA) )}
}


get_genes_in_all_peak <- function( peak_data_bed, protein_coding_genes){

  genes_in_peak <- foverlaps(peak_data_bed, protein_coding_genes, type="any")
  genes_in_peak <- genes_in_peak[!is.na(GeneName)]

  ##--write.table(genes_in_peak, paste0("Result_",data_set,"_genes_in_peak.txt"))
  #---
  return(genes_in_peak)

}

get_significant_peak_info <- function(peaks_ranges_table, peak_data_info, threshold){

  peaks_ranges_table[ , peak.ID := do.call(paste, .(chrom, seq_len, sep = "_"))]
  significant_peak_info <- peaks_ranges_table[peak_data_info[N_sample > 2 & sample_recurrence.qvalue.BH < threshold, .(peak.ID, range_s, range_d)], on=.( peak.ID), ]  #peaks_ranges_table_MET[peak.ID %in% c("chr17_1"), ]
  colnames(significant_peak_info)[2] <- "chr"
  significant_peak_info <- significant_peak_info[,-c(1)]
  colnames(significant_peak_info)[2:3] <- c("BP.start", "BP.end")#---
  setkey(significant_peak_info, chr, BP.start, BP.end)

  return(significant_peak_info)
}


get_element_in_peak <- function(significant_peak_info, peak_data_info, genome_element_bed, element){
    # add the selected genes within broad peak

  element_in_peak <- foverlaps(significant_peak_info, genome_element_bed, type="any")
  if(element != "enhancer"){
    element_in_peak <- unique(element_in_peak[!is.na(GeneName), .(start, end, GeneName, peak.ID)])
  }else{ element_in_peak <- element_in_peak[!is.na(GeneName), .(start, end, GeneName, peak.ID, Enhancer_ID)] # early version has , i.start, i.end
}
  element_in_peak <- peak_data_info[element_in_peak, on=.( peak.ID), ]

  ##--write.table(genes_in_peak, paste0("Result_",data_set,"_",element,"_in_peak.txt"))
  #---

  return(element_in_peak)

}

get_enhancer_bed <- function(enhancer_data, enhancer_mark){

  enhancer_bed <- foverlaps(enhancer_mark, enhancer_data, type="any", mult = "last", nomatch=NA)
  enhancer_bed <- enhancer_bed[!is.na(GeneName)]

  return(enhancer_bed)

}

get_significant_best_candidate_in_BP_peak <- function(significant_peak_info, peak_data_info, genome_element_bed, element) {
  # add the selected genes within broad peak

  significant_in_BP_peak <- significant_peak_info
  setkey(significant_in_BP_peak, chr, BP.start, BP.end)
  significant_in_BP_peak <- foverlaps(significant_in_BP_peak, genome_element_bed, type="any")
  significant_in_BP_peak <- significant_in_BP_peak[!is.na(GeneName)]

  significant_in_BP_peak[, N_sample := uniqueN(donor_id), by = peak.ID]
  significant_in_BP_peak[, N_sample_GENE := uniqueN(donor_id) , by = .(GeneName, peak.ID)]
  significant_in_BP_peak[, N_SV := uniqueN(sv_id), by = peak.ID]
  significant_in_BP_peak[, N_SV_GENE := uniqueN(sv_id), by = .(GeneName, peak.ID)]
  significant_in_BP_peak[, Rear_score := (N_sample_GENE/N_sample)*(N_SV_GENE/(end-start)), by = .(GeneName, peak.ID)]

  significant_best_candidate_in_BP_peak <- significant_in_BP_peak[significant_in_BP_peak[, .I[Rear_score == max(Rear_score)], by=peak.ID]$V1]

  significant_best_candidate_in_BP_peak <- unique(significant_best_candidate_in_BP_peak[N_sample_GENE > 1 ,.(start, end, GeneName, element_type, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])
  significant_best_candidate_in_BP_peak <- peak_data_info[significant_best_candidate_in_BP_peak, on=.( peak.ID), ]

  ##--write.table(significant_insulators_best_candidate_in_BP_peak, paste0("Result_",data_set,"_significant_",element,"_best_candidate_in_BP_peak.txt"))

  return(significant_best_candidate_in_BP_peak)
  }


get_significant_best_candidate_in_Rearrange_peak <- function(significant_peak_info, peak_data_info, genome_element_bed, element){
#  browser()
significant_in_Rearrange_peak <- significant_peak_info
significant_in_Rearrange_peak <- significant_in_Rearrange_peak[, BP_pair := uniqueN(BP.start), by = .(sv_id, peak.ID)]

significant_in_Rearrange_peak[ (BP_pair == 1) & (BP_order == 1) & (sv_class != "TRA"), c("Rearrange.start", "Rearrange.end") := .(BP.start, range_d) ]
significant_in_Rearrange_peak[ (BP_pair == 1) & (BP_order == 2) & (sv_class != "TRA"), c("Rearrange.start", "Rearrange.end") := .(range_s, BP.end) ]
significant_in_Rearrange_peak[ sv_class == "TRA" , c("Rearrange.start", "Rearrange.end") := .(as.integer(BP.start-5000) ,as.integer(BP.end+5000) ) ]

#Both_BP_temp <- significant_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 1), c(1:46)][ significant_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 2), c("sv_id", "BP.end")], on=.( sv_id)  ]
Both_BP_temp <- significant_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 1), ][ significant_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 2), c("sv_id", "BP.end")], on=.( sv_id)  ][,-c("Rearrange.end")]
#colnames(Both_BP_temp)[47] <- "Rearrange.end"
colnames(Both_BP_temp)[dim(Both_BP_temp)[2]] <- "Rearrange.end"
Both_BP_temp[ , Rearrange.start := BP.start ]

significant_in_Rearrange_peak <- rbindlist(list(significant_in_Rearrange_peak[(BP_pair == 1)], Both_BP_temp), use.names=TRUE, fill=TRUE, idcol=NULL)

setkey(significant_in_Rearrange_peak, chr, Rearrange.start, Rearrange.end)
significant_in_Rearrange_peak <- foverlaps(significant_in_Rearrange_peak, genome_element_bed, type="any") #ehancer_bed   #peak_data_ehancer[sv_class == "DUP",]
significant_in_Rearrange_peak <- significant_in_Rearrange_peak[!is.na(GeneName)]

significant_in_Rearrange_peak[, N_sample := uniqueN(donor_id), by = peak.ID]
significant_in_Rearrange_peak[, N_sample_GENE := uniqueN(donor_id) , by = .(GeneName, peak.ID)]
significant_in_Rearrange_peak[, N_SV := uniqueN(sv_id), by = peak.ID]
significant_in_Rearrange_peak[, N_SV_GENE := uniqueN(sv_id), by = .(GeneName, peak.ID)]
significant_in_Rearrange_peak[, N_BP_GENE := length(duplicated(sv_id)), by = .(GeneName, peak.ID)]
significant_in_Rearrange_peak[, Rear_score := (N_sample_GENE/N_sample)*(N_SV_GENE/(end-start)), by = .(GeneName, peak.ID)]

# if there is 0 alt least take out the warning
significant_best_candidate_in_Rearrange_peak <- significant_in_Rearrange_peak[significant_in_Rearrange_peak[, .I[Rear_score == max(Rear_score)], by=peak.ID]$V1]

significant_best_candidate_in_Rearrange_peak <- unique(significant_best_candidate_in_Rearrange_peak[N_sample_GENE > 1 ,.(start, end, GeneName, element_type, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])
significant_best_candidate_in_Rearrange_peak <- peak_data_info[significant_best_candidate_in_Rearrange_peak, on=.( peak.ID), ]

#write.table(significant_best_candidate_in_Rearrange_peak, paste0("Result_",data_set,"_significant_",element,"_best_candidate_in_Rearrange_peak.txt"))


return(significant_best_candidate_in_Rearrange_peak)

}


result_plot <- function( cohort = NA,
                     input,              # data frame input from scanone
                     input_peak = NA,              # data frame input from scanone
                     candidate = NA,  # optional data to candidate labels
                     chrs = NA,          # chromosomes to display
                     rug = FALSE,        # plot marker positions as rug?
                     ncol = NA,          # number of columns for facetting
                     threshold = 0.25,
                     legend = "none"
) {
  #browser()
  # if not all chromosomes should be displayed, subset input
  if (!is.na(chrs)[1]) {
    input <- input[as.character(input$chrom) %in% chrs, ]
    if (!is.na(input_peak)[1])
      input_peak <- input_peak[as.character(input_peak$chrom) %in% chrs, ]
    if (!is.na(candidate)[1])
      candidate <- candidate[as.character(candidate$chrom) %in% chrs, ]
  }
  # if no number of columns for facetting is defined, plot all in one row
  if (is.na(ncol)) {
    ncol <- length(unique(input$chrom))
  }
   # plot input data frame position and LOD score
  plot <- ggplot(input, aes(x = global_indx, y = y_c)) + {
  } +  {
    # plot rug on bottom, if TRUE
    if (rug) geom_rug( aes(colour= input$cancer_type), size = 0.05, sides = "b", length = unit(0.05, "npc"))
  } + {
    ####### if input has column method but not group, plot line and color by method
    if (!is.null(input$y_c) & is.null(input$group))
      geom_smooth( aes(global_indx, y_c, group=chrom), color= "grey50", method = "loess", formula = y ~ x, se = TRUE, span=0.25, size=0.3, show.legend= FALSE)
  } + {
    ####### if input has column method but not group, plot line and color by method
#    if (!is.na(input_peak) & !is.null(input$y_c) & is.null(input$group))
#      geom_segment(data= input_peak ,#[N_sample > 1 & sample_recurrence.qvalue.BH < 0.1, .(chrom, global_indx, curve.c, cancer_type, cancer_type_enrichment.qvalue.BH)],
#                   aes( x = global_indx, y = curve.c, xend = global_indx, yend = -Inf, color= cancer_type), size=0.2)
  } + {
    scale_colour_brewer(name = "Cancer type", palette = "Dark2" ) #"Paired" Accent
  } + {
    new_scale("colour")
  } + {
    if (!is.na(input_peak) & !is.null(input$y_c) & is.null(input$group))
      geom_point( data= input_peak ,#[, .(chrom, global_indx, curve.c, sample_recurrence.qvalue.BH)],
                  aes(global_indx, curve.c, colour = sample_recurrence.qvalue.BH),
                  size = ifelse((input_peak$N_sample < 3 & (input_peak$peak.pvalue < threshold & input_peak$chrom.qvalue.BH < threshold)) | ((input_peak$outliers_N_BP_sample > 0 & input_peak$outliers_N_BP_sample < 3) & (input_peak$peak.pvalue < threshold & input_peak$chrom.qvalue.BH < threshold)), 1.5, 1),
                  shape = ifelse((input_peak$N_sample < 3 & (input_peak$peak.pvalue < threshold & input_peak$chrom.qvalue.BH < threshold)) | ((input_peak$outliers_N_BP_sample > 0 & input_peak$outliers_N_BP_sample < 3) & (input_peak$peak.pvalue < threshold & input_peak$chrom.qvalue.BH < threshold)), 8, 19)
      ) ##cut(sample_recurrence.qvalue.BH, seq(0,1,0.1), right=FALSE)#, labels=c(0:10)
  } + {
    scale_colour_gradient2("Recurrence\np-value", midpoint = 0.1, low = "brown4", mid = "brown", high = "bisque")
  }+{
    if(!is.na(input_peak) )
      if (dim(input_peak[ (N_sample < 3 & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)) | ((outliers_N_BP_sample > 0 & outliers_N_BP_sample < 3) & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)), .(chrom, global_indx, curve.c, peak.qvalue.BH ,chrom.qvalue.BH)])[1] != 0)
        geom_point( data= input_peak [(N_sample < 3 & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)) | ((outliers_N_BP_sample > 0 & outliers_N_BP_sample < 3) & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)), .(chrom, global_indx, curve.c, peak.qvalue.BH ,chrom.qvalue.BH)],
                    aes(global_indx, curve.c), color = "black", shape = 8, size=1.5)
  } + {
    new_scale("colour")
  } + {
    if( (!is.na(candidate)[1])) # & (dim(candidate[N_sample_GENE > 1 & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, GeneName, file)])[1] != 0))
      geom_label_repel( data=  merge(candidate, input_peak, by = c("peak.ID", "chrom"), all.x = TRUE)[N_sample_GENE > 1  & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, GeneName, file)],
                        aes(global_indx, curve.c, label = GeneName, colour = file),
                        seed = 1,
                        size = 3,
                        force = 2,
                        point.padding = 0.25,
                        box.padding = 0.25,
                        nudge_y = 0.1,
      )
  } + {
    scale_colour_hue("Driver", h = c(0, 360), l=20)
    scale_colour_manual(values = c( 'Enhancer' = "red", 'Gene' = "brown", "Insulator" = "blue", 'LncRNA' = "orange", "purple"))
  } +
    # facet by chromosome
    facet_wrap(~ chrom, ncol = ncol, scales = "free_x") +
    # minimal plotting theme
    theme_void() + #_minimal()
    # increase strip title size
    theme(strip.text = element_text(face = "bold", size = 8.5),
          panel.grid.minor = element_blank(),
          #axis.text.x = element_text(color="black",  size=8, angle=90, hjust = 1, vjust = 0.5), # face="bold",
          axis.text.x=element_blank(),# Remove x axis tick labels
          #axis.title.x=element_blank(),
          axis.text.y = element_text(color="black",  size=9, hjust = 1, vjust = 1), # face="bold",
          legend.key.size = unit(0.4, "cm"),
          legend.key.width = unit(0.4,"cm"),
          legend.position = legend,
          axis.title.x = element_text(color="black",  size=9, hjust = 0.5, vjust = 0.5), # face="bold",
          axis.title.y = element_text(color="black",  size=9, angle=90, hjust = 0.5, vjust = 0.5) # face="bold",
          #plot.margin = unit(c(3,3,3,3), "lines")
    ) +
    labs(title = cohort, x = "breakpoint order-index", y = "BPpc", color = "", linetype = "")
  # use RcolorBrewer palette
  # scale_color_brewer(palette = "Set1") +
  # Change plot labels
  #print(plot)
return(plot)
}







###################

# my.plot.gam <- function (x, residuals = FALSE, rug = NULL, se = TRUE, pages = 0,
#                 select = NULL, scale = -1, n = 100, n2 = 40, n3 = 3, pers = FALSE,
#                 theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL,
#                 main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE,
#                 shade = FALSE, shade.col = "gray80", shift = 0, trans = I,
#                 seWithMean = FALSE, unconditional = FALSE, by.resids = FALSE,
#                 scheme = 0, ...)
# {
#   sub.edf <- function(lab, edf) {
#     pos <- regexpr(":", lab)[1]
#     if (pos < 0) {
#       pos <- nchar(lab) - 1
#       lab <- paste(substr(lab, start = 1, stop = pos), ")", sep = "")
#             #paste(substr(lab, start = 1, stop = pos), ",", round(edf, digits = 2), ")", sep = "")
#     }
#     else {
#       lab1 <- substr(lab, start = 1, stop = pos - 2)
#       lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
#       lab <- paste( lab1, ")", sep = "")
#              #paste( lab1,  ",", round(edf, digits = 2), lab2, sep = "")
#     }
#     lab
#   }
# ##
#   sub.edf.main <- function(lab, edf) {
#     pos <- regexpr(":", lab)[1]
#     if (pos < 0) {
#       pos <- nchar(lab) - 1
#       lab <- paste(substr(lab, start = 1, stop = pos), ",", round(edf, digits = 2), ")", sep = "")
#     }
#     else {
#       lab1 <- substr(lab, start = 1, stop = pos - 2)
#       lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
#       lab <- paste( lab1,  ",", round(edf, digits = 2), lab2, sep = "")
#     }
#
#     pos <- regexpr("_", lab)[1]
#     if (substr(lab, start =  nchar(lab), stop = nchar(lab)) == "0") {
#       lab <- paste( "Out of", substr(lab, start =  pos+1, stop = nchar(lab)-2), sep = " ")
#     }
#     else {
#       lab <- paste( "Within", substr(lab, start =  pos+1, stop = nchar(lab)-2), sep = " ")
#     }
#
#    lab
#   }
# ##
#   if (pers)
#     warning("argument pers is deprecated, please use scheme instead")
#   if (is.null(rug))
#     rug <- if (nrow(x$model) > 10000)
#       FALSE
#   else TRUE
#   if (unconditional) {
#     if (is.null(x$Vc))
#       warning("Smoothness uncertainty corrected covariance not available")
#     else x$Vp <- x$Vc
#   }
#   w.resid <- NULL
#   if (length(residuals) > 1) {
#     if (length(residuals) == length(x$residuals))
#       w.resid <- residuals
#     else warning("residuals argument to plot.gam is wrong length: ignored")
#     partial.resids <- TRUE
#   }
#   else partial.resids <- residuals
#   m <- length(x$smooth)
#   if (length(scheme) == 1)
#     scheme <- rep(scheme, m)
#   if (length(scheme) != m) {
#     warn <- paste("scheme should be a single number, or a vector with",
#                   m, "elements")
#     warning(warn)
#     scheme <- rep(scheme[1], m)
#   }
#   order <- if (is.list(x$pterms))
#     unlist(lapply(x$pterms, attr, "order"))
#   else attr(x$pterms, "order")
#   if (all.terms)
#     n.para <- sum(order == 1)
#   else n.para <- 0
#   if (se) {
#     if (is.numeric(se))
#       se2.mult <- se1.mult <- se
#     else {
#       se1.mult <- 2
#       se2.mult <- 1
#     }
#     if (se1.mult < 0)
#       se1.mult <- 0
#     if (se2.mult < 0)
#       se2.mult <- 0
#   }
#   else se1.mult <- se2.mult <- 1
#   if (se && x$Vp[1, 1] < 0) {
#     se <- FALSE
#     warning("No variance estimates available")
#   }
#   if (partial.resids) {
#     if (is.null(w.resid)) {
#       if (is.null(x$residuals) || is.null(x$weights))
#         partial.resids <- FALSE
#       else {
#         wr <- sqrt(x$weights)
#         w.resid <- x$residuals * wr
#       }
#     }
#     if (partial.resids)
#       fv.terms <- predict(x, type = "terms")
#   }
#   pd <- list()
#   i <- 1
#   if (m > 0)
#     for (i in 1:m) {
#       first <- x$smooth[[i]]$first.para
#       last <- x$smooth[[i]]$last.para
#       edf <- sum(x$edf[first:last])
#       term.lab <- sub.edf(x$smooth[[i]]$label, edf)
#       attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
#       term.main <- sub.edf.main(x$smooth[[i]]$label, edf) ###
#       P <- plot(x$smooth[[i]], P = NULL, data = x$model,
#                 partial.resids = partial.resids, rug = rug, se = se,
#                 scale = scale, n = n, n2 = n2, n3 = n3, pers = pers,
#                 theta = theta, phi = phi, jit = jit, xlab = xlab,
#                 ylab = ylab, main = term.main, label = term.lab, ylim = ylim,
#                 xlim = xlim, too.far = too.far, shade = shade,
#                 shade.col = shade.col, se1.mult = se1.mult, se2.mult = se2.mult,
#                 shift = shift, trans = trans, by.resids = by.resids,
#                 scheme = scheme[i], ...)
#       if (is.null(P))
#         pd[[i]] <- list(plot.me = FALSE)
#       else if (is.null(P$fit)) {
#         p <- x$coefficients[first:last]
#         offset <- attr(P$X, "offset")
#         if (is.null(offset))
#           P$fit <- P$X %*% p
#         else P$fit <- P$X %*% p + offset
#         if (!is.null(P$exclude))
#           P$fit[P$exclude] <- NA
#         if (se && P$se) {
#           if (seWithMean && attr(x$smooth[[i]], "nCons") >
#               0) {
#             if (length(x$cmX) < ncol(x$Vp))
#               x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
#             if (seWithMean == 2)
#               x$cmX[-(1:x$nsdf)] <- 0
#             X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp),
#                          byrow = TRUE)
#             meanL1 <- x$smooth[[i]]$meanL1
#             if (!is.null(meanL1))
#               X1 <- X1/meanL1
#             X1[, first:last] <- P$X
#             se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) *
#                                              X1)))
#           }
#           else se.fit <- sqrt(pmax(0, rowSums((P$X %*%
#                                                  x$Vp[first:last, first:last, drop = FALSE]) *
#                                                 P$X)))
#           if (!is.null(P$exclude))
#             se.fit[P$exclude] <- NA
#         }
#         if (partial.resids) {
#           P$p.resid <- fv.terms[, length(order) + i] +
#             w.resid
#         }
#         if (se && P$se)
#           P$se <- se.fit * P$se.mult
#         P$X <- NULL
#         P$plot.me <- TRUE
#         pd[[i]] <- P
#         rm(P)
#       }
#       else {
#         if (partial.resids) {
#           P$p.resid <- fv.terms[, length(order) + i] +
#             w.resid
#         }
#         P$plot.me <- TRUE
#         pd[[i]] <- P
#         rm(P)
#       }
#     }
#   n.plots <- n.para
#   if (m > 0)
#     for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me)
#   if (n.plots == 0)
#     stop("No terms to plot - nothing for plot.gam() to do.")
#   if (pages > n.plots)
#     pages <- n.plots
#   if (pages < 0)
#     pages <- 0
#   if (pages != 0) {
#     ppp <- n.plots%/%pages
#     if (n.plots%%pages != 0) {
#       ppp <- ppp + 1
#       while (ppp * (pages - 1) >= n.plots) pages <- pages -
#           1
#     }
#     c <- r <- trunc(sqrt(ppp))
#     if (c < 1)
#       r <- c <- 1
#     if (c * r < ppp)
#       c <- c + 1
#     if (c * r < ppp)
#       r <- r + 1
#     oldpar <- par(mfrow = c(r, c))
#   }
#   else {
#     ppp <- 1
#     oldpar <- par()
#   }
#   if (scale == -1 && is.null(ylim)) {
#     k <- 0
#     if (m > 0)
#       for (i in 1:m) if (pd[[i]]$plot.me && pd[[i]]$scale) {
#         if (se && length(pd[[i]]$se) > 1) {
#           ul <- pd[[i]]$fit + pd[[i]]$se
#           ll <- pd[[i]]$fit - pd[[i]]$se
#           if (k == 0) {
#             ylim <- c(min(ll, na.rm = TRUE), max(ul,
#                                                  na.rm = TRUE))
#             k <- 1
#           }
#           else {
#             if (min(ll, na.rm = TRUE) < ylim[1])
#               ylim[1] <- min(ll, na.rm = TRUE)
#             if (max(ul, na.rm = TRUE) > ylim[2])
#               ylim[2] <- max(ul, na.rm = TRUE)
#           }
#         }
#         else {
#           if (k == 0) {
#             ylim <- range(pd[[i]]$fit, na.rm = TRUE)
#             k <- 1
#           }
#           else {
#             if (min(pd[[i]]$fit, na.rm = TRUE) < ylim[1])
#               ylim[1] <- min(pd[[i]]$fit, na.rm = TRUE)
#             if (max(pd[[i]]$fit, na.rm = TRUE) > ylim[2])
#               ylim[2] <- max(pd[[i]]$fit, na.rm = TRUE)
#           }
#         }
#         if (partial.resids) {
#           ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
#           if (ul > ylim[2])
#             ylim[2] <- ul
#           ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
#           if (ll < ylim[1])
#             ylim[1] <- ll
#         }
#       }
#     ylim <- trans(ylim + shift)
#   }
#   if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
#       pages > 1 && dev.interactive())
#     ask <- TRUE
#   else ask <- FALSE
#   if (!is.null(select)) {
#     ask <- FALSE
#   }
#   if (m > 0)
#     for (i in 1:m) if (pd[[i]]$plot.me && (is.null(select) ||
#                                            i == select)) {
#       plot(x$smooth[[i]], P = pd[[i]], partial.resids = partial.resids,
#            rug = rug, se = se, scale = scale, n = n, n2 = n2,
#            n3 = n3, pers = pers, theta = theta, phi = phi,
#            jit = jit, xlab = xlab, ylab = ylab, main = main,
#            ylim = ylim, xlim = xlim, too.far = too.far,
#            shade = shade, shade.col = shade.col, shift = shift,
#            trans = trans, by.resids = by.resids, scheme = scheme[i],
#            ...)
#       if (ask) {
#         oask <- devAskNewPage(TRUE)
#         on.exit(devAskNewPage(oask))
#         ask <- FALSE
#       }
#     }
#   if (n.para > 0) {
#     class(x) <- c("gam", "glm", "lm")
#     if (is.null(select)) {
#       attr(x, "para.only") <- TRUE
#       termplot(x, se = se, rug = rug, col.se = 1, col.term = 1, ylabs = "Partial coeficient",
#                main = "", attr(x$pterms, "term.labels"), ...)
#     }
#     else {
#       if (select > m) {
#         select <- select - m
#         term.labels <- attr(x$pterms, "term.labels")
#         term.labels <- term.labels[order == 1]
#         if (select <= length(term.labels)) {
#           termplot(x, terms = term.labels[select], se = se, ylabs = "Partial coeficient",
#                    rug = rug, col.se = 1, col.term = 1, ...)
#         }
#       }
#     }
#   }
#   if (pages > 0)
#     par(oldpar)
#   invisible(pd)
# }
#
# environment(my.plot.gam) <- asNamespace('mgcv')
# assignInNamespace("plot.gam", my.plot.gam, ns = "mgcv")
