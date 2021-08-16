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
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- genes(txdb)

kp <- plotKaryotype(genome="hg19")
kp <- kpPlotDensity(kp, all.genes)

density <- cbind(as.data.frame(kp$latest.plot$computed.values$windows ), as.data.frame( kp$latest.plot$computed.values$density) )
colnames(density)[c(1,6)] <- c("chrom", "Gene.density")

chr_arm_lengths <- as.data.table(readRDS("/Users/alm2069/Documents/Work/Projects/SV/R-code/ReCa-Driver/ReCaDriver/ints/extdata/gennome_annotation/chr_arm_lengths.rds"))

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


get_result.plot <- function(work_BP_table, peak_data_info, significant_element_best_candidate_in_peak, threshold){


result.plot <- ggplot(data = work_BP_table) + #y.fit
  geom_hline(yintercept= 0, size = 0.5, linetype = "solid") +
  geom_smooth( aes(global_indx, y_c, group=chrom), color= "grey50", method = "loess", formula = y ~ x, se = TRUE, span=0.25, size=0.3, show.legend= FALSE) + #scale_fill_brewer() +
  geom_rug(aes(global_indx, colour= cancer_type ) ) +
  # { if(uniqueN(peak_data_info[N_sample > 1 & sample_recurrence.qvalue.BH < 0.1, .(cancer_type)]) > 1)
  { if(dim(peak_data_info[N_sample > 1 & sample_recurrence.qvalue.BH < threshold, .(cancer_type)])[1] != 0)
    geom_segment(data= peak_data_info [N_sample > 1 & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, cancer_type, cancer_type_enrichment.qvalue.BH)],
                 aes( x = global_indx, y = curve.c, xend = global_indx, yend = -Inf, colour= cancer_type))
    #       xintercept = global_indx, curve.c, label = " ", colour= cancer_type))
    #   geom_text_repel( data= peak_data_info [N_sample > 1 & sample_recurrence.qvalue.BH < 0.1, .(chrom, global_indx, curve.c, cancer_type, cancer_type_enrichment.qvalue.BH)],
    #                     aes(global_indx, curve.c, label = " ", colour= cancer_type), size = 2.5, show.legend = FALSE, nudge_y = -0.15 - peak_data_info [N_sample > 1 & sample_recurrence.qvalue.BH < 0.1, .(chrom, global_indx, curve.c, cancer_type, cancer_type_enrichment.qvalue.BH)]$curve.c,
    #                     nudge_x = 0.05, direction = "y", angle = 90, vjust = 0, segment.size = 0.5, force = 0 )
  } +
  scale_colour_brewer(palette = "Paired") +
  #scale_colour_viridis_d("Cancer\nType", aesthetics = "colour",  option = "C", direction =-1) +
  new_scale("colour") +
  geom_point( data= peak_data_info [, .(chrom, global_indx, curve.c, sample_recurrence.qvalue.BH)], aes(global_indx, curve.c, colour = sample_recurrence.qvalue.BH), size = 1 ) +
  scale_colour_gradient2("p-value\nSP-scores", midpoint = 0.15, low = "brown4", mid = "brown", high = "bisque") +
  new_scale("colour") +
  { if(dim(significant_element_best_candidate_in_peak[N_sample_GENE > 1 & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, GeneName, file)])[1] != 0)
    geom_text_repel( data= significant_element_best_candidate_in_peak[N_sample_GENE > 1 & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, GeneName, file)],
                     aes(global_indx, curve.c, label = GeneName, colour = file ),
                     size = 3, force = 0.5, box.padding = 0.05,# color = file,
                     nudge_y = 0.2 + significant_element_best_candidate_in_peak [N_sample > 1 & sample_recurrence.qvalue.BH < threshold, .(chrom, global_indx, curve.c, file)]$curve.c,
                     nudge_x = 0, seed = 1,
                     min.segment.length = 0, direction = "both", angle = 90, vjust = 1, hjust = 0.5, segment.size = 0.2, segment.color="black")
  } +
  scale_colour_hue("Driver", h = c(0, 360), l=20) +
  { if(dim(peak_data_info[ (N_sample < 5 & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)) | ((outliers_N_BP_sample > 0 & outliers_N_BP_sample < 5) & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)), .(chrom, global_indx, curve.c, peak.qvalue.BH ,chrom.qvalue.BH)])[1] != 0)
    geom_text_repel( data= peak_data_info [(N_sample < 5 & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)) | ((outliers_N_BP_sample > 0 & outliers_N_BP_sample < 5) & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)), .(chrom, global_indx, curve.c, peak.qvalue.BH ,chrom.qvalue.BH)],
                     aes(global_indx, curve.c, label = "ss-spot"), size = 3, show.legend = FALSE, nudge_y = -0.2 - peak_data_info [(N_sample < 5 & (peak.pvalue < threshold & chrom.qvalue.BH < threshold))| ((outliers_N_BP_sample > 0 & outliers_N_BP_sample < 5) & (peak.pvalue < threshold & chrom.qvalue.BH < threshold)), .(chrom, global_indx, curve.c, peak.qvalue.BH, chrom.qvalue.BH)]$curve.c,
                     nudge_x = 0.05, direction = "y", angle = 90, vjust = 0, segment.size = 0.25, force = 0 )
  } +
  # scale_fill_tron()+
  #  scale_fill_manual(values = cols)  +
  theme_bw(base_size = 8.5) +
  scale_x_continuous(expand = c(0,0),
                     minor_breaks = round(work_BP_table[, max(global_indx), by = chrom]$V1),
                     breaks = round(work_BP_table[, mean(global_indx), by = chrom]$V1),
                     labels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
  ) +
  coord_cartesian(ylim=c(-0.20, 0.25)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + #labels=scaleFUN
  ##     scale_y_continuous(breaks = seq(-0.25, 0.25, 0.1)) +
  theme(
    legend.position = "top", #right
    legend.box = "horizontal", #
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
    #    axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
    axis.text.y = element_text(color="black",  size=8.5, hjust = 1, vjust = 1), # face="bold",
    axis.text.x = element_text(color="black",  size=8.5, angle=90, hjust = 1, vjust = 0.5) # face="bold",
    #  axis.title.x=element_blank(),
    #     axis.text.x=element_blank(),# Remove x axis tick labels
    #     axis.ticks.x=element_blank()# Remove ticks
  ) +
  labs(x = "BP index (sorted position)", y = "Adjusted RP-signal ")

return(result.plot)


}



###################
#get_significant_genes_best_candidate_in_Rearrange_peak <- function(significant_peak_info, peak_data_info, protein_coding_genes){
#  #--#
#  significant_genes_in_Rearrange_peak <- significant_peak_info
#  significant_genes_in_Rearrange_peak <- significant_genes_in_Rearrange_peak[, BP_pair := uniqueN(BP.start), by = .(sv_id, peak.ID)]
#
#  significant_genes_in_Rearrange_peak[ (BP_pair == 1) & (BP_order == 1) & (sv_class != "TRA"), c("Rearrange.start", "Rearrange.end") := .(BP.start, range_d) ]
#  significant_genes_in_Rearrange_peak[ (BP_pair == 1) & (BP_order == 2) & (sv_class != "TRA"), c("Rearrange.start", "Rearrange.end") := .(range_s, BP.end) ]
#  significant_genes_in_Rearrange_peak[ sv_class == "TRA" , c("Rearrange.start", "Rearrange.end") := .(as.integer(BP.start-5000), as.integer(BP.end+5000)) ]
#
#  Both_BP_temp <- significant_genes_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 1), c(1:46)][ significant_genes_in_Rearrange_peak[ (BP_pair == 2) & (BP_order == 2), c("sv_id", "BP.end")], on=.( sv_id)  ]
#  colnames(Both_BP_temp)[47] <- "Rearrange.end"
#  Both_BP_temp[ , Rearrange.start := BP.start ]
#
#  significant_genes_in_Rearrange_peak <- rbindlist(list(significant_genes_in_Rearrange_peak[(BP_pair == 1)], Both_BP_temp), use.names=TRUE, fill=TRUE, idcol=NULL)
#
#  setkey(significant_genes_in_Rearrange_peak, chr, Rearrange.start, Rearrange.end)
#  significant_genes_in_Rearrange_peak <- foverlaps(significant_genes_in_Rearrange_peak, protein_coding_genes, type="any") #ehancer_bed   #peak_data_ehancer[sv_class == "DUP",]
#  significant_genes_in_Rearrange_peak <- significant_genes_in_Rearrange_peak[!is.na(GeneName)]
#
#  significant_genes_in_Rearrange_peak[, N_sample := uniqueN(donor_id), by = peak.ID]
#  significant_genes_in_Rearrange_peak[, N_sample_GENE := uniqueN(donor_id) , by = .(GeneName, peak.ID)]
#  significant_genes_in_Rearrange_peak[, N_SV := uniqueN(sv_id), by = peak.ID]
#  significant_genes_in_Rearrange_peak[, N_SV_GENE := uniqueN(sv_id), by = .(GeneName, peak.ID)]
#  significant_genes_in_Rearrange_peak[, N_BP_GENE := length(duplicated(sv_id)), by = .(GeneName, peak.ID)]
#  significant_genes_in_Rearrange_peak[, Rear_score := (N_sample_GENE/N_sample)*(N_SV_GENE/(end-start)), by = .(GeneName, peak.ID)]
#
#  significant_genes_best_candidate_in_Rearrange_peak <- significant_genes_in_Rearrange_peak[significant_genes_in_Rearrange_peak[, .I[Rear_score == max(Rear_score)], by=peak.ID]$V1]
#
#  significant_genes_best_candidate_in_Rearrange_peak <- unique(significant_genes_best_candidate_in_Rearrange_peak[N_sample_GENE > 1 ,.(start, end, GeneName, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])
#  significant_genes_best_candidate_in_Rearrange_peak <- peak_data_info[significant_genes_best_candidate_in_Rearrange_peak, on=.( peak.ID), ]
#
#  ##--write.table(significant_genes_best_candidate_in_peak, paste0("Result_",data_set,"_significant_genes_best_candidate_in_peak.txt"))
#
#
#  return(significant_genes_best_candidate_in_Rearrange_peak)
#
#
#}


my.plot.gam <- function (x, residuals = FALSE, rug = NULL, se = TRUE, pages = 0,
                select = NULL, scale = -1, n = 100, n2 = 40, n3 = 3, pers = FALSE,
                theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL,
                main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE,
                shade = FALSE, shade.col = "gray80", shift = 0, trans = I,
                seWithMean = FALSE, unconditional = FALSE, by.resids = FALSE,
                scheme = 0, ...)
{
  sub.edf <- function(lab, edf) {
    pos <- regexpr(":", lab)[1]
    if (pos < 0) {
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab, start = 1, stop = pos), ")", sep = "")
            #paste(substr(lab, start = 1, stop = pos), ",", round(edf, digits = 2), ")", sep = "")
    }
    else {
      lab1 <- substr(lab, start = 1, stop = pos - 2)
      lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
      lab <- paste( lab1, ")", sep = "")
             #paste( lab1,  ",", round(edf, digits = 2), lab2, sep = "")
    }
    lab
  }
##
  sub.edf.main <- function(lab, edf) {
    pos <- regexpr(":", lab)[1]
    if (pos < 0) {
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab, start = 1, stop = pos), ",", round(edf, digits = 2), ")", sep = "")
    }
    else {
      lab1 <- substr(lab, start = 1, stop = pos - 2)
      lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
      lab <- paste( lab1,  ",", round(edf, digits = 2), lab2, sep = "")
    }

    pos <- regexpr("_", lab)[1]
    if (substr(lab, start =  nchar(lab), stop = nchar(lab)) == "0") {
      lab <- paste( "Out of", substr(lab, start =  pos+1, stop = nchar(lab)-2), sep = " ")
    }
    else {
      lab <- paste( "Within", substr(lab, start =  pos+1, stop = nchar(lab)-2), sep = " ")
    }

   lab
  }
##
  if (pers)
    warning("argument pers is deprecated, please use scheme instead")
  if (is.null(rug))
    rug <- if (nrow(x$model) > 10000)
      FALSE
  else TRUE
  if (unconditional) {
    if (is.null(x$Vc))
      warning("Smoothness uncertainty corrected covariance not available")
    else x$Vp <- x$Vc
  }
  w.resid <- NULL
  if (length(residuals) > 1) {
    if (length(residuals) == length(x$residuals))
      w.resid <- residuals
    else warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  }
  else partial.resids <- residuals
  m <- length(x$smooth)
  if (length(scheme) == 1)
    scheme <- rep(scheme, m)
  if (length(scheme) != m) {
    warn <- paste("scheme should be a single number, or a vector with",
                  m, "elements")
    warning(warn)
    scheme <- rep(scheme[1], m)
  }
  order <- if (is.list(x$pterms))
    unlist(lapply(x$pterms, attr, "order"))
  else attr(x$pterms, "order")
  if (all.terms)
    n.para <- sum(order == 1)
  else n.para <- 0
  if (se) {
    if (is.numeric(se))
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0)
      se1.mult <- 0
    if (se2.mult < 0)
      se2.mult <- 0
  }
  else se1.mult <- se2.mult <- 1
  if (se && x$Vp[1, 1] < 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  if (partial.resids) {
    if (is.null(w.resid)) {
      if (is.null(x$residuals) || is.null(x$weights))
        partial.resids <- FALSE
      else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals * wr
      }
    }
    if (partial.resids)
      fv.terms <- predict(x, type = "terms")
  }
  pd <- list()
  i <- 1
  if (m > 0)
    for (i in 1:m) {
      first <- x$smooth[[i]]$first.para
      last <- x$smooth[[i]]$last.para
      edf <- sum(x$edf[first:last])
      term.lab <- sub.edf(x$smooth[[i]]$label, edf)
      attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
      term.main <- sub.edf.main(x$smooth[[i]]$label, edf) ###
      P <- plot(x$smooth[[i]], P = NULL, data = x$model,
                partial.resids = partial.resids, rug = rug, se = se,
                scale = scale, n = n, n2 = n2, n3 = n3, pers = pers,
                theta = theta, phi = phi, jit = jit, xlab = xlab,
                ylab = ylab, main = term.main, label = term.lab, ylim = ylim,
                xlim = xlim, too.far = too.far, shade = shade,
                shade.col = shade.col, se1.mult = se1.mult, se2.mult = se2.mult,
                shift = shift, trans = trans, by.resids = by.resids,
                scheme = scheme[i], ...)
      if (is.null(P))
        pd[[i]] <- list(plot.me = FALSE)
      else if (is.null(P$fit)) {
        p <- x$coefficients[first:last]
        offset <- attr(P$X, "offset")
        if (is.null(offset))
          P$fit <- P$X %*% p
        else P$fit <- P$X %*% p + offset
        if (!is.null(P$exclude))
          P$fit[P$exclude] <- NA
        if (se && P$se) {
          if (seWithMean && attr(x$smooth[[i]], "nCons") >
              0) {
            if (length(x$cmX) < ncol(x$Vp))
              x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
            if (seWithMean == 2)
              x$cmX[-(1:x$nsdf)] <- 0
            X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp),
                         byrow = TRUE)
            meanL1 <- x$smooth[[i]]$meanL1
            if (!is.null(meanL1))
              X1 <- X1/meanL1
            X1[, first:last] <- P$X
            se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) *
                                             X1)))
          }
          else se.fit <- sqrt(pmax(0, rowSums((P$X %*%
                                                 x$Vp[first:last, first:last, drop = FALSE]) *
                                                P$X)))
          if (!is.null(P$exclude))
            se.fit[P$exclude] <- NA
        }
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        if (se && P$se)
          P$se <- se.fit * P$se.mult
        P$X <- NULL
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
      else {
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
    }
  n.plots <- n.para
  if (m > 0)
    for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me)
  if (n.plots == 0)
    stop("No terms to plot - nothing for plot.gam() to do.")
  if (pages > n.plots)
    pages <- n.plots
  if (pages < 0)
    pages <- 0
  if (pages != 0) {
    ppp <- n.plots%/%pages
    if (n.plots%%pages != 0) {
      ppp <- ppp + 1
      while (ppp * (pages - 1) >= n.plots) pages <- pages -
          1
    }
    c <- r <- trunc(sqrt(ppp))
    if (c < 1)
      r <- c <- 1
    if (c * r < ppp)
      c <- c + 1
    if (c * r < ppp)
      r <- r + 1
    oldpar <- par(mfrow = c(r, c))
  }
  else {
    ppp <- 1
    oldpar <- par()
  }
  if (scale == -1 && is.null(ylim)) {
    k <- 0
    if (m > 0)
      for (i in 1:m) if (pd[[i]]$plot.me && pd[[i]]$scale) {
        if (se && length(pd[[i]]$se) > 1) {
          ul <- pd[[i]]$fit + pd[[i]]$se
          ll <- pd[[i]]$fit - pd[[i]]$se
          if (k == 0) {
            ylim <- c(min(ll, na.rm = TRUE), max(ul,
                                                 na.rm = TRUE))
            k <- 1
          }
          else {
            if (min(ll, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(ll, na.rm = TRUE)
            if (max(ul, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(ul, na.rm = TRUE)
          }
        }
        else {
          if (k == 0) {
            ylim <- range(pd[[i]]$fit, na.rm = TRUE)
            k <- 1
          }
          else {
            if (min(pd[[i]]$fit, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(pd[[i]]$fit, na.rm = TRUE)
            if (max(pd[[i]]$fit, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(pd[[i]]$fit, na.rm = TRUE)
          }
        }
        if (partial.resids) {
          ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
          if (ul > ylim[2])
            ylim[2] <- ul
          ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
          if (ll < ylim[1])
            ylim[1] <- ll
        }
      }
    ylim <- trans(ylim + shift)
  }
  if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
      pages > 1 && dev.interactive())
    ask <- TRUE
  else ask <- FALSE
  if (!is.null(select)) {
    ask <- FALSE
  }
  if (m > 0)
    for (i in 1:m) if (pd[[i]]$plot.me && (is.null(select) ||
                                           i == select)) {
      plot(x$smooth[[i]], P = pd[[i]], partial.resids = partial.resids,
           rug = rug, se = se, scale = scale, n = n, n2 = n2,
           n3 = n3, pers = pers, theta = theta, phi = phi,
           jit = jit, xlab = xlab, ylab = ylab, main = main,
           ylim = ylim, xlim = xlim, too.far = too.far,
           shade = shade, shade.col = shade.col, shift = shift,
           trans = trans, by.resids = by.resids, scheme = scheme[i],
           ...)
      if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
        ask <- FALSE
      }
    }
  if (n.para > 0) {
    class(x) <- c("gam", "glm", "lm")
    if (is.null(select)) {
      attr(x, "para.only") <- TRUE
      termplot(x, se = se, rug = rug, col.se = 1, col.term = 1, ylabs = "Partial coeficient",
               main = "", attr(x$pterms, "term.labels"), ...)
    }
    else {
      if (select > m) {
        select <- select - m
        term.labels <- attr(x$pterms, "term.labels")
        term.labels <- term.labels[order == 1]
        if (select <= length(term.labels)) {
          termplot(x, terms = term.labels[select], se = se, ylabs = "Partial coeficient",
                   rug = rug, col.se = 1, col.term = 1, ...)
        }
      }
    }
  }
  if (pages > 0)
    par(oldpar)
  invisible(pd)
}

environment(my.plot.gam) <- asNamespace('mgcv')
assignInNamespace("plot.gam", my.plot.gam, ns = "mgcv")
