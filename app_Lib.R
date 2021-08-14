


#source("lib.R")


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


get_GC_covarie <- function( target_DT) {
    GC <-  as.data.table( apply(target_DT[,c("chrom", "start", "end")], 1, function(row) GCcontent(Hsapiens, GRanges(as.character(row[1]), IRanges(as.numeric(row[2])-50, as.numeric(row[3])+50)))) )
    colnames(GC) <- c("GC")
    return(GC)
}

get_covaries <- function( covarie_DT, target_DT) { #cov,

     BP_covariate <- foverlaps(target_DT, covarie_DT, type="any", mult = "last", nomatch=NA)
     BP_covariate <- BP_covariate[,-c(2,3)]
     BP_covariate <- BP_covariate[, c(2) , with=FALSE] # c(1,3:dim(BP_covariate)[2], 2)
     #colnames(BP_covariate)[c(2,3)] <- c("start", "end")
     covarie_name <- names(BP_covariate)
     if( covarie_name == "LAD" ){ BP_covariate <-  BP_covariate[is.na(get(covarie_name)) , LAD:= 0.65 ] }
     if( covarie_name == "RT" ){ BP_covariate <-  BP_covariate[is.na(get(covarie_name)) , RT:= 0 ] }
     if(covarie_name == "ChromMark"){
       BP_covariate <- BP_covariate[, as.vector(unique(BP_covariate$ChromMark)):= as.list(rep(1, uniqueN(BP_covariate$ChromMark)) )]
       for (j in unique(BP_covariate$ChromMark)) set(BP_covariate, j=j, value= ifelse(BP_covariate[["ChromMark"]]==j, "1", "0") )
     }
     else { if( covarie_name == "TAD" ){
                                        BP_covariate <-  BP_covariate[is.na(get(covarie_name)) , TAD := "Boundary" ]
                                        BP_covariate$TAD <- factor(BP_covariate$TAD, levels = c('Boundary', 'Heterochromatin', 'Repressed', 'Low', 'Low-active', 'Active'))
                                        BP_covariate$TAD <- as.factor(BP_covariate$TAD)
                                        }
            else{
              #BP_covariate[, names(BP_covariate) := lapply(.SD, function(x) {x[is.na(x)] <- "0" ; "1"})]
              BP_covariate[ !is.na(get(covarie_name)), (covarie_name) := "1" ] # eval(paste("bool", covarie_name, sep ="_")
              BP_covariate[  is.na(get(covarie_name)), (covarie_name) := "0"] # eval(paste("bool", covarie_name, sep ="_")
              ### check this BP_covariate$eval(covarie_name) <- as.factor(BP_covariate$eval(covarie_name))
            }
      }
  return(BP_covariate)
}
# make predictor as factor:


# get_GAM_nodel <- function( target_DT) {
#     gam_signal <- gam(curve.o ~ as.factor(H3K27ac_H3K4me1_H3K4me3) + as.factor(DNase) + as.factor(FS) + as.factor(TAD) +
#                     s( RT, k = c(10), bs = "cr", by = as.factor(FS)) + #
#                     s( BP.ps, k = c(10), bs = "cr", by = as.factor(FS) ) + #
#                     s( gc100Base, k = c(10), bs = "cr", by = as.factor(FS)) + #
#                     s( as.factor(TAD), , bs = "cr", by = as.factor(FS)) + #
#                     ti( RT, gc100Base, k = c(10, 10), bs = c("cr", "cr"), by = as.factor(FS)) + #
#                     ti( RT, BP.ps, k = c(10, 10), bs = c("cr", "cr"), by = as.factor(DNase)) + #
#                     ti( gc100Base, BP.ps, k = c(10, 10), bs = c("cr", "cr"), by = as.factor(DNase)), #
#                     data = target_DT,
#                     family = Gamma(link = "identity"), method="REML")#, select=TRUE  inverse
#     return(gam_signal)
# }

# work_BP_table$gc100Base <-  apply(work_BP_table[,c("chrom", "start", "end")], 1, function(row) GCcontent(Hsapiens, GRanges(as.character(row[1]), IRanges(as.numeric(row[2])-50, as.numeric(row[3])+50))))
# # input a file
# #Covaries_Files <- list.files(path="/Users/alm2069/Researchs/SV_drivers/input/input/covaries", pattern="^.+bed$", full.names=TRUE)
#
# Covaries_Files <- list.files(path="/Users/alm2069/Researchs/SV_drivers/input/pcawg_data/pcawg_consensus_1.6.161116.somatic_svs/by_cancer_type/by_organ/prostate", pattern="^.+bed$", full.names=TRUE)
# Covaries_Files[5] <- "/Users/alm2069/Researchs/SV_drivers/input/pcawg_data/pcawg_consensus_1.6.161116.somatic_svs/by_cancer_type/by_organ/common_covaries/RT.bed"
# Covaries_Files[6] <- "/Users/alm2069/Researchs/SV_drivers/input/pcawg_data/pcawg_consensus_1.6.161116.somatic_svs/by_cancer_type/by_organ/common_covaries/FS.bed"
#
#
# #Covaries_Files <- list.files(path="/Users/alm2069/Researchs/SV_drivers/input/icgc-dataset-ALL_SV/covaries/Blood", pattern="^.+bed$", full.names=TRUE)
# #Covaries_Files[5] <- "/Users/alm2069/Researchs/SV_drivers/input/input/covaries/FS.bed"
# #Covaries_Files[6] <-"/Users/alm2069/Researchs/SV_drivers/input/icgc-dataset-ALL_SV/covaries/all_RT.bed"
#
# for ( i in Covaries_Files) {
#   covaries_DT <- fread( file=i, sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE )
#   covarie_name <- unlist(strsplit(i, "/|[.]",  perl = TRUE))[lengths (strsplit(i, "/|[.]",  perl = TRUE) )-1]
#   colnames(covaries_DT) <- c("chrom", "start", "end", covarie_name)
#   setkey(covaries_DT, chrom, start, end)
#   work_BP_table <- foverlaps(work_BP_table, covaries_DT, type="any", mult = "last", nomatch=NA)
#   work_BP_table <- work_BP_table[,-c(2,3)]
#   work_BP_table <- work_BP_table[, c(1,3:dim(work_BP_table)[2],2) , with=FALSE]
#   colnames(work_BP_table)[c(2,3)] <- c("start", "end")
#   # if( covarie_name == "RT" ){
#   #   work_BP_table <-  work_BP_table[is.na(covarie_name) , RT:= 0 ]
#   # }
#   # else {
#   #    work_BP_table[ !is.na(covarie_name),  eval(paste("bool", covarie_name, sep ="_")) := 1]
#   #    work_BP_table[  is.na(covarie_name),  eval(paste("bool", covarie_name, sep ="_")) := 0]
#   #  }
# }
#
# #colnames(work_BP_table)[,c(26:33)] <- c("DNase", "H3K27ac", "H3K4m1", "H3K4m3")
# work_BP_table[is.na(RT) , RT:= 0 ]
#
# work_BP_table[is.na(FS), bool_FS:= 0 ]
# work_BP_table[!is.na(FS), bool_FS:= 1 ]
#
# work_BP_table[is.na(DNase), bool_DNase:= 0 ]
# work_BP_table[!is.na(DNase), bool_DNase:= 1 ]
#
# work_BP_table[is.na(H3K27ac) , bool_H3K27ac:= 0 ]
# work_BP_table[!is.na(H3K27ac) , bool_H3K27ac:= 1 ]
#
# work_BP_table[is.na(H3K4me1) , bool_H3K4me1:= 0 ]
# work_BP_table[!is.na(H3K4me1) , bool_H3K4me1:= 1 ]
#
# work_BP_table[is.na(H3K4me3) , bool_H3K4me3:= 0 ]
# work_BP_table[!is.na(H3K4me3) , bool_H3K4me3:= 1 ]
#
#
# work_BP_table[ ,bool_H3K27ac_H3K4me1_H3K4me3 :=as.factor(paste(bool_H3K27ac, bool_H3K4me1, bool_H3K4me3 , sep = "|")) ]
#
