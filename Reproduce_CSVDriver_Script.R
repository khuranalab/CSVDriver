source("./lib/lib.R")

# Cancer cohort in the study
cohort <- c("bone", "brain", "breast", "colon", "esophagus", "kidney", "liver", "lung",  "lymph-nodes", "ovary", "pancreas", "prostate", "skin", "stomach", "uterus")
# Set the cohort to be analyzed
cohort_name <- cohort[1]


# list of SV dataset for a given cohort
Input_Files_list <- list.files(path=paste0("./Input/cancer_organ/",cohort_name), pattern=".unique.bedpe$", full.names=TRUE)

#reading and combinding several file
ISV_DT <- rbindlist( lapply(Input_Files_list, fread,  sep="\t", sep2="", dec=".", quote="\"", skip=0, header=TRUE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE) )

#split SV in BP
BP_DT <- get_BP_DT(ISV_DT)

# Generating the observed BPpc
BPpc_working_table <- get_signal_BP_DT(BP_DT)


##ADD LAD
BPpc_working_table <- get_LAD_covariate(LAD(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## ADD RT
BPpc_working_table <- get_RT_covariate(Rep.Time(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## GENE DENSITY
BPpc_working_table <-  get_Gene_density_covariate(Gene.Density(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## TAD SEGMENT
BPpc_working_table <-  get_TADsegment_covariate(TAD.Segment(cohort_name), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## TAD RECURRENCE
BPpc_working_table <-  get_TAD_recurr_covariate(TAD.Recurrence(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## TAD BOUNDARIE RECURRENCE
BPpc_working_table <-  get_TAD.B_recurr_covariate(TADB.Recurrence(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## ChromHMM Marks
BPpc_working_table <-  get_ChromMark_covariate(ChromHMM.Mark(cohort_name), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## Repeat Class
BPpc_working_table <-  get_RepClass_covariate(Repeat.Class(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## FS
BPpc_working_table <-  get_FS_covariate(Fragile.Site(), BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

## GC
BPpc_working_table <-  get_GC_covariate(BPpc_working_table) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

#  write.table(signal_BP_Covariate, paste0("Input_Table_",cohort_name,".txt"))


## GAM covariate modeling
#gam_model <- get_GAM( BPpc_working_table)
gam_model <- gam(curve.o ~ FS +  LAD.factor + repClass + ChromMark +
                   s( RT.mean,  bs = "cr", by = LAD.factor) +
                   s( GC,  bs = "cr", by = LAD.factor) +
                   s( BP.TAD.recurrence, bs = "cr", by = LAD.factor) +
                   s( BP.TAD_boundary.recurrence, bs = "cr", by = LAD.factor) +
                   s( Gene.density, bs = "cr", by = LAD.factor) +
                   te( Gene.density, RT.mean, bs = c("cr", "cr"), by = TAD_segment.class) +
                   te( Gene.density, BP.TAD.recurrence, bs = c("cr", "cr"), by = TAD_segment.class) +
                   te( RT.mean, BP.TAD.recurrence, bs = c("cr", "cr"), by = TAD_segment.class)
                 , data = BPpc_working_table, na.action=na.exclude,
                 family = Gamma(link = "log"), method="REML" )

#  saveRDS(gam_model, file = paste0("GAM_",cohort_name) )

# Report of the GAM modeling
summary(gam_model)
par(mfrow=c(1,2))
gam.check(gam_model, type="response") #c("deviance","pearson","response")
plot(gam_model, pages=3, scale=0, shift = coef(gam_model)[1], residuals = FALSE, rug= FALSE, las = 2, cex.axis=1.0, cex.lab = 1.0, pch=1.0, cex=1.0, scheme=2, se = TRUE, shade=TRUE, shade.col='gray90', all.terms =	TRUE, seWithMean=TRUE)#,


#adding the predicted values and computing the corrected values
BPpc_working_table[, curve.p := gam_model$fitted.values]
BPpc_working_table[, y_c := curve.o - curve.p]

#write.table(work_BP_table, paste0("Table_Result_",cohort_name,".txt"))


## Compute the peaks

Sample_count.info <- BPpc_working_table[ , uniqueN(donor_id),by = cancer_type]
Sample_count.info$Total.Sample <- BPpc_working_table[ , uniqueN(donor_id) ]
Sample_count.info[, cancer_type.sample_fract := V1/Total.Sample]
colnames(Sample_count.info)[2] <- "Cancer_type.Sample"

chr <- unique(as.character(BPpc_working_table[,chrom]))
temp_work_BP_table <- lapply(chr, get_peaks_signal, target_DT=BPpc_working_table, window = 50, span = 0.25)
work_BP_table_temp <- rbindlist(temp_work_BP_table, use.names=TRUE, fill=TRUE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value
work_BP_table_temp$chrom <- factor(work_BP_table_temp$chrom, levels = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX'))

# summit threshold
PF <- 0.75
list_peaks_ranges <- lapply( chr, get_peaks_range, target_DT=work_BP_table_temp, PF=0.75) # here I have the peaks ranges 0.75 --> 25; 0.09 --> 0.1
peaks_ranges_table <- rbindlist(list_peaks_ranges , use.names=TRUE, fill=TRUE, idcol=FALSE) # signal_BP_DT is a copy of BP_DT with the chromosome's info plus the center.peak and fit.RDscore value

#write.table(peaks_ranges_table, paste0("Table_Result_",cohort_name,"_peak_range_data.txt"))

list_peaks_ranges_table <- lapply(PF, get_peaks_range_by_cutoff, target_DT=work_BP_table_temp) # here I have the peaks ranges 0.75 --> 25
list_peak_data_info <- get_peaks_info(1, target_DT=list_peaks_ranges_table)

# to create the qqplot of the significant
peak_data_info <- list_peak_data_info[[1]]

#write.table(peak_data_info, paste0("Table_Result_",cohort_name,"_peak_data_info.txt"))

## Plot of empirical and theoretical distributions and QQ-plot for significant peaks
recurrent.gamma_dist <- list_peak_data_info[[2]]
plotdist(data.frame(Value =  peak_data_info$Recu_Posit_Sele )[,c("Value")], histo = TRUE, demp = TRUE, distr = "gamma", para = list(shape = recurrent.gamma_dist$estimate[1], rate = recurrent.gamma_dist$estimate[2]))

peak_data_bed <- peak_data_info [, c("chrom", "range_s", "range_d", "peak.ID", "sample_recurrence.qvalue.BH")]
setkey(peak_data_bed, chrom, range_s, range_d)

#####
# Getting the Driver candidates
#####

## gene annotation
protein_coding_genes <- fread( file="./Input/gennome_annotation/protein_coding_genes_gencode.v29lift37.annotation.bed",
                                 sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE,
                                 blank.lines.skip=TRUE, data.table=TRUE )
colnames(protein_coding_genes) <-c("chr", "start", "end", "GeneName")
setkey(protein_coding_genes, chr, start, end)
protein_coding_genes <- protein_coding_genes[, element_type:= "CDS"]

## Obtaining significant peaks
significant_peak_info <- get_significant_peak_info(peaks_ranges_table, peak_data_info, 0.20)


genes_in_peak <- get_element_in_peak(significant_peak_info, peak_data_info, protein_coding_genes, "gene")

significant_genes_best_candidate_in_BP_peak <- get_significant_best_candidate_in_BP_peak(significant_peak_info, peak_data_info, protein_coding_genes, "gene")

significant_genes_best_candidate_in_Rearrange_peak <- get_significant_best_candidate_in_Rearrange_peak(significant_peak_info, peak_data_info, protein_coding_genes, "gene")

significant_element_best_candidate_in_peak1<- unique(significant_genes_best_candidate_in_BP_peak[,.(GeneName, element_type, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])
significant_element_best_candidate_in_peak2 <- unique(significant_genes_best_candidate_in_Rearrange_peak[,.(GeneName, element_type, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])


significant_genes_best_candidate_in_peak <- rbindlist(list( gene_BP = significant_genes_best_candidate_in_BP_peak,
                                                            genes_R = significant_genes_best_candidate_in_Rearrange_peak),
                                                      idcol = "file")

#write.table(significant_genes_best_candidate_in_peak, paste0("Table_Result_",cohort_name,"_significant_Genes_best_candidate_in_peak.txt"))


## enhancer annotation

cCRE_enhancer_data <- as.data.table (read.delim2("./Input/gennome_annotation/tissue_ehnacer_marks/Encode_Linked_Gene.txt", header=FALSE))
colnames(cCRE_enhancer_data) <- c("Enhancer_ID", "GeneName", "cell_line", "Method")

cCRE_enhancer_mark_files <- list.files(path = "./Input/gennome_annotation/tissue_ehnacer_marks/", pattern = paste0(cohort_name,"_ccre_enhancers_data_hg19.bed$"), full.names = TRUE) # f_name = data_set
cCRE_enhancer_mark <- fread( file= cCRE_enhancer_mark_files , sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE )
colnames(cCRE_enhancer_mark) <- c("chrom", "start", "end", "Enhancer_ID")
cCRE_enhancer_mark <- cCRE_enhancer_mark[cCRE_enhancer_data[,.(Enhancer_ID, GeneName)], on = .(Enhancer_ID), nomatch = NULL]
setkey(cCRE_enhancer_mark, chrom, start, end)
cCRE_enhancer_mark <- cCRE_enhancer_mark[, element_type:= "Enhancer"]
enhancer_bed <- cCRE_enhancer_mark

enhancer_in_peak <- get_element_in_peak(significant_peak_info, peak_data_info, enhancer_bed, "enhancer")

significant_enhancer_best_candidate_in_Rearrange_peak <- get_significant_best_candidate_in_Rearrange_peak(significant_peak_info, peak_data_info, enhancer_bed, "enhancer")

significant_element_best_candidate_in_peak3 <- unique(significant_enhancer_best_candidate_in_Rearrange_peak[,.(GeneName, element_type, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])

#write.table(significant_enhancer_best_candidate_in_Rearrange_peak, paste0("Table_Result_",cohort_name,"_significant_Enhancer_best_candidate_in_peak_info.txt"))


## insulator annotation
# insulator_bed <- fread( file= "./Input/gennome_annotation/CTCF_cohesin_insulator_6_out_of_7_ChIA_PET.intersectOverlap.April_12_2018_named.bed",
#                        sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE,
#                        blank.lines.skip=TRUE, data.table=TRUE )
# colnames(insulator_bed) <-c("chr", "start", "end", "GeneName")
# setkey(insulator_bed, chr, start, end)
# insulator_bed <- insulator_bed[, element_type:= "Insulator"]

#insulator_in_peak <- get_element_in_peak(significant_peak_info, peak_data_info, insulator_bed, "insulator")

#significant_insulators_best_candidate_in_BP_peak <- get_significant_best_candidate_in_BP_peak(significant_peak_info, peak_data_info, insulator_bed, "insulator")

#significant_element_best_candidate_in_peak4<- unique(significant_insulators_best_candidate_in_BP_peak[,.(GeneName, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])

#write.table(significant_element_best_candidate_in_peak4, paste0("Table_Result_",cohort_name,"_significant_Insulator_best_candidate_in_peak.txt"))


## LncRNA annotation
# LncRNA_bed <- fread( file="/Users/alm2069/Documents/Work/Projects/SV/R-code/ReCa-Driver/ReCaDriver/ints/extdata/gennome_annotation/lincRNA_gencode.v29lift37.long_noncoding_RNAs.bed",
#                                   sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, strip.white=TRUE, fill=FALSE,
#                                   blank.lines.skip=TRUE, data.table=TRUE ) #.input
# colnames(LncRNA_bed) <-c("chr", "start", "end", "GeneName")
# setkey(LncRNA_bed, chr, start, end)
#LncRNA_bed <- LncRNA_bed[, element_type:= "LncRNA"]

#LncRNA_in_peak <- get_element_in_peak(significant_peak_info, peak_data_info, LncRNA_bed, "LncRNA_bed")

#significant_LncRNA_best_candidate_in_BP_peak <- get_significant_best_candidate_in_BP_peak(significant_peak_info, peak_data_info, LncRNA_bed, "LncRNA_bed")

#significant_element_best_candidate_in_peak5 <- unique(significant_LncRNA_best_candidate_in_BP_peak[,.(GeneName, N_SV_GENE, N_sample_GENE, Rear_score, peak.ID)])

#write.table(significant_element_best_candidate_in_peak5, paste0("Table_Result_",cohort_name,"_significant_LncRNA_best_candidate_in_peak.txt"))


##  Plotting the results

if( nrow(significant_genes_best_candidate_in_peak) != 0){
significant_gene <- rbind(significant_genes_best_candidate_in_peak[file == "gene_BP",],
                          significant_genes_best_candidate_in_peak[file == "genes_R", ][ !(GeneName %in% as.vector( unique(significant_genes_best_candidate_in_peak[file == "gene_BP",]$GeneName) )), ]
                    )
}else{significant_gene  <- data.table()}

if( nrow(significant_enhancer_best_candidate_in_Rearrange_peak) != 0){
  significant_enhancer <- significant_enhancer_best_candidate_in_Rearrange_peak[!is.na(chrom),]
  if( nrow(significant_enhancer) != 0) {
    significant_enhancer[,file := "Enhancer"]
    significant_enhancer <- significant_enhancer[ !(GeneName %in% as.vector( unique(significant_gene$GeneName) )) ]
  }
}else{significant_enhancer  <- data.table()}

significant_element = data.table()
if(  nrow(significant_gene) != 0)
  significant_element <- rbind(significant_element, significant_gene[ sample_recurrence.qvalue.BH <0.2 ,.(peak.ID, chrom, GeneName, file, N_sample_GENE)])
if( nrow(significant_enhancer) != 0)
  significant_element <- rbind(significant_element, unique(significant_enhancer[sample_recurrence.qvalue.BH <0.2 ,.(peak.ID, chrom, GeneName, file, N_sample_GENE)]) )

#  pdf(file=paste0(cohort_name,"_BPpc.pdf"))
result_plot(cohort = cohort_name, input = BPpc_working_table, input_peak = peak_data_info, threshold = 0.20, rug = FALSE, legend="right" )
#  dev.off()

#  pdf(file=paste0(cohort_name,"_Driver-Cadiadtes.pdf"))
result_plot(cohort = cohort_name, input = BPpc_working_table, input_peak = peak_data_info, candidate = significant_element, pcawg_regions = NA, threshold = 0.20, rug = FALSE, chrs = as.vector(unique(significant_element$chrom)), ncol = 4 ) #
#  dev.off()

