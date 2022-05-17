## CSVDriver
Bioinformatics tool to identify SV drivers from whole-genome sequences.

> A ShyniAppTo of the the method is being implemented to improve usability for the users that are not expert in the R programming language.

### To reproduce the result of the study it can be run Reproduce_CSVDriver_Script.R (R version 3.6.2)

> The script could run on a desktop computer. Reproduce_CSVDriver_Script.R will perform every step of CSVDriver using the SVs call from a WGS cancer cohort
 
### Running the demo  (time 30 min)

> Before running the Reproduce_CSVDriver_Script.R with the input demo gennome_annotation.zip file must be decompressed.

* To run the demo just set the cohort_name = "demo" 

> In the Reproduce_CSVDriver_Script.R  cohort_name <- cohort[1]
 
> Then run the complete script to get the full set of results  

### Input (input folder)

Demo datasets of SV and tissue-specific genomic feature annotations to evaluate the method 
 
* demo.somatic.sv.bedpe 
 
* ChromHMM/demo.bed                                                                             
 
* TAD_segment_demo_ChromMark_class.txt
 
* demo_ccre_enhancers_data_hg19.bed 

### Output (output folder)

* GAM model
GAM_demo

* Peak data info
Table_Result_demo_peak_data_info

* Significant gene driver candidate within significant rearranged peak 
Table_Result_demo_significant_Genes_best_candidate_in_peak 

* Gene linked to significant enhancer driver candidate within significant rearranged peak
 Table_Result_demo_significant_Enhancer_best_candidate_in_peak

* Set of plots for model performance and the contribution of each covariate to the GAM model
 demo_Model_covariate_contribution.pdf
 demo_Model_Fitting_and_Residual

* QQ-plot of peaks rearrangement score PRe for significant peaks
 demo_PRe_QQplot

* Graphic of BPpc
 demo_BPpc

* Graphic of driver candidates from the BPpc analysis
 demo_Driver-Cadiadtes


### Running in a full SVs dataset from cancer cohort  

> While runing on a full SV dataset from cancer cohorts the script could take several hours, because the GAM modeling the algorithm is computationally heavy.

* To run CSVDriver in a cohort the SV input files must be placed in the corresponding folders likewise the demo.

>  Before running the Reproduce_CSVDriver_Script.R, some files must be decompressed for each cancer cohort. 
  * TAD_segment_[cohort]_ChromMark_class.zip
  * ccre_cohort_Archive.zip  
  * chromHMM_cohort_Archive.zip
  * tissue_ehnacer_marks 
  
* To run the a particular cohort set the cohort_name from the list
> cohort <- c("demo", "bone", "brain", "breast", "colon", "esophagus", "kidney", "liver", "lung",  "lymph-nodes", "ovary", "pancreas", "prostate", "skin", "stomach", "uterus")
 
> cohort_name <- cohort[ cohort_index]


### CSVDriver data sources

Genes:	GENCODE Release 29 (GRCh37) Comprehensive gene annotation

Enhancer:	GENCODE 3 annotations

Replication Timing (RT):	ENCODE Repli-seq data; mean of 1 Mb window, average of 8 different tissue type cell line

Fragile sites (FS):	HumCFS database

Gene density: 	Function kpPlotDensity from the R Package kpPlotDensity (compute in 1 Mb windows)

Chromatin mark (ChromMark):	Roadmap Epigenomics Mapping Consortium

Repeat Calss (RepClass):	reeatmapsker data from UCSC Genome Browser

LAD:	Akdemir, K. C. et al. 2020

GC content (GC):	Function GCcontent from the R Package biovizBase

TAD: 3D Genome Browser
