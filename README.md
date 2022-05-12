# CSVDriver
Bioinformatics tool to identify SV drivers from whole-genome sequences.
The code is fully functional althought currently the algorithm is being implemented in an R shyniApp to improve usability of the finaly release for the users that are not expert in the R programming language.
To reproduce the result of the study it can be run Reproduce_CSVDriver_Script.R

# CSVDriver data sources

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
