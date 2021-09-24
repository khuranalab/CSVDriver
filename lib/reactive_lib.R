

###############################
####### loading SV calls ######
###############################

ISV_DT <- eventReactive(input$Input_Files_list, {
  req(input$Input_Files_list)
  tryCatch({
    Input_Files_list <- rbindlist( lapply(input$Input_Files_list$datapath, fread,  sep="\t", sep2="", dec=".", quote="\"", skip=0, header=TRUE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE) )
  }, error = function(e) {
    # return a safeError if a parsing error occurs
    stop(safeError(e))
  }
  )
  Input_Files_list
})

#########################################
#######         Transform to BP DT ######
#########################################

BP_DT <- reactive( get_BP_DT(ISV_DT()) )

#########################################
####### compute BPpc signal        ######
#########################################

signal_BP_DT <- reactive( { get_signal_BP_DT( BP_DT() ) } )

#########################################
####### loading Genomic Covariates ######
#########################################

signal_BP_Covariates_DT <- reactive( {

  signal_BP_Covariate <- signal_BP_DT()

  ##ADD LAD
  signal_BP_Covariate <- get_LAD_covariate(LAD(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## ADD RT
  signal_BP_Covariate <- get_RT_covariate(Rep.Time(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## GENE DENSITY
  signal_BP_Covariate <-  get_Gene_density_covariate(Gene.Density(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## TAD SEGMENT
  signal_BP_Covariate <-  get_TADsegment_covariate(TAD.Segment(input$cohort), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## TAD RECURRENCE
  signal_BP_Covariate <-  get_TAD_recurr_covariate(TAD.Recurrence(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## TAD BOUNDARIE RECURRENCE
  signal_BP_Covariate <-  get_TAD.B_recurr_covariate(TADB.Recurrence(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## ChromHMM Marks
  signal_BP_Covariate <-  get_ChromMark_covariate(ChromHMM.Mark(input$cohort), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## Repeat Class
  signal_BP_Covariate <-  get_RepClass_covarie(Repeat.Class(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## FS
  signal_BP_Covariate <-  get_FS_covariate(Fragile.Site(), signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

  ## GC
  signal_BP_Covariate <-  get_GC_covariate(signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)

#  write.table(signal_BP_Covariate, paste0("Input_Table_",input$cohort,".txt"))
signal_BP_Covariate
})

GAM.model <- reactive({
  gam_model <- get_GAM( signal_BP_Covariates_DT() )
#  saveRDS(gam_model, file = paste0("GAM_",input$cohort) )
gam_model
})

##

signal_covariates_BP_DT <- eventReactive(input$Genomic_covariates, {
  # input$Genomic_covariates will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
  req(input$Genomic_covariates)
  Genomic_covariates = list()
  covarie_name= list()
  signal_BP_Covariate = signal_BP_DT()
  for(nr in 1:length(input$Genomic_covariates[, 1])){
    # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
    tryCatch({
      covarie_name <- sub(".bed$", "", basename(input$Genomic_covariates[[nr, 'name']])) #$name

      ## ADD RT
        if( covarie_name == "RT" ){
          RT <- fread( file= input$Genomic_covariates[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", covarie_name[[nr]]), key = c("chrom", "start", "end") ) # , colClasses=c("character", "numeric", "numeric", "factor")
          colnames(RT) <- c("chrom", "start", "end", "RT.mean", "early", "late")
          RT[,RTodd := early/late]
          setkey(RT, chrom, start, end, physical = TRUE)
          signal_BP_Covariate <- get_RT_covariate(RT, signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)
        }
      ## ADD FS
        if( covarie_name == "FS" ){
          FS <- fread( file= input$Genomic_covariates[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "FS"), colClasses=c("character", "integer", "integer", "factor"), key = c("chrom", "start", "end") )
          signal_BP_Covariate <-  get_FS_covariate(FS, signal_BP_Covariate) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)
        }


      #  browser()
      Genomic_covariates[[ nr ]] <- fread( file= input$Genomic_covariates[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", covarie_name[[nr]]), key = c("chrom", "start", "end") ) # , colClasses=c("character", "numeric", "numeric", "factor")
      Genomic_covariates[[ nr ]] <- get_covaries(Genomic_covariates[[ nr ]], signal_BP_DT()) #foverlaps(signal_BP_DT(), Genomic_covariates[[ nr ]], type="any", mult = "last", nomatch=NA)
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      stop(safeError(e))
    }
    )
  }

  # GC content
  Genomic_covariates[[ length(input$Genomic_covariates[, 1])+1 ]] <- get_GC_covarie( signal_BP_DT())

  Covaries_list <- as.data.table( do.call("cbind", Genomic_covariates) )
  Covaries_list[ , ("H3K27ac|H3K4me1|H3K4me3") := as.factor(paste( H3K27ac, H3K4me1, H3K4me3 , sep = "|")) ]

  signal_covariates_BP_DT <- cbind( signal_BP_DT(), Covaries_list )
  write.table(signal_covariates_BP_DT, "Result_bone_02_06_2021.txt")
  signal_covariates_BP_DT

})



# BP_signal_table <- reactive({ BP_signal_table <- cbind( signal_BP_DT(), Covaries_list(), GC() )
#                               write.table(BP_signal_table, "Result_brain.txt")
#                               BP_signal_table
#                             })


