#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
#library(shinyFiles)

source("app_Lib.R")


# Define UI for application that draws a histogram
ui <- fluidPage(
navbarPage("Rearranged Cancer Drivers",
    tabPanel("App",
    # Application title
   # titlePanel("Rearranged Cancer Drivers"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(

            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(type = "tabs",
                        tabPanel("Available Data",
                                 hr(),
                                 helpText("SV calls from PCAWG."),
                                 selectInput("CancerType", "Cancer type",
                                             choices=c("",
                                                       "Brain",
                                                       "Breast",
                                                       "Esophagus",
                                                       "Lung",
                                                       "Lymph-nodes",
                                                       "Ovary",
                                                       "Pancreas",
                                                       "Prostate",
                                                       "Skin",
                                                       "Uterus")
                                             ) #colnames(WorldPhones)
                        ),
                        tabPanel("New Data",
                                 tags$hr(),
                                 h3("Structural variation data set"),
                                 fileInput("Input_Files_list", "SVs call", multiple = TRUE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 hr(),
                                 h3("Covariates"),
                                 h4("Breakpoint's genomic covariates"),
                                 #hr(),
                                 fileInput("Genomic_covariates", "DNase, H3K4me1, H3K4me3, H3K27ac, Chromatin Marks, TAD boundaries, Fragile Site, Replication Time", multiple = TRUE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),

                                 #hr(),
                                 # fileInput("H3K4me1", "H3K4me1", multiple = FALSE,
                                 #           accept = c(
                                 #               "text/csv",
                                 #               "text/comma-separated-values,text/plain",
                                 #               ".csv")
                                 # ),
                                 # #hr(),
                                 # fileInput("H3K4me3", "H3K4me3", multiple = FALSE,
                                 #           accept = c(
                                 #               "text/csv",
                                 #               "text/comma-separated-values,text/plain",
                                 #               ".csv")
                                 # ),
                                 # #hr(),
                                 # fileInput("H3K27ac", "H3K27ac", multiple = FALSE,
                                 #           accept = c(
                                 #               "text/csv",
                                 #               "text/comma-separated-values,text/plain",
                                 #               ".csv")
                                 # ),
                                 #hr(),
                                 # h4("Common Covariates"),
                                 # fileInput("Common_Covariates", "Fragile Site, Replication Time", multiple = TRUE,
                                 #           accept = c(
                                 #               "text/csv",
                                 #               "text/comma-separated-values,text/plain",
                                 #               ".csv")
                                 # ),
                                 # fileInput("Replication_Time", "Replication Time", multiple = FALSE,
                                 #           accept = c(
                                 #               ".bed", "text/csv",
                                 #               "text/comma-separated-values,text/plain",
                                 #               ".csv")
                                 # ),
                                 tags$hr(),
                                 h3("Genome annotation"),
                                 fileInput("protein_coding_genes", "Protein coding genes", multiple = FALSE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 fileInput("lincRNA", "Long noncoding RNAs", multiple = FALSE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 fileInput("CTCF_cohesin_insulator", "CTCFCohesin loop's insulator", multiple = FALSE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 fileInput("Enhanser_genes", "Enhancer-Gene links", multiple = FALSE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 h4("Tisue especific Enhancer"),
                                 #hr(),
                                 fileInput("Enhancer", "Enhancer", multiple = FALSE,
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),

                        )

            ),

        actionButton('run', 'run')

        ),

        # Main panel for displaying outputs ----
        mainPanel(
            # Show a plot of the generated distribution
            #plotOutput("distPlot"),
            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(type = "tabs",
                        tabPanel("Rearrangemet signal", withSpinner(plotOutput("distPlot"))),
                        tabPanel("GAM Model",
                                 h4("Model function"),
                                 withMathJax(),
                              #   helpText('and a fact about \\(\\pi\\):
                              #              $$\\frac2\\pi = \\frac{\\sqrt2}2 \\cdot
                              #              \\frac{\\sqrt{2+\\sqrt2}}2 \\cdot
                              #              \\frac{\\sqrt{2+\\sqrt{2+\\sqrt2}}}2 \\cdots$$'),
                              #   uiOutput('ex1'),
                              #  uiOutput('ex2'),
                              uiOutput("model_function"),
                              #verbatimTextOutput("model_function"),


                                 withSpinner(plotOutput("gamPlot", height = "900px",width = "900px"))
                                 ),
                        tabPanel("Summary", withSpinner(tableOutput("summary"))), #verbatimTextOutput("summary")
                        tabPanel("Table", withSpinner(tableOutput("table")))
                        )
            )
        )
    ),
    tabPanel("About")
)

)
# Define server logic required to draw a histogram
server <- function(input, output) {

  options(shiny.maxRequestSize=100*1024^2)

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

    BP_DT <- reactive( get_BP_DT(ISV_DT()) )

    signal_BP_DT <- reactive( { get_signal_BP_DT( BP_DT() ) } )


#########################################
####### loading Genomic Covariates ######
#########################################


    signal_covariates_BP_DT <- eventReactive(input$Genomic_covariates, {
      # input$Genomic_covariates will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
      req(input$Genomic_covariates)
      # validate(
      #   need( length(input$Genomic_covariates[, 1]) != 2, label = "Two data set corresponding to FS and RT")
      # )
      # if ( length(input$Genomic_covariates[, 1] == 2) ){
      Genomic_covariates = list()
      covarie_name= list()
      for(nr in 1:length(input$Genomic_covariates[, 1])){
        # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
        tryCatch({
          #Genomic_covariates <- lapply(input$Genomic_covariates$datapath, fread, sep="\t", sep2="", dec=".", quote="\"", skip=0, header=TRUE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE, data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "Enhancer_ID"), key = c("chrom", "start", "end"))
         # Genomic_covariates[[nr]] <- get_add_covarie( covariate_path = input$Genomic_covariates[[nr, 'datapath']], signal_BP_DT)
          #covarie_name <- unlist(strsplit(input$Genomic_covariates[[nr, 'datapath']], "/|[.]",  perl = TRUE))[lengths (strsplit(input$Genomic_covariates[[nr, 'datapath']], "/|[.]",  perl = TRUE) )-1]
          covarie_name[[nr]] <- sub(".bed$", "", basename(input$Genomic_covariates[[nr, 'name']])) #$name
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
      # } else {
      #         error = function(e) {
      #           # return a safeError if a parsing error occurs
      #           stop(safeError(e))
      #         }
      #   }
     # rbindlist(Genomic_covariates, use.names=FALSE, fill=FALSE, idcol=TRUE)

     # GC content
      Genomic_covariates[[ length(input$Genomic_covariates[, 1])+1 ]] <- get_GC_covarie( signal_BP_DT())

      Covaries_list <- as.data.table( do.call("cbind", Genomic_covariates) )
      Covaries_list[ , ("H3K27ac|H3K4me1|H3K4me3") := as.factor(paste( H3K27ac, H3K4me1, H3K4me3 , sep = "|")) ]

      signal_covariates_BP_DT <- cbind( signal_BP_DT(), Covaries_list )
      write.table(signal_covariates_BP_DT, "Result_bone_02_06_2021.txt")
      signal_covariates_BP_DT

      })


   # GC <-  reactive({ get_GC_covarie( signal_BP_DT() ) })

    # BP_signal_table <- reactive({ BP_signal_table <- cbind( signal_BP_DT(), Covaries_list(), GC() )
    #                               write.table(BP_signal_table, "Result_brain.txt")
    #                               BP_signal_table
    #                             })

  #  gam_model <- reactive({ gam_model <- get_GAM_nodel( BP_signal_table() )
  #  saveRDS(gam_model, file ="gam_signal_prostate")
  #  gam_model
  #  })

################################
####### Genome annoations ######
################################

    protein_coding_genes <- eventReactive(input$protein_coding_genes, {
        # input$protein_coding_genes will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
        req(input$protein_coding_genes)
        protein_coding_genes = list()
        for(nr in 1:length(input$protein_coding_genes[, 1])){
            # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
            tryCatch({
                protein_coding_genes[[nr]] <- fread( file= input$protein_coding_genes[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE,
                                                     data.table=TRUE, select = c(1:4), col.names = c("chr", "start", "end", "GeneName"), key = c("chr", "start", "end") )
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
            )
        }
        rbindlist(protein_coding_genes, use.names=FALSE, fill=FALSE, idcol=FALSE)
    })

    enhancer_data <- eventReactive(input$Enhanser_genes, {
        # input$Enhanser_genes will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
        req(input$Enhanser_genes)
        Enhanser_genes = list()
        for(nr in 1:length(input$Enhanser_genes[, 1])){
            # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
            tryCatch({
                Enhanser_genes[[nr]] <- fread( file= input$Enhanser_genes[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE,
                                                     data.table=TRUE, select = c(1:4), col.names = c("chr", "start", "end", "GeneName"), key = c("chr", "start", "end") )
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
            )
        }
        rbindlist(Enhanser_genes, use.names=FALSE, fill=FALSE, idcol=FALSE)
    })

    enhancer_mark_files <- eventReactive(input$Enhancer, {
        # input$Enhancer will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
        req(input$Enhancer)
        Enhancer = list()
        for(nr in 1:length(input$Enhancer[, 1])){
            # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
            tryCatch({
                Enhancer[[nr]] <- fread( file= input$Enhancer[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE,
                                               data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "Enhancer_ID"), key = c("chrom", "start", "end") )
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
            )
        }
        rbindlist(Enhancer, use.names=FALSE, fill=FALSE, idcol=FALSE)
    })

    insulator_bed <- eventReactive(input$CTCF_cohesin_insulator, {
      # input$CTCF_cohesin_insulator will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
      req(input$CTCF_cohesin_insulator)
      CTCF_cohesin_insulator = list()
      for(nr in 1:length(input$CTCF_cohesin_insulator[, 1])){
        # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
        tryCatch({
          CTCF_cohesin_insulator[[nr]] <- fread( file= input$CTCF_cohesin_insulator[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE,
                                   data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "Enhancer_ID"), key = c("chrom", "start", "end") )
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
        )
      }
      rbindlist(CTCF_cohesin_insulator, use.names=FALSE, fill=FALSE, idcol=FALSE)
    })

    LncRNA_bed <- eventReactive(input$lincRNA, {
      # input$lincRNA will be NULL initially. After the user selects and uploads a file, head of that data file by default, or all rows if selected, will be shown.
      req(input$lincRNA)
      lincRNA = list()
      for(nr in 1:length(input$lincRNA[, 1])){
        # when reading tab separated files, having a comma/semicolon separator causes `read.csv` to error
        tryCatch({
          lincRNA[[nr]] <- fread( file= input$lincRNA[[nr, 'datapath']], sep="\t", sep2="", dec=".", quote="\"", skip=0, header=FALSE, stringsAsFactors=TRUE, na.strings="NA", logical01=TRUE, strip.white=TRUE, fill=FALSE, blank.lines.skip=TRUE,
                                                 data.table=TRUE, select = c(1:4), col.names = c("chrom", "start", "end", "Enhancer_ID"), key = c("chrom", "start", "end") )
        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
        )
      }
     rbindlist(lincRNA, use.names=FALSE, fill=FALSE, idcol=FALSE)
    })


################################
####### Method processing ######
################################

    #   output$distPlot <- renderPlot({
    #                                 withProgress(message = 'Creating plot', style = "notification", value = 0.1, {
    #                                     for (i in 1:100) {
    #                                         incProgress(1/100)
    #                                         Sys.sleep(0.05)
    #                                     }
    #                                 })
    #                              get_result.plot(work_BP_table, peak_data_info, significant_element_best_candidate_in_peak)
    #
    #                             })
    #
    #
    # output$gamPlot <- renderPlot({
    #                             plot(gam_plot, pages=1, scale=0, rug=FALSE, las = 2, cex.axis=1.0, cex.lab = 1.0, pch=8, cex=1.0, scheme=2, se = TRUE, shade=TRUE, shade.col='gray90', all.terms = TRUE, seWithMean=TRUE)
    #                             })


    output$table <- renderTable({
      #BP_DT()
      #  significant_element_best_candidate_in_peak[, .(chrom, range_s, range_d, peak.ID, range_length, N_SV, cancer_type, file, GeneName)]
        #xtable
      signal_covariates_BP_DT()
    })

   # output$summary <- renderTable({
    #  BP_signal_table()
    #})

    #output$ex1 <- renderText({ input$Input_Files_list$datapath})
    #output$ex2 <- renderText({ Input_Files_list})

    output$model_function <- renderUI({
        if (!input$run) return()
        withMathJax(
            helpText('$$e^{i \\pi} + 1 = 0$$')
        )
        print( as.character(length(Covaries_Files_list)) )
    })


#    output$distPlot <- renderPlot({
#        # generate bins based on input$bins from ui.R
#        x    <- faithful[, 2]
#        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
#        hist(x, breaks = bins, col = 'darkgray', border = 'white')
#    })
}

# Run the application
shinyApp(ui = ui, server = server)
