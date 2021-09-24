
output$ui_import_data_sidebar_left <- renderUI({
#  output$ui_sel_resbtns_import <- renderUI({
    cc <- c("Testdata/internal data", "R-dataset (.rdata)", "SPSS-file (.sav)", "SAS-file (.sasb7dat)",
            "CSV-file (.csv, .txt)", "STATA-file (.dta)")
#    out <- fluidRow(column(12, h4("Select data source")))
#    for (i in 1:length(cc)) {
#      id <- paste0("btn_import_data_", i)
#      if (obj$cur_selection_import==id) {
        style <- "primary"
#      } else {
#        style <- "default"
#      }

       # btn <- actionButton("run",label=("Run CSVDriver"), "primary")
        btn_run <-  actionButton('RunCSVDriver', 'Run CSVDriver')
        observeEvent(input$btn_run, {
          #appendTab(session, inputId = "tabs",
          #          tab = tabPanel("Dynamic", "This a dynamically-added tab"),
          #          #target = "Table",
          #          select = TRUE,
          #)
          updateTabsetPanel(session, "tabs", selected = paste0("xxx"))

        })

        out <- list( fluidRow( #out,
        # column(12, bsButton(id, label=cc[i], block=TRUE, size="extra-small", style=style), tags$br())
       # column(12, bsButton("idxxx", label=cc[1], block=TRUE, size="extra-small", style=style)),
#          column(12,
#          sidebarPanel(
            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(id = "Input_tab", type = "tabs",
                        tabPanel("Input Data",
                                  #column(12, tags$hr()),
                                  #column(12, helpText("Organ of the cancer cohort")),
                                  column(12, selectInput("cohort", helpText("Cancer cohort name"),
                                                        choices=c("",
                                                                  "brain",
                                                                  "breast",
                                                                  "esophagus",
                                                                  "lung",
                                                                  "lymph-nodes",
                                                                  "ovary",
                                                                  "pancreas",
                                                                  "prostate",
                                                                  "skin",
                                                                  "uterus")
                                              )
                                  ),
                                 #column(12, h4("Structural variation")),
                                  column(12, fileInput("Input_Files_list", helpText("SVs call"), multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, tags$hr()),
                                  #column(12, h4("Covariates")),
                                  column(12, fileInput("Genomic_covariates", "Genomic covariates", multiple = TRUE,
                                                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("LAD", "LAD", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("ChromMarks", "Chromatin Marks", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("RepeatClass", "Repeat class", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("Replication_Time", "Replication Time", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, tags$hr()),
                                  column(12, h4("Genome annotation")),
                                  column(12, fileInput("protein_coding_genes", "Protein coding genes", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("Enhanser", "Enhancer", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("LincRNA", "Long noncoding RNAs", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, fileInput("Insulator", "Loop's insulator", multiple = TRUE,
                                          accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                  ),
                                  column(12, p(btn_run, align="center")),
                        )
              )
#            )
#          )
        ))
#    }
    out
#  })
  # required observers that update the color of the active button!
#  eval(parse(text=genObserver_menus(pat="btn_import_data_", n=1:6, updateVal="cur_selection_import")))
#  return(uiOutput("ui_sel_resbtns_import"))
})

# specific (gui)-options for tabs
output$ui_tabs <- renderUI({
  out <- fluidRow(
          column(width = 10,
             tabsetPanel(id = "tabs", type = "tabs",

                         tabPanel("Rearrangemet signal"),#, withSpinner(plotOutput("distPlot"))),
                         tabPanel("Table", withSpinner(tableOutput("table")))

                         )
            )
        )
  out
})

output$ui_import_data_main <- renderUI({
    out <- fluidRow(
    column(width = 10, uiOutput("ui_tabs") , class="tab"), #
    column(width = 10, offset = 0, h3("Uploading SVs dataset and genomic annotations for a cancer cohort"), class="wb-header"),
    column(width = 10, offset = 0, p("Load the requaried dataset to compute the BPpc and predict the potential driver elements ."), class="wb-header-hint"),
    )
    #out <- list(out, uiOutput("ui_tabs"))

})

output$ui_app <- renderUI({
  fluidRow(
    # Sidebar layout with input and output  ----
#      sidebarLayout(
      # Sidebar panel for inputs ----
          column(width = 2, uiOutput("ui_import_data_sidebar_left"), class="wb_sidebar"),
#          mainPanel(
            column(width = 10, uiOutput("ui_import_data_main"), class="wb-maincolumn")
#          )
#    )
  )
})

output$table <- renderTable({
  #BP_DT()
  #  significant_element_best_candidate_in_peak[, .(chrom, range_s, range_d, peak.ID, range_length, N_SV, cancer_type, file, GeneName)]
  #xtable
  #signal_BP_DT() #
  #signal_covariates_BP_DT()
  signal_BP_Covariates_DT()

})


