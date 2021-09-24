#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
CSVDriver_ui <- shinyUI(
    fluidPage(

            navbarPage("CSVDriver", id="nav.main", #theme = paste0("sdcwww/", getShinyOption(".guitheme")),
                       tabPanel("About/Help", id = "landing_t", uiOutput("ui.about")),
                       tabPanel("App", id = "app_t", uiOutput("ui_app"))#,
                       #tabPanel("Anonymize", uiOutput("ui_anonymize")),
                       #tabPanel("Risk/Utility", uiOutput("ui_results")),
                       #tabPanel("Export Data", uiOutput("ui_export")),
                       #tabPanel("Reproducibility", uiOutput("ui_script")),
                       #tabPanel("Undo", uiOutput("ui_undo")),

                       #tags$head(tags$script(
                        #   src = paste0("sdcwww/", getShinyOption(".guijsfile"))
                       #))
            )
        # # Application title
        # titlePanel("Old Faithful Geyser Data"),
        #
        # # Sidebar with a slider input for number of bins
        # sidebarLayout(
        #     sidebarPanel(
        #         sliderInput("bins",
        #                     "Number of bins:",
        #                     min = 1,
        #                     max = 50,
        #                     value = 30)
        #     ),
        #
        # # Show a plot of the generated distribution
        #     mainPanel(
        #         plotOutput("distPlot")
        #     )
        # )

    )
)
