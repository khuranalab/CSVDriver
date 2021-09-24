#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Define server logic required to draw a histogram
CSVDriver_server <- shinyServer(function(input, output, session) {

#    options(shiny.maxRequestSize=100*1024^2)
    source('lib/ui_about.R', local = TRUE)
    source('lib/ui_app.R', local = TRUE)
    source('lib/reactive_lib.R', local = TRUE)



    output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')

    })

})
