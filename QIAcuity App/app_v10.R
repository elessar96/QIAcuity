
options(RCurlOptions=list(proxy="www-int2.inet.dkfz-heidelberg.de:80", http.version=1))

Renviron <- readLines('~/.Renviron')

http_pos <-  grep("^http_proxy", Renviron)
edit_necessary <- FALSE

if(!Renviron[[http_pos]] == "http_proxy=http://www-int2.dkfz-heidelberg.de:3128/"){
  Renviron[[http_pos]] <- "http_proxy=http://www-int2.dkfz-heidelberg.de:3128/"
  edit_necessary <- TRUE
}

https_pos <- grep("^https_proxy", Renviron)

if(!Renviron[[https_pos]] == "https_proxy=http://www-int2.dkfz-heidelberg.de:3128/"){
  Renviron[[https_pos]] <- "https_proxy=http://www-int2.dkfz-heidelberg.de:3128/"
  edit_necessary <- TRUE
}

if(edit_necessary==TRUE){
  #writeLines(Renviron, con = "~/.Renviron")
  print(".Renviron file needs to be altered to include DKFZ proxy server info. App restart required.")
  stopApp()
}

if (!require("shiny", quietly = TRUE))
  install.packages("shiny")
  library(shiny)

if (!require("shinyFiles", quietly = TRUE))
  install.packages("shinyFiles")
  library(shinyFiles)

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
  library(dplyr)

if (!require("bslib", quietly = TRUE))
  install.packages("bslib")
  library(bslib)

if (!require("shinycssloaders", quietly = TRUE))
  install.packages("shinycssloaders")
  library(shinycssloaders)

source("../QIAcuityAnalysis/R/QIAcuity_package_V3.R")

ui <- fluidPage( # Application title
  titlePanel("QIAcuity data analysis"),
  mainPanel(
    
    fluidRow(
    p("Please paste the path to a folder which contains raw data exported from the QIAcuity software."),
    textInput("dir", "Path:", value=getwd()),
    actionButton("setdir", "Select")),
    fluidRow(strong("Currently selected path:")),
    fluidRow(verbatimTextOutput("dir", placeholder = TRUE)), 
    strong("Plotting options:"),
    checkboxInput("detailPlots", label="Generate 2D scatterplots"),
    checkboxInput("scat1dPlots", label="Generate 1D scatterplots"),
    strong("Coupled channels:"),
    checkboxGroupInput("coupled",
                       label = "Select coupled channels (targets are on same chromosome or are othwerwise expected to be in the same partitions):",
                       choices = c("C R", "C Y", "C O", "C G", "R Y", "R G", "R O", "Y G", "Y O", "O G")),
    fluidRow(
      actionButton("analyse", "Start analysis"),
      actionButton('stop', 'Exit')),
    shinycssloaders::withSpinner(verbatimTextOutput("done", placeholder = TRUE), type=5, proxy.height =50, size=0.5)
  ),
  theme = bs_theme(version = 4, bootswatch = "darkly",
                   bg = "#dcdcdc", fg = "#003fa4", primary = "#4875c3", secondary ="#f55000", success="#c9c9c9",
                   base_font = font_google("Open Sans"),
                   code_font = font_google("Open Sans")))

server <- function(input, output) {
  
  global <- reactiveValues(datapath = getwd(), status = "")
  
  dir <- reactive(input$dir)
  
  output$dir <- renderText({
    global$datapath
  })
  
  output$done <- renderText({
    if(input$analyse==TRUE){
      path <- global$datapath %>% normalizePath()
      if(length(input$coupled)==0){
        coupled <- data.frame(ch1=character(length=0), ch2=character(length=0))
      }else{
        coupled <- input$coupled %>% strsplit(., " ") %>% data.frame() %>% t() %>% data.frame()
        colnames(coupled) <- c("ch1", "ch2")
        rownames(coupled) <- NULL
      }
      

      setup()
      result <- QIAcuityAnalysis(path, coupled_channels = coupled)
      summarizeQIAcuityResult(result, detailed_plots=input$detailPlots, scatterplots_1d=input$scat1dPlots)
      "Analysis completed"
    }
    
  })
  
  observeEvent(
               eventExpr = {
                 input$setdir
               },
               handlerExpr = {
                 print("Select")
                 global$status <- "Ready"
                 global$datapath <- input$dir
               })
  
  observeEvent(input$stop,{
    print("App closed")
    stopApp()
  })

}

# Run the application
shinyApp(ui = ui, server = server)