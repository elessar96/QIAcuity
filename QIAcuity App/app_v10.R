

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
    shinyDirButton(id='directory_input', label='Select an input folder', title='Please select an input folder')),
    fluidRow(strong("Currently selected output path:")),
    fluidRow(verbatimTextOutput("indir", placeholder = TRUE)), 
    fluidRow(shinyDirButton(id='directory_output', label='Select an output folder', title='Please select an output folder')),
    fluidRow(strong("Currently selected output path:")),
    fluidRow(verbatimTextOutput("outdir", placeholder = TRUE)), 
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
  root <- getwd() %>% dirname() %>% dirname() %>% dirname()
  
  shinyDirChoose(input, id = 'directory_input', roots=c(wd=root))
  shinyDirChoose(input, id = 'directory_output', roots=c(wd=root))
  global <- reactiveValues(datapath = getwd(), status = "")
  
  output$indir <- renderText({
    input$directory_input %>% unlist() %>% paste(collapse="/") %>% gsub("wd", "", .) %>% paste(".", ., sep="") %>% gsub("^\\.\\d+$", "./", .)
  })
  
  observe(outdir <- reactive(input$directory_output))
  
  output$outdir <- renderText({
    input$directory_output %>% unlist() %>% paste(collapse="/") %>% gsub("wd", "", .) %>% paste(".", ., sep="") %>% gsub("^\\.\\d+$", "./", .)
  })
  
  output$done <- renderText({
    if(input$analyse==TRUE){
      root <- getwd() %>% dirname() %>% dirname() %>% dirname()
      
      input_path <- input$directory_input%>% unlist() %>% paste(collapse="/") %>% gsub("wd", "", .) %>% paste(root, ., sep="") # %>% normalizePath()
      output_path <- input$directory_output%>% unlist() %>% paste(collapse="/") %>% gsub("wd", "", .) %>% paste(root, ., sep="") # %>% normalizePath()
      
      if(length(input$coupled)==0){
        coupled <- data.frame(ch1=character(length=0), ch2=character(length=0))
      }else{
        coupled <- input$coupled %>% strsplit(., " ") %>% data.frame() %>% t() %>% data.frame()
        colnames(coupled) <- c("ch1", "ch2")
        rownames(coupled) <- NULL
      }

      setup()
      
      result <- QIAcuityAnalysis(input_path=input_path, output_path=output_path, coupled_channels = coupled)
      summarizeQIAcuityResult(result, detailed_plots=input$detailPlots, scatterplots_1d=input$scat1dPlots, output_path=output_path)
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