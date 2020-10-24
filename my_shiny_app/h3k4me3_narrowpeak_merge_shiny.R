suppressMessages({
library(shiny)
library(shinyFiles)
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(grid)
require(gridExtra)
library(data.table)
options(stringsAsFactors = FALSE)})
#source("myapp/h3k4me3_narrowpeak_mergeshiny.R")

## Merge h3k4me3 narrow peaks from all samples in each celltype. The union of broad peaks across samples are generated.
celltype_unionprofile <- function(peakdir, resdir, celltypes, qValue_cutoff ){
  extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric")
  celltypes2 <- as.list(read.csv(celltypes,header = FALSE))
  for (ct in celltypes2){
    p = paste0(peakdir,"/", ct)
    files = list.files(p)
    i=0
    if (length(files) != 0){
      x = GRanges()
      cat ("Processing cell type ", ct, "\n" )
      for (f in files){
        cat ("Sample: ", f[1])
        pks = import(paste0(peakdir,"/", ct, "/",f[1]), format = "bed", extraCols = extraCols_broadPeak)
        pks_cutoff = pks[elementMetadata(pks)$qValue  >= qValue_cutoff]
        cat("\t Number of peaks: ", length(pks_cutoff), "\n" )
        if (length(pks_cutoff) >= 10000){
          x = c(x, pks_cutoff) 
        }
      }
      if (length(x) !=0){
        x = reduce(x)
        export.bed(x, paste0(resdir, "/",  ct, "-union.bed"))
      }
    }
  }
}


# See above for the definitions of ui and server
ui <- fluidPage(

# App title ----
  titlePanel("Merging narrow H3K4Me3 peaks from multiple samples!"),
  shinyDirButton("dir", "Peak directory", "Upload"),
  verbatimTextOutput("dir", placeholder = TRUE),
  numericInput("num", "qValue_cutoff", value = 0.05),
  fileInput('cell', 'celltypes',accept=c('csv', 'comma-separated- values,text/plain', '.csv')),
  verbatimTextOutput("qvalue"),
  tableOutput(outputId = 'cell_csv'),
  tableOutput(outputId = 'union_bed')
  
)

server <- function(input, output) {
  # directory path for celltype_unionprofile function
  shinyDirChoose(
    input,
    'dir',
    roots = c(home = '~'),
    filetypes = c('', ".bed", "bed.gz"))

  global <- reactiveValues(datapath = getwd())
  dir <- reactive(input$dir)
  output$dir <- renderText({global$datapath})

  observeEvent(ignoreNULL = TRUE,
    eventExpr = {
        input$dir},
        handlerExpr = {
        if (!"path" %in% names(dir())) return()
        home <- normalizePath("~")
        global$datapath <-
        file.path(home, paste(unlist(dir()$path[-1]), collapse = .Platform$file.sep))
        })

  # qvalue input for celltype_unionprofile function
  output$qvalue <- renderText({ input$num })
  
  # celltypes input for celltype_unionprofile function from csv file
  output$cell_csv <- renderTable({
    inFile <- input$cell
        tbl <- read.csv(inFile$datapath)
        return(tbl)})

  # celltype_unionprofile function output
  bed_output <- reactive(celltype_unionprofile(global$datapath, global$datapath, input$cell$datapath, input$num))
  output$union_bed <- renderTable(bed_output())
}

shinyApp(ui = ui, server = server)