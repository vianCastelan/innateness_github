# Libraries -------------------------------------------------------------------------

# load libraries
library(tidyverse)
library(ggplot2)
library(shiny)
library(shinyjs)
library(glue)
library(data.table)
library(reticulate)

# Data ------------------------------------------------------------------------------

source("R/load_data.R")

# Parameters ------------------------------------------------------------------------

## default gene
one_gene_symbol = "Tbx21"
## use python
# use_python("../.venv/bin") ## This should be set in .Rprofile
## get python function from script
source_python("python/innate_score.py") ## long load time, optimize....

# Functions -------------------------------------------------------------------------

source("R/plot_boxplot.R")
source("R/save_figure.R")
source("R/optimize_png.R")

which_numeric_cols <- function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}


# User interface -------------------------------------------------------------------
ui <- fluidPage(
    useShinyjs(),
    ## tags 
    tags$head(
      tags$link(
        rel = "stylesheet", type = "text/css", href = "app.css"
      ),
    ),
    ##------------------------------------------------ Application layout in ui_tab.R
    navbarPage(
        title = "Innateness Dashboard 1.0v",
        ##------------------------------------------------------- layout in ui script
        source(file.path("R","ui_tab.R"), local=TRUE)$value
    ),
    HTML(
    "<footer class='myfooter page-footer'>
      <div class='text-center'>
      This website was created by <a href='https://www.github.com/gabrielascui'>Gabriel Ascui</a>
      </div>
    </footer>"
    )
)

# Server -------------------------------------------------------------------------
server <- function(input, output, session) {
    ## Innateness Beta interfase -------------------------------------------------
    this_gene <- one_gene_symbol
    up_table <- NULL
    ### Dynamic table ------------------------------------------------------------
    output$beta_table <- DT::renderDataTable({
        beta_table <- beta %>% 
            arrange(pval) %>%
            head(nrow(beta))
        rownames(beta_table) <- 1:nrow(beta_table)
        beta_table_genes <<- beta_table$gene
        numeric_cols <- colnames(beta_table)[which_numeric_cols(beta_table)]
        # Javascript-enabled table
        DT::datatable(
        data = beta_table,
        options = list(
            "lengthChange" = FALSE,
            "orderClasses" = TRUE
        ),
        selection = "single"
        ) %>%
        DT::formatSignif(columns = numeric_cols, digits = 2)
    }, server = TRUE)
    ### Dynamic graph per gene selected in table --------------------------------
    output$rnaseq_one_gene <- renderText({
        beta_table_rowid <- input$beta_table_rows_selected
        if (length(beta_table_rowid)) {
            this_gene <- beta_table_genes[beta_table_rowid]
        }
        # Querry
        if (this_gene %in% beta$gene) {
            plot_boxplot(this_gene)
        }
        retval <- "<div></div>"
        if (this_gene %in% beta$gene) {
          retval <- save_figure(
            filename = glue("rnaseq_boxplot_{marker}.png", marker = this_gene),
            width = 6, height = 5, dpi = 300,
            html_alt = this_gene,
            ggplot_function = function() { plot_boxplot(this_gene) }
          )
        }
        retval
    })
    ## Innateness table processing ----------------------------------------------- 
    tableInput <- reactive({
      req(input$file1)
      # Determine the original encoding of the file and convert to UTF-8
      encoding <- readr::guess_encoding(input$file1$datapath)$encoding[1]
      file_contents <- readLines(input$file1$datapath, encoding = encoding)
      user_data <- fread(text = iconv(file_contents, from = encoding, to = "UTF-8"), encoding = "UTF-8")
      user_data <- as.data.frame(user_data)
      ## process file 
      file_content <- run_innate(user_data)
    })

    ## Generate the CSV file link
    output$download_output <- downloadHandler(
      filename = function() {
        paste("measurements_", Sys.Date(),".csv", sep = "")
      },
      content = function(file) {
        write.csv(tableInput(),file, row.names = FALSE)
      }
    )

    output$done_data <- DT::renderDataTable({
      # Javascript-enabled table
      DT::datatable(
      data = tableInput(),
      options = list(
          "lengthChange" = FALSE,
          "orderClasses" = TRUE
      )
      )
    }, server = TRUE)

    ## Generate download link for the Python script
    output$downloadPYTHON <- downloadHandler(
      filename = function() {
        "innate_score.py" 
      },
      content = function(file) {
        # Read the content of the Python script file
        script_content <- readLines("python/innate_score.py")
        # Write the content to the file for download
        writeLines(script_content, file)
    })

    ## download example matrix 

}


# Launch App ---------------------------------------------------------------------

shinyApp(ui = ui, server = server)

# shiny::runApp("shiny_app/app.R")