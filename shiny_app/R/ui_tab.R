tabPanel(
    title = "Home",
    value = "data",

    fluidPage(
        HTML("
<h1 class='display-4 font-weight-normal'>Innateness of mouse T cells</h1>
<p class='lead font-weight-normal'>Gene expression on lymphocyte innateness gradient.</p>
    "),
    hr(),
    fluidRow(
        column(
            width = 6,
            h3("ImmGEN ULI RNA-seq with sorted splenocyte population"),
            htmlOutput("rnaseq_one_gene", style = "min-height:450px;")
        ),
        column(
            width = 6,
            h3("beta levels of association of genes to innateness gradient"),
            DT::dataTableOutput("beta_table"), 
            style = "min-height: 674px:"
        )
    ),
    br(),
    hr(),
    div(id = "geneinfo"),

    fluidRow(
        column(
            width = 12,
            h3("Use your own file (size restriction, 50 MB) to calculate relative innateness gradients"),
            ## make file interactive part here
            fileInput("file1", 
                      "Upload countmatrix .CSV file",
                      accept = c("text/txt","text/comma-separated-values","text/plain",".txt",".csv")
            ),
            hr(),
            # fileInput("file2", 
            #           "Upload col_data .CSV file",
            #           accept = c("text/txt","text/comma-separated-values","text/plain",".txt",".csv")
            # ),
            DT::dataTableOutput("done_data"),
            hr(),
            downloadButton("download_output", "Download Innateness Scores"),
            br(),
            downloadButton("downloadPYTHON","Download Python script")
        )
    ),
    hr(),
    HTML(
        "
        <h2>Read the documentation</h2>
        <div class='mypaper'>
        <h3>
        </h3>
        </div>
        "
    )
    )
)