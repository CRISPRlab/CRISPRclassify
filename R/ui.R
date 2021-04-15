presentModal <- function(title, msg){
  showModal(modalDialog(
    title = title,
    msg,
    easyClose = TRUE,
    footer = modalButton("Dismiss")
  ))
}

getCSS <- function(){
  return("
.fileSelectDiv {
  margin-top: 10px;
  background-color: #f5f5f5;
  padding: 25px;
  border-radius: 4px;
  box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
  webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
}
.instructionsCol {
  padding-left: 0px;
}
.instructionsDiv {
  text-align: center;
  margin-top:72px;
  background-color: #f5f5f5;
  padding:25px;
  border-radius: 4px;
  box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
  webkit-box-shadow: inset 0 1px 1px rgba(0,0,0,.05);
  font-size: 19px;
  height: 188px;
}
.resultsDiv {
  margin-top: 25px;
  margin-bottom: 10px;
}
.shiny-input-container:not(.shiny-input-container-inline) {
  width: auto;
}
.nowrap {
  white-space: nowrap;
}
.shiny-notification {
  position: fixed; top: 5px ;right: 13px;
  min-width: 350px;
}
.logo {
  position: relative;
  left: -5px;
  top: 5px;
}
")
}

getShinyUI <- function(css){
  return(fluidPage(

    fluidRow(
      column(4,
             tags$head(tags$style(HTML(css))),
             img(src ="assets/logo.png", height=62,width=300) %>% tagAppendAttributes(class = 'logo'),
             tags$div(
               fileInput("inFile", h4("Choose Genomic .FASTA File"),
                         multiple = FALSE,
                         accept = c(".fasta",
                                    ".fna",
                                    ".fa")),
               uiOutput("dynamicUI"),
               br(),
               uiOutput("dynamicUIDownload"),
             ) %>% tagAppendAttributes(class = 'fileSelectDiv')
      ),
      column(8, tags$div(
        tags$span(
          "Welcome to CRISPRclassify! To get started, choose any assembled genome or metagenome file for analysis.
              Once it has been uploaded, the Classify button will appear. Click it to begin analysis! Results will appear
              below once the analysis completes, along with a Download button to save your work. For more details, check out
              "
        ),
        tags$a(href="https://github.com/CRISPRlab/CRISPRclassify", "GitHub.")
      ) %>% tagAppendAttributes(class = 'instructionsDiv')

      ) %>% tagAppendAttributes(class = 'instructionsCol')
    ),
    fluidRow(

      column(12,
             tags$div(
               plotOutput("cas_dist_plot"),
               DT::dataTableOutput("classified_results")
             ) %>% tagAppendAttributes(class = 'resultsDiv')
      )
    )

  ))
}

getFeatureImportance <- function(type){
  feat_imp_file <- CRISPRclassify:::getFeatImpFile(CRISPRclassify:::getWrkDir())
  read_csv(feat_imp_file, col_names = F) %>%
    filter(X5 == type,
           !X1 %in% c("gc", "palIdx", "length")) %>%
    group_by(X5) %>%
    mutate(norm_gain = as.numeric(scale(X2))) %>%
    ungroup() %>%
    filter(norm_gain > 0.675) %>% pull(X1) #75th percentile
}

ngramToHtmlBlue <- function(ngram) {
  #c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
  str_c('<span class="no-indent" style="color:Blue;font-weight:bold;">',ngram,'</span>', collapse = '')
}
ngramToHtmlOrange<- function(ngram) {
  #c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
  str_c('<span class="no-indent" style="color:Orange;font-weight:bold;">',ngram,'</span>', collapse = '')
}

getSyntaxHighlightedSequence <- function(seq, type){
  getFeatureImportanceCache <- memoise(getFeatureImportance)
  feat_vect <- getFeatureImportanceCache(type)
  formatted_seq <- seq
  for (feat in feat_vect) {
    formatted_seq <- str_replace_all(formatted_seq,CRISPRclassify:::getReverseCompliment(feat), ngramToHtmlBlue)
    formatted_seq <- str_replace_all(formatted_seq,feat, ngramToHtmlOrange)
  }

  return(formatted_seq)
}
