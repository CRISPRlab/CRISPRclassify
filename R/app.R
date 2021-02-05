
## -------------------------- ##

#' @examples
#' launchApp()
#' @export
launchApp <- function() {
  require(shiny)
  require(dplyr)
  require(stringr)
  require(magrittr)
  require(readr)
  require(xgboost)
  require(DT)
  require(memoise)
  require(tidyr)
  require(stringi)
  require(stringdist)
  require(ggplot2)
  require(fs)
  #TODO these require statements might not be needed after qualifying all function names

  wrkDir <- fs::path_package("CRISPRclassify")
  print(str_c("working directory is ", wrkDir))
  selectedFullModel <- '1'
  selectedSeqModel <- '1'

  options(java.parameters = "-Xmx8G")
  GB_size <- 200
  options(shiny.maxRequestSize=GB_size*1000*1024^2)
  ## Manage shiny script permissions ##
  system(str_c("chmod +x ",wrkDir,"/lib/crVSkeleton.sh"))
  system(str_c("chmod +x ",wrkDir,"/lib/java/minced/minced"))

  css <- "
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
"
getFasta <- function(file_path){
  read_csv(file_path, col_names = F) %>%
    mutate(X2 = lead(X1)) %>%
    dplyr::rename(id = X1, sequence = X2) %>%
    filter(!str_detect(sequence, ">")) %>%
    mutate(id = str_replace_all(id, ">", "")) %>% return()
}
getSpacerId <- function(seq){
  parsed_vect <- seq %>% str_replace_all( c("CRISPR_" = "|", "_spacer_" = "|")) %>% str_split("\\|",  simplify = T)
  return(parsed_vect[1,3])
}
getSpacerLocus <- function(seq){
  parsed_vect <- seq %>% str_replace_all( c("CRISPR_" = "|", "_spacer_" = "|")) %>% str_split("\\|",  simplify = T)
  return(parsed_vect[1,2])
}
getRepeatId <- function(seq){
  parsed_vect <- seq %>% str_replace_all( c("locus" = "|", "\\_" = "|")) %>% str_split("\\|",  simplify = T)
  parsed_vect[1,4]
}
getRepeatLocus <- function(seq){
  parsed_vect <- seq %>% str_replace_all( c("locus" = "|", "\\_" = "|")) %>% str_split("\\|",  simplify = T)
  parsed_vect[1,3]
}
getGC <- function(seq){
  return((length((str_extract_all(seq, "G"))[[1]]) + length((str_extract_all(seq, "C"))[[1]]))/nchar(seq)*100)
}
getLength <- function(seq){
  return(nchar(seq))
}
getPalIdx <- function(seq){
  sum(str_split(seq, "")[[1]] == str_split(stri_reverse(seq),"")[[1]])/nchar(seq)
}

getContigName <- function(id){
  return(str_split(id, "_CRISPR")[[1]][1])
}

getContigLabels <- function(dat_vect){
  contig_name <- NA
  contig_name_vect <- as.character()
  for (dat in dat_vect){
    if(str_detect(dat, "Sequence")){
      contig_name <- dat
    }
    contig_name_vect <- c(contig_name_vect, contig_name)
  }
  return(contig_name_vect)
}

getLocusNamesLabels <- function(dat_vect){
  locus_name <- NA
  locus_name_vect <- as.character()
  for (dat in dat_vect){
    if(str_detect(dat, "CRISPR")){
      locus_name <- dat
    }
    locus_name_vect <- c(locus_name_vect, locus_name)
  }
  return(locus_name_vect)
}

getAllFastaDtl <- function(datapath){

  crispr_file = paste(datapath, ".crisprs", sep="")

  crispr_df <- read_csv(crispr_file, col_names = F) %>%  rename(dat = 1)
  all_info_fasta_df <- tibble(
    contig = getContigLabels(crispr_df$dat),
    locus = getLocusNamesLabels(crispr_df$dat),
    record = crispr_df$dat) %>%
    filter(!str_detect(record, "POSITION|CRISPR|Sequence|-----|Average|Time")) %>%
    mutate(contig = str_replace_all(contig, "Sequence |'", "")) %>%
    separate(contig, c("contig_name", "contig_id"), sep = " \\(") %>%
    separate(record, c("position", NA, "repeat_sequence", "spacer_sequence", NA), sep = "\t")  %>%
    mutate(contig_id = trimws(str_replace_all(contig_id," bp\\)", "")),
           locus_id = trimws(str_sub(locus, 7, 11))) %>%
    separate(locus, c(NA, "range"), sep = "Range: ") %>%
    mutate(range = str_replace_all(range, "\\s" , ""))  #%>% group_by(contig_name, contig_id, range, repeat_sequence, locus_id) %>% tally() %>% rename(spacer_n = n) %>% ungroup()

  return(all_info_fasta_df)

}

getRepeatDf <- function(repeat_file, datapath){
  repeat_df <- getFasta(repeat_file)
  repeat_locus_vect <- sapply(repeat_df$id, getRepeatLocus)
  repeat_id_vect <- sapply(repeat_df$id, getRepeatId)

  repeat_gc_vect <- sapply(repeat_df$sequence, getGC)
  repeat_palidx_vect <- sapply(repeat_df$sequence, getPalIdx)
  repeat_length_vect <- sapply(repeat_df$sequence, getLength)

  repeat_expanded_df <- repeat_df %>% mutate(
    locus = repeat_locus_vect,
    id = repeat_id_vect,
    gc = repeat_gc_vect,
    palIdx = repeat_palidx_vect,
    length = repeat_length_vect)

  all_fasta_dtl_df <- getAllFastaDtl(datapath) %>% dplyr::select(repeat_sequence, contig_name, locus_id) %>% distinct()

  repeat_expanded_df <-  repeat_expanded_df %>% distinct() %>%  inner_join(all_fasta_dtl_df, by = c("sequence" = "repeat_sequence", "locus" = "locus_id")) #Note: No contig name to join on
  repeat_expanded_df <- repeat_expanded_df %>%
    group_by(locus, contig_name) %>%
    mutate(repeat_count = n(), average_repeat_length = mean(length)) %>%
    ungroup()  %>%
    group_by(locus, contig_name, sequence) %>%
    mutate(distinct_repeat_count = n()) %>%
    ungroup()  %>%
    group_by(locus,contig_name, gc, palIdx, length, sequence, repeat_count, average_repeat_length, distinct_repeat_count) %>%
    summarize(id = max(id)) %>% ungroup()
  return(repeat_expanded_df)
}

getSpacerDf <- function(spacer_file, datapath){
  spacer_df <- getFasta(spacer_file)
  contig_name_vect =sapply(spacer_df$id, getContigName)
  spacer_locus_vect <- sapply(spacer_df$id, getSpacerLocus)
  spacer_id_vect <- sapply(spacer_df$id, getSpacerId)
  spacer_length_vect <- sapply(spacer_df$sequence, getLength)

  spacer_df <- spacer_df %>% mutate(locus = spacer_locus_vect,
                                    #                                    id = spacer_id_vect,
                                    length = spacer_length_vect,
                                    contig_name = contig_name_vect) %>%
    select(-id) %>%
    distinct()


  all_fasta_dtl_df <- getAllFastaDtl(datapath) %>% select(spacer_sequence, contig_name, locus_id,range) %>% distinct()

  spacer_df <- spacer_df %>% distinct() %>%  inner_join(all_fasta_dtl_df, by = c("sequence" = "spacer_sequence" , "contig_name", "locus" = "locus_id"))

  spacer_df <- spacer_df %>% group_by(locus, contig_name, range) %>%
    summarise(spacer_n = n(), avg_spacer_length = mean(length)) %>%
    ungroup()
  return(spacer_df)
}
# getSequenceAnalyticDf <- function(crispr_file = "data/meta.fasta.crisprs", spacer_file = "data/meta.fasta_spacers.fa", repeat_file = "data/meta.fasta_repeats.fa"){



getSequenceAnalyticDf <- function(datapath){

  spacer_file = paste(datapath, "_spacers.fa", sep="")
  repeat_file = paste(datapath, "_repeats.fa", sep="")

  all_fasta_df <- getAllFastaDtl(datapath)

  #get spacer count and average spacer length
  #getting repeat features
  repeat_df <- getRepeatDf(repeat_file, datapath)
  spacer_df <- getSpacerDf(spacer_file, datapath)

  all_fasta_df %>% left_join(repeat_df)

  sequence_analytic_df <- repeat_df %>% inner_join(spacer_df, by = c("locus", "contig_name")) %>%
    mutate(repeat_spacer_ratio = average_repeat_length/avg_spacer_length,
           unique_id = str_c(locus, "_", id, "_", contig_name)) %>%
    rename(seq = sequence )

  return(sequence_analytic_df)
}
#repeat counts should be count of identical repeats or count of repeats in locus (even if different)



xg_matrix_format_file <- str_c(wrkDir,"/csvRef/matrix_format_2020-12-04_nchar_5.csv")
feat_imp_file <- str_c(wrkDir,"/csvRef/xg_feature_importance_5_2020-12-04.csv")
cas_vect <- c('I-A','I-B','I-C','I-D','I-E','I-F','I-U','II-A','II-B','II-C','III-A','III-B','III-C','III-D','III-E','III-F','IV-A','IV-C','V-A','V-B1','V-B2','V-F','V-J','V-U1','V-U2','V-U4','VI-A','VI-B','VI-C','VI-D')
n_char_cnt <- 5

getCompliment <- function(nuc){
  c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
}

#reverses string and applies compliment
getReverseCompliment <- function(sequence){
  rc_vect <- as.vector(character())
  for (nuc in str_split(stri_reverse(toupper(sequence)) ,'')[[1]]) {
    rc_vect <- c(rc_vect, getCompliment(nuc))
  }
  return(str_c(rc_vect, collapse = ''))
}

getFeatureImportance <- function(type){
  read_csv(feat_imp_file, col_names = F) %>%
    filter(X5 == type,
           !X1 %in% c("gc", "palIdx", "length")) %>%
    group_by(X5) %>%
    mutate(norm_gain = as.numeric(scale(X2))) %>%
    ungroup() %>%
    filter(norm_gain > 0.675) %>% pull(X1) #75th percentile
}
getFeatureImportanceCache <- memoise(getFeatureImportance)


ngramToHtmlBlue <- function(ngram) {
  #c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
  str_c('<span class="no-indent" style="color:Blue;font-weight:bold;">',ngram,'</span>', collapse = '')
}
ngramToHtmlOrange<- function(ngram) {
  #c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
  str_c('<span class="no-indent" style="color:Orange;font-weight:bold;">',ngram,'</span>', collapse = '')
}


getSyntaxHighlightedSequence <- function(seq, type){
  feat_vect <- getFeatureImportanceCache(type)
  # print(feat_vect)
  formatted_seq <- seq
  for (feat in feat_vect) {
    formatted_seq <- str_replace_all(formatted_seq,getReverseCompliment(feat), ngramToHtmlBlue)
    formatted_seq <- str_replace_all(formatted_seq,feat, ngramToHtmlOrange)
  }

  return(formatted_seq)
}


get_tokens <- function(seq, id) {
  token_vect <- as.vector(as.character())
  n_char_for_iter <- n_char_cnt - 1

  rc_seq <- getReverseCompliment(seq)

  for (i in 1:nchar(seq)) {
    if (i + n_char_for_iter <= nchar(seq)) {
      token <- str_sub(seq, i, i + n_char_for_iter)
      token_vect <- c(token, token_vect)
    }
  }

  for (i in 1:nchar(rc_seq)) {
    if (i + n_char_for_iter <= nchar(rc_seq)) {
      token <- str_sub(rc_seq, i, i + n_char_for_iter)
      token_vect <- c(token, token_vect)
    }
  }

  token_agg <- table(token_vect)
  return(tibble(
    id = id,
    token_count = token_agg,
    feature = names(token_agg)
  ))
}
getModelFilePath <- function(cas) {
  return(str_c(wrkDir,"/models/",cas,"_2020-12-04_nchar_5"))
}
getXGModel <- function(cas){
  return(xgb.load(getModelFilePath(cas)))
}
getXGModelCache <- getXGModel

getTokenDF <- function(sequence_analytic_df){

  token_df <- do.call("rbind", mapply(get_tokens, sequence_analytic_df$seq, sequence_analytic_df$unique_id, SIMPLIFY = F)) %>%
    rename(value = token_count) %>%
    mutate(value = as.numeric(value)) %>%
    dplyr::select(id, feature, value) %>%
    rename(unique_id = id) %>%
    pivot_wider(names_from = feature, values_from = value, values_fill = 0)
  return(token_df)
}

strainLkp <- function(repeatSeq, masterLkpDf){
  idx <- amatch(repeatSeq, masterLkpDf$repeat_seq, nomatch=(nrow(masterLkpDf)), maxDist = 10)
  org <- masterLkpDf[idx,]
  return(org)
}

getXGBoostPredictions <- function(sequence_analytic_df, progress){
  print("...Getting xgboost predictions...")
  all_feat_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% dplyr::select(X3) %>% distinct() %>% pull()
  token_df <- getTokenDF(sequence_analytic_df)
  bio_feat_df <- sequence_analytic_df %>% inner_join(token_df, by = "unique_id")

  missing_col_vect <- setdiff(all_feat_vect , colnames(bio_feat_df))
  missing_feature_df <- matrix( rep( 0, len=length(missing_col_vect)*nrow(bio_feat_df)), nrow = nrow(bio_feat_df)) %>% as.data.frame()
  names(missing_feature_df) <- missing_col_vect

  bio_feat_df <- bio_feat_df %>% bind_cols(missing_feature_df)
  #bio_feat_df %>% write_csv("seq_analytic_df_audit.csv")

  result_df <- tibble(unique_id = as.character(), seq = as.character(), cas_type = as.character(), probability = as.double())

  # cas <- 'I-E'
  allCasN <- length(cas_vect)
  adjustmentInc <- (allCasN + (.1 * allCasN)) ## since we start on .1, account for that ##
  for (cas in cas_vect){

    print(str_c("running: ", cas))
    progress$inc(1/adjustmentInc, message = str_c("Checking ", cas, "."))

    if (file.exists(getModelFilePath(cas))){
      cas_feature_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% filter(X1 == cas) %>% dplyr::select(X3) %>% pull()
      feature_df <- bio_feat_df %>% select_at(cas_feature_vect) #%>% dplyr::select(-unique_id)
      feature_mat <- feature_df %>%  as.matrix()
      bstSparse <- getXGModelCache(cas)
      feature_sparse_mat <- as(feature_mat, "dgCMatrix")

      if(nrow(feature_sparse_mat) > 0){

        pred <- predict(bstSparse, feature_sparse_mat)
        temp_result_df <- sequence_analytic_df %>% bind_cols(probability = pred) %>% dplyr::select(contig_name, unique_id, locus, id, seq, probability, range, distinct_repeat_count) %>% distinct() %>% mutate(cas_type = cas)

        result_df <- result_df %>% bind_rows(temp_result_df)

      } else {
        print(str_c("predictions not possible for subtype ", cas))
      }

    }

    else {print("file doesnt exist")}
  }

  progress$set(message = "Looking up closest strains. ", detail = 'Just a moment...', value = .98)
  # result_df %>% write_csv("_tmp_full_output.csv")
  max_result_df <- result_df %>% group_by(unique_id) %>% filter(probability == max(probability)) %>% ungroup() %>%  select(seq,cas_type, distinct_repeat_count , probability,contig_name, range, locus)
  lkpMaster <- read_csv(str_c(wrkDir, '/csvRef/master_repeat_lkp.csv'))
  max_result_df_lkp <- transform(max_result_df, closestStrain = strainLkp(seq, lkpMaster))
  progress$set(message = "Wrapping up. ", detail = 'Just a moment...', value = 1)
  return(max_result_df_lkp)
}



## ------------------------- SHINY APP BEGIN --------------------- ##
  shinyApp(


ui <- fluidPage(

  fluidRow(
    column(4,
           tags$head(tags$style(HTML(css))),
           img(src ="assets/logo.png", height=62,width=300) %>% tagAppendAttributes(class = 'logo'),
           tags$div(
             fileInput("inFile", h4("Choose Genomic .FASTA File"),
                              multiple = FALSE,
                              accept = c("text/comma-separated-values,text/plain",
                                         ".fasta",
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

),

server <- function(input, output){

  values <- reactiveValues()
  values$fileStatus <- 'empty'

  # onclick("inFile", resetUI(values))

  observeEvent(input$inFile, {
    # print('Testing upload')
    on.exit(values$fileStatus <- 'uploaded')
  })

  ## Load Classification Model based on Selection Panel ##
  observeEvent(input$fullFileModels, {
    # load model based on radio selection
    print(input$fullFileModels)
    selectedFullModel <- input$fullFileModels
  })

  observeEvent(input$sequenceModels, {
    # load model based on radio selection
    print(input$sequenceModels)
    selectedSeqModel <- input$sequenceModels
  })


  ## ----------- Manage Classify Btn click ---------- ##

  ## << Genome File Upload >> ##
  parseResults <- eventReactive(input$fileUploadBtn,{
    if (is.null(input$inFile)){
      showNotification("Please select a file...")
      return(NULL)
    }

    on.exit(values$fileStatus <- 'empty')

  progress <- shiny::Progress$new()
  on.exit(progress$close(), add=TRUE)

  # Extract repeats from .crisprs file
  print(input$intFile$datapath)

  progress$set(message = "Finding CRISPR repeats.", detail = 'This may take a few moments...', value = .05)
  system(sprintf("%s/lib/crVSkeleton.sh -c -f %s",wrkDir, input$inFile$datapath), intern = TRUE) #TRUE

  progress$set(message = "Tokenizing repeats.", detail = 'This may take a few moments...', value = .1)

  analytic_df <- getSequenceAnalyticDf(input$inFile$datapath)
  progress$set(message = "Predicting subtypes.", detail = "This may take a few moments...",  value = .1)
  getXGBoostPredictions(analytic_df, progress) %>%
    mutate(cas_type = case_when(cas_type == 'V-J' ~ 'V-K', cas_type == "I-U" ~ "I-G", T ~ cas_type)) %>%
    return()

  })


  output$classified_results <- DT::renderDataTable({
    parseResults() %>%
      rowwise() %>%
      mutate(
             seq = getSyntaxHighlightedSequence(seq, cas_type),
             probability = round(probability, 3)) %>%
      ungroup() %>%
      mutate(cas_type = as.factor(cas_type)) %>%
     rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Contig` = contig_name, `Range` = range, `Locus` = locus, `Repeat Count` = distinct_repeat_count, `Closest Strain` = closestStrain.org_desc) %>%
      arrange(Locus) %>%
      dplyr::select(Locus, `Contig`,  Range, `Subtype`, `Repeat`, `Repeat Count`, `Probability`, `Closest Strain`)
  }, filter = 'top', escape = FALSE, options = list(
    columnDefs = list(
      list(className = "nowrap", targets = "_all")
    )
  ))

  output$cas_dist_plot <- renderPlot({
    parseResults() %>%
      rename(Subtype = cas_type) %>%
      group_by(Subtype) %>%
      select(Subtype, locus) %>% distinct() %>%
      tally() %>%
      ggplot(aes(x = n, y = Subtype, fill=Subtype)) + geom_histogram(stat = 'identity') + xlab("Locus Count by Subtype") + ylab("Subtype")
  })

  output$downloadData <- downloadHandler(

    filename = function() {
      paste("subtype_detail.csv", sep = "")
    },
    content = function(file) {
      parseResults() %>%
        mutate(
          probability = round(probability, 6)) %>%
        ungroup() %>%
        mutate(cas_type = as.factor(str_sub(cas_type, 1,15))) %>%
        rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Contig Name` = contig_name, `Range` = range, `Locus` = locus, `Repeat Count` = distinct_repeat_count, `Closest Strain` = closestStrain.org_desc, `Accession` = closestStrain.accession, `Location` = closestStrain.location, `Type` = closestStrain.type, `Closest Repeat` = closestStrain.repeat_seq) %>%
        arrange(Locus) %>%
        dplyr::select(Locus, `Contig Name`,  Range, `Subtype`, `Repeat`, `Repeat Count`, `Probability`, `Closest Strain`, `Accession`, `Location`, `Type`, `Closest Repeat`) %>%
      write.csv(file, row.names = FALSE)
    }
  )



  output$dynamicUI <- renderUI({
    if(values$fileStatus == 'uploaded'){
      tagList(
        actionButton(inputId = "fileUploadBtn", label = "Classify")

      )
    }

  })

    output$dynamicUIDownload <- renderUI({
      #if (input$file_choice == ".fasta"){
      if(!is.null(parseResults())){
        tagList(
          downloadButton("downloadData", "Download")

        )
      }
  })


}


)}

