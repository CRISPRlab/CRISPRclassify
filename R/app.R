
## -------------------------- ##

#' @examples
#' launchApp()
#' @export
launchApp <- function() {

  ## Java check ##
  CRISPRclassify:::validateJava()

  ## Load HC vars ##
  wrkDir <- CRISPRclassify:::getWrkDir()
  xg_matrix_format_file <- CRISPRclassify:::getMatrixFormatFile(wrkDir)
  cas_vect <- CRISPRclassify:::getCasVect()

  print(str_c("Package working directory is ", wrkDir))
  selectedFullModel <- '1'
  selectedSeqModel <- '1'

  ## Manage jvm params ##
  options(java.parameters = "-Xmx8G")
  GB_size <- 200
  options(shiny.maxRequestSize=GB_size*1000*1024^2)

  ## Manage shiny script permissions ##
  system(str_c("chmod +x ",wrkDir,"/lib/crVSkeleton.sh"))
  system(str_c("chmod +x ",wrkDir,"/lib/java/minced/minced"))

  ## Load CSS and UI methods ##
  css <- CRISPRclassify:::getCSS()
  getXGModelCache <- CRISPRclassify:::getXGModel

## ------------------------- SHINY APP BEGIN --------------------- ##
  shinyApp(

## Define UI and Server ##
ui <- CRISPRclassify:::getShinyUI(css),
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

    ## Check file extension ##
    print(paste0('Input file: ', input$inFile$datapath))
    fileExt <- CRISPRclassify:::getFileExt(input$inFile$datapath)
    if (!fileExt %in% c('fasta','fna','fa')){
      CRISPRclassify:::presentModal('File Type Error', 'Invalid file format provided. Please provide a properly formatted .fasta file with a .fasta, .fna, or .fa extension.')
      return(CRISPRclassify:::createDummyTable())
    }

    on.exit(values$fileStatus <- 'empty')
    progress <- shiny::Progress$new()
    on.exit(progress$close(), add=TRUE)

    # Extract repeats from .crisprs file
    print(input$inFile$datapath)
    progress$set(message = "Finding CRISPR repeats.", detail = 'This may take a few moments...', value = .05)
    system(sprintf("%s/lib/crVSkeleton.sh -c -f %s",wrkDir, input$inFile$datapath), intern = TRUE) #TRUE

    print('Checking CRISPRs exist: ')
    if (CRISPRclassify:::checkCrisprExists(input$inFile$datapath)){
      print('A valid CRISPRS file was generated.')
    } else {
      print('No valid CRISPRS file was found.')
      CRISPRclassify:::presentModal('No CRISPR loci detected', 'No detectable CRISPR loci could be found in this genome file. Upload a different genome/metagenome to continue searching for CRISPR loci.')
      return(CRISPRclassify:::createDummyTable())
    }

    progress$set(message = "Tokenizing repeats.", detail = 'This may take a few moments...', value = .1)
    analytic_df <- CRISPRclassify:::getSequenceAnalyticDf(input$inFile$datapath)

    if (nrow(analytic_df) == 1){
      firstRow <- analytic_df %>% slice_head(n=1)

      if (firstRow$seq == ""){
        # all potential CRISPR loci were filtered out during false positive filtering
        CRISPRclassify:::presentModal('No valid CRISPR loci found', 'No valid CRISPR loci could be found after filtering for false positives. Upload a different genome/metagenome to continue searching for CRISPR loci.')
        return(analytic_df)
      }

    }

    progress$set(message = "Predicting subtypes.", detail = "This may take a few moments...",  value = .1)
    rename = FALSE
    CRISPRclassify:::getXGBoostPredictions(analytic_df, progress, xg_matrix_format_file, cas_vect, wrkDir, rename) %>%
      return()

  })

  output$classified_results <- DT::renderDataTable({
    parseResults() %>%
      rowwise() %>%
      mutate(
             seq = CRISPRclassify:::getSyntaxHighlightedSequence(seq, cas_type),
             probability = round(probability, 3)) %>%
      ungroup() %>%
      mutate(cas_type = as.factor(cas_type)) %>%
      rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Contig` = contig_name, `Range` = range, `Locus` = locus, `Repeat Count` = distinct_repeat_count, `Closest Strain` = closestStrain.org_desc, `Edit Dist` = edit_dist) %>%
      arrange(Locus) %>%
      dplyr::select(Locus, `Contig`,  Range, `Subtype`, `Repeat`, `Repeat Count`, `Probability`, `Closest Strain`, `Edit Dist`)
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
        rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Contig Name` = contig_name, `Range` = range, `Locus` = locus, `Repeat Count` = distinct_repeat_count, `Closest Strain` = closestStrain.org_desc, `Closest Strain Accession` = closestStrain.accession, `Closest Strain Location` = closestStrain.location, `Closest Strain Subtype` = closestStrain.type, `Closest Strain Repeat` = closestStrain.repeat_seq, `Edit Dist` = edit_dist) %>%
        arrange(Locus) %>%
        dplyr::select(Locus, `Contig Name`,  Range, `Subtype`, `Repeat`, `Repeat Count`, `Probability`, `Closest Strain`, `Closest Strain Accession`, `Closest Strain Location`, `Closest Strain Subtype`, `Closest Strain Repeat`, `Edit Dist`) %>%
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

#' @examples
#' classifyRepeats("/path/to/repeat.txt")
#' @export
classifyRepeats <- function(repeatFile){

  ## Java check ##
  CRISPRclassify:::validateJava()

  ## Check if file exists ##
  CRISPRclassify:::checkFileExists(repeatFile)

  ## Check file extension ##
  validateFileArg(repeatFile, c('csv','txt','tsv'))

  ## Load HC vars ##
  wrkDir <- CRISPRclassify:::getWrkDir()
  xg_matrix_format_file <- CRISPRclassify:::getMatrixFormatFile(wrkDir)
  cas_vect <- CRISPRclassify:::getCasVect()

  ## Read in .txt file ##
  repeat_df <- read_csv(repeatFile, col_names=FALSE, col_types = cols()) %>%
    rename(sequence = X1)

  ## Derive bio features ##
  repeat_expanded_df <- deriveBioFeatures(repeat_df)

  ## Classify ##
  writedf <- CRISPRclassify:::getXGBoostPredictionsForRPPL(repeat_expanded_df, xg_matrix_format_file, cas_vect, wrkDir)
  writedf %>% write.csv(paste0(repeatFile,'.crclass'), row.names = FALSE)

  return(writedf)
}

#' @examples
#' classifyFasta('/path/to/file.fasta')
#' @export
classifyFasta <- function(fastaFile){

  ## Java check ##
  CRISPRclassify:::validateJava()

  ## Check if file exists ##
  CRISPRclassify:::checkFileExists(fastaFile)

  ## Check file extension ##
  validateFileArg(fastaFile, c('fasta','fna','fa'))

  ## Load HC vars ##
  wrkDir <- CRISPRclassify:::getWrkDir()
  xg_matrix_format_file <- CRISPRclassify:::getMatrixFormatFile(wrkDir)
  cas_vect <- CRISPRclassify:::getCasVect()

  ## Read in .fasta file ##
  pathToFasta <- CRISPRclassify:::processFastaFile(fastaFile)

  print('Checking CRISPR exists: ')
  if (CRISPRclassify:::checkCrisprExists(pathToFasta)){
    print('A valid CRISPRS file was generated.')
  } else {
    stop('ERROR: No detectable CRISPR loci could be found in this genome file. Upload a different genome/metagenome to continue searching for CRISPR loci.')
  }

  analytic_df <- CRISPRclassify:::getSequenceAnalyticDf(pathToFasta)

  if (nrow(analytic_df) == 1){
    firstRow <- analytic_df %>% slice_head(n=1)
    if (firstRow$seq == ""){
      # all potential crispr loci were filtered out during false positive filtering
      stop("ERROR: No valid CRISPR loci could be found after filtering for false positives. Upload a different genome/metagenome to continue searching for CRISPR loci.")
    }
  }

  rename = TRUE
  writedf <- CRISPRclassify:::getXGBoostPredictions(analytic_df, NULL, xg_matrix_format_file, cas_vect, wrkDir, rename)
  writedf %>% write.csv(paste0(fastaFile,'.crclass'), row.names = FALSE)

  return(writedf)
}




