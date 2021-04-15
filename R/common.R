checkFileExists <- function(fileWithPath){
  if (!file.exists(fileWithPath)){
    stop(paste0('ERROR: the specified file or path could not be located on your system: "',fileWithPath,'"'))
  }
}

getFileExt <- function(fileWithPath){
  ext <- str_split(fileWithPath, '[.]')
  if (length(ext[[1]]) < 2){ # file doesnt have period to parse
    return("invalid")
  } else {
    return(ext[[1]][length(ext[[1]])])
  }
}

validateFileArg <- function(fileWithPath, expected){
  fileExt <- getFileExt(fileWithPath)
  if (fileExt == "invalid" || !fileExt %in% expected){
    stop(paste0('ERROR: Invalid file format provided. You passed in a "',fileExt,'" file, but a "',paste(expected,collapse=" or "), '" was expected.'))
  }
}

getCompliment <- function(nuc){
  c("T","A", "C", "G")[which(str_starts(nuc, c("A","T", "G", "C")))]
}

getReverseCompliment <- function(sequence){
  rc_vect <- as.vector(character())
  for (nuc in str_split(stri_reverse(toupper(sequence)) ,'')[[1]]) {
    rc_vect <- c(rc_vect, getCompliment(nuc))
  }
  return(str_c(rc_vect, collapse = ''))
}

getModelFilePath <- function(cas, wrkDir) {
  return(str_c(wrkDir,"/models/",cas,"_2020-12-04_nchar_5"))
}

getXGModel <- function(cas, wrkDir){
  return(xgb.load(getModelFilePath(cas, wrkDir)))
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

getContigName <- function(id){
  return(str_split(id, "_CRISPR")[[1]][1])
}

checkCrisprExists <- function(datapath){
  crispr_file = paste(datapath, ".crisprs", sep="")
  crispr_df <- read_csv(crispr_file, col_names = F)

  if (nrow(crispr_df) <= 0){
    return(FALSE)
  }

  return(TRUE)
}

get_tokens <- function(seq, id) {
  n_char_cnt <- CRISPRclassify:::getCharCount()
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
  idx <- amatch(repeatSeq, masterLkpDf$repeat_seq, nomatch=(nrow(masterLkpDf)), maxDist = 20)
  org <- masterLkpDf[idx,]
  return(org)
}

minEditDist <- function(seq1, seq2){
  seq2 <- na_if(seq2, "Unknown")
  fwdDist <- mapply(adist, seq1, seq2)
  revDist <- mapply(adist, seq1, getReverseCompliment(seq2))
  minList <- mapply(min, fwdDist, revDist)
  return(minList)
}

createDummyTable <- function(){
  dummy_df <- data.frame(seq=character(),
                         cas_type=character(),
                         distinct_repeat_count=integer(),
                         probability=double(),
                         contig_name=character(),
                         range=character(),
                         locus=integer(),
                         closestStrain.org_desc=character(),
                         closestStrain.accession=character(),
                         closestStrain.location=character(),
                         closestStrain.type=character(),
                         closestStrain.repeat_seq=character(),
                         edit_dist=integer(),
                         stringsAsFactors=FALSE)

  dummy_df <- dummy_df %>% add_row(seq='',cas_type='',distinct_repeat_count=0,probability=0,contig_name='',range='',locus=0,closestStrain.org_desc='',closestStrain.accession='',closestStrain.location='',closestStrain.type='',closestStrain.repeat_seq='',edit_dist=0)
  return(dummy_df)
}

getRepeatDf <- function(repeat_file, datapath){
  print('...getting repeat df...')
  repeat_df <- getFasta(repeat_file)
  repeat_locus_vect <- sapply(repeat_df$id, getRepeatLocus)
  repeat_id_vect <- sapply(repeat_df$id, getRepeatId)

  repeat_gc_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getGC)
  repeat_palidx_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getPalIdx)
  repeat_length_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getLength)

  repeat_expanded_df <- repeat_df %>% mutate(
    locus = repeat_locus_vect,
    id = repeat_id_vect,
    gc = repeat_gc_vect,
    palIdx = repeat_palidx_vect,
    length = repeat_length_vect)

  repeat_expanded_df <- CRISPRclassify:::filterFalsePositives(repeat_expanded_df, .7, "repeat")

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
  print('...getting spacer df...')
  spacer_df <- getFasta(spacer_file)
  contig_name_vect =sapply(spacer_df$id, getContigName)
  spacer_locus_vect <- sapply(spacer_df$id, getSpacerLocus)
  spacer_id_vect <- sapply(spacer_df$id, getSpacerId)
  spacer_length_vect <- sapply(spacer_df$sequence, getLength)

  spacer_df <- spacer_df %>% mutate(locus = spacer_locus_vect,
                                    id = spacer_id_vect,
                                    length = spacer_length_vect,
                                    contig_name = contig_name_vect)

  spacer_df <- CRISPRclassify:::filterFalsePositives(spacer_df, .6, "spacer") #.6
  spacer_df <- spacer_df %>% select(-id) %>%
    distinct()

  all_fasta_dtl_df <- getAllFastaDtl(datapath) %>% select(spacer_sequence, contig_name, locus_id,range) %>% distinct()

  spacer_df <- spacer_df %>% distinct() %>%  inner_join(all_fasta_dtl_df, by = c("sequence" = "spacer_sequence" , "contig_name", "locus" = "locus_id"))

  spacer_df <- spacer_df %>% group_by(locus, contig_name, range) %>%
    summarise(spacer_n = n(), avg_spacer_length = mean(length)) %>%
    ungroup()
  return(spacer_df)
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

percentIden <- function(seq1, seq2){
  dist <- stringdist(seq1, seq2, method = "lv")
  similarity <- 1 - dist / nchar(as.character(seq1))
  return(similarity)
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
    rowwise() %>%
    mutate(contig_id = trimws(str_replace_all(contig_id," bp\\)", "")),
           locus_id = trimws(str_sub(locus, 7, nchar(strsplit(locus, 'Range')[[1]][1])))) %>%
    separate(locus, c(NA, "range"), sep = "Range: ") %>%
    mutate(range = str_replace_all(range, "\\s" , ""))  #%>% group_by(contig_name, contig_id, range, repeat_sequence, locus_id) %>% tally() %>% rename(spacer_n = n) %>% ungroup()

  return(all_info_fasta_df)

}

getSequenceAnalyticDf <- function(datapath){

  spacer_file = paste(datapath, "_spacers.fa", sep="")
  repeat_file = paste(datapath, "_repeats.fa", sep="")

  print('...getting all fasta df...')

  all_fasta_df <- getAllFastaDtl(datapath)

  #get spacer count and average spacer length
  #getting repeat features
  repeat_df <- getRepeatDf(repeat_file, datapath)
  spacer_df <- getSpacerDf(spacer_file, datapath)

  all_fasta_df %>% left_join(repeat_df)

  if (nrow(repeat_df) <= 0 || nrow(spacer_df) <= 0){
    return (CRISPRclassify:::createDummyTable())
  } else {
    sequence_analytic_df <- repeat_df %>% inner_join(spacer_df, by = c("locus", "contig_name")) %>%
      mutate(repeat_spacer_ratio = average_repeat_length/avg_spacer_length,
             unique_id = str_c(locus, "_", id, "_", contig_name)) %>%
      rename(seq = sequence )

    return(sequence_analytic_df)
  }
}

getXGBoostPredictions <- function(sequence_analytic_df, progress, xg_matrix_format_file, cas_vect, wrkDir, rename){

  print("...Getting xgboost predictions...")

  getXGModelCache <- getXGModel

  all_feat_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% dplyr::select(X3) %>% distinct() %>% pull()
  token_df <- CRISPRclassify:::getTokenDF(sequence_analytic_df)
  bio_feat_df <- sequence_analytic_df %>% inner_join(token_df, by = "unique_id")

  missing_col_vect <- setdiff(all_feat_vect , colnames(bio_feat_df))
  missing_feature_df <- matrix( rep( 0, len=length(missing_col_vect)*nrow(bio_feat_df)), nrow = nrow(bio_feat_df)) %>% as.data.frame()
  names(missing_feature_df) <- missing_col_vect

  bio_feat_df <- bio_feat_df %>% bind_cols(missing_feature_df)

  result_df <- tibble(unique_id = as.character(), seq = as.character(), cas_type = as.character(), probability = as.double())

  allCasN <- length(cas_vect)
  adjustmentInc <- (allCasN + (.1 * allCasN)) ## since we start on .1, account for that ##
  for (cas in cas_vect){

    print(str_c("running: ", cas))
    if (!is.null(progress)){
      progress$inc(1/adjustmentInc, message = str_c("Checking ", cas, "."))
    }

    if (file.exists(CRISPRclassify:::getModelFilePath(cas, wrkDir))){
      cas_feature_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% filter(X1 == cas) %>% dplyr::select(X3) %>% pull()
      # feature_df <- bio_feat_df %>% select_at(cas_feature_vect)
      feature_df <- bio_feat_df %>% select(one_of(cas_feature_vect))
      feature_mat <- feature_df %>%  as.matrix()
      bstSparse <- getXGModelCache(cas, wrkDir)
      feature_sparse_mat <- as(feature_mat, "dgCMatrix")

      if(nrow(feature_sparse_mat) > 0){

        pred <- predict(bstSparse, feature_sparse_mat)
        temp_result_df <- sequence_analytic_df %>% bind_cols(probability = pred) %>% dplyr::select(contig_name, unique_id, locus, id, seq, probability, range, distinct_repeat_count) %>% distinct() %>% mutate(cas_type = cas)

        result_df <- result_df %>% bind_rows(temp_result_df)
      } else {
        print(str_c("predictions not possible for subtype ", cas))
      }

    }
    else {print("Subtype excluded due to low volume")}
  }

  if (!is.null(progress)){
    progress$set(message = "Looking up closest strains. ", detail = 'Just a moment...', value = .98)
  }
  max_result_df <- result_df %>% group_by(unique_id) %>% filter(probability == max(probability)) %>% ungroup() %>%  select(seq,cas_type, distinct_repeat_count , probability,contig_name, range, locus)

  lkpMaster <- read_csv(str_c(wrkDir, '/csvRef/master_repeat_lkp.csv'))

  # Lookup closest matching strain by edit distance #
  max_result_df_lkp <- transform(max_result_df, closestStrain = CRISPRclassify:::strainLkp(seq, lkpMaster)) %>%
    mutate(edit_dist = CRISPRclassify:::minEditDist(seq, closestStrain.repeat_seq))

  if (!is.null(progress)){
    progress$set(message = "Wrapping up. ", detail = 'Just a moment...', value = 1)
  }
  print(max_result_df_lkp)

  if (rename){
    return(max_result_df_lkp %>%
             mutate(cas_type = case_when(cas_type == 'V-J' ~ 'V-K', cas_type == "I-U" ~ "I-G", T ~ cas_type)) %>%
             rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Contig Name` = contig_name, `Range` = range, `Locus` = locus, `Repeat Count` = distinct_repeat_count, `Closest Strain` = closestStrain.org_desc, `Closest Strain Accession` = closestStrain.accession, `Closest Strain Location` = closestStrain.location, `Closest Strain Subtype` = closestStrain.type, `Closest Strain Repeat` = closestStrain.repeat_seq, `Edit Dist` = edit_dist)
           )
  } else {
    return(max_result_df_lkp %>%
             mutate(cas_type = case_when(cas_type == 'V-J' ~ 'V-K', cas_type == "I-U" ~ "I-G", T ~ cas_type)))
  }
}
