
xg_matrix_format_file <- "matrix_format_2020-12-04_nchar_5.csv"
feat_imp_file <- "xg_feature_importance_5_2020-12-04.csv"
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
    # print(".....................")
    # print(formatted_seq)
    formatted_seq <- str_replace_all(formatted_seq,getReverseCompliment(feat), ngramToHtmlBlue)
    formatted_seq <- str_replace_all(formatted_seq,feat, ngramToHtmlOrange)
    
    # print(getReverseCompliment(feat))
    # print(formatted_seq)
    # print(formatted_seq)
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
  return(str_c("models/",cas,"_2020-12-04_nchar_5"))
}
getXGModel <- function(cas){
  return(xgb.load(getModelFilePath(cas)))
}
getXGModelCache <- getXGModel

#getRepeatDF <- function(){
#  wrkDir <- getwd()
#  
#   read_csv(sprintf('%s/data/meta.fasta_repeats.fa', wrkDir), col_names = F) %>%
#    mutate(X2 = lead(X1)) %>%
#    dplyr::rename(id = X1, sequence = X2) %>%
#    filter(!str_detect(sequence, ">")) %>%
#    mutate(id = str_replace_all(id, ">", "")) %>% return()
#}

getTokenDF <- function(sequence_analytic_df){
 # c('test_id')
 # c('GTCGCTCCCCGTGCGGGAGCGTGGATAGCCGA')
 #   token_df <- do.call("rbind", mapply(get_tokens,  c('GTCGCTCCCCGTGCGGGAGCGTGGATAGCCGA'), c('test_id'), SIMPLIFY = F)) %>% 
 #     rename(value = token_count) %>% 
 #     mutate(value = as.numeric(value)) %>% 
 #     dplyr::select(id, feature, value) %>% 
 #     rename(unique_id = id) %>% 
 #     pivot_wider(names_from = feature, values_from = value, values_fill = 0)
 #   
 #   token_df %>% write_csv("cr_local_test_record.csv")
    
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

getXGBoostPredictions <- function(sequence_analytic_df){
  print("...Getting xgboost predictions...")
  all_feat_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% dplyr::select(X3) %>% distinct() %>% pull()
  token_df <- getTokenDF(sequence_analytic_df) 
  bio_feat_df <- sequence_analytic_df %>% inner_join(token_df, by = "unique_id")
  
  missing_col_vect <- setdiff(all_feat_vect , colnames(bio_feat_df))
  missing_feature_df <- matrix( rep( 0, len=length(missing_col_vect)*nrow(bio_feat_df)), nrow = nrow(bio_feat_df)) %>% as.data.frame() 
  names(missing_feature_df) <- missing_col_vect

  bio_feat_df <- bio_feat_df %>% bind_cols(missing_feature_df)
  bio_feat_df %>% write_csv("seq_analytic_df_audit.csv")
  
  result_df <- tibble(unique_id = as.character(), seq = as.character(), cas_type = as.character(), probability = as.double())
  
   # cas <- 'I-E'
for (cas in cas_vect){
  
  print(str_c("running: ", cas))
  
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
}
  
#result_df %>% write_csv("_tmp_full_output.csv")
max_result_df <- result_df %>% group_by(unique_id) %>% filter(probability == max(probability)) %>% ungroup() %>%  select(seq,cas_type, distinct_repeat_count , probability,contig_name, range, locus)
lkpMaster <- read_csv('master_repeat_lkp.csv')
max_result_df_lkp <- transform(max_result_df, closestStrain = strainLkp(seq, lkpMaster))
return(max_result_df_lkp)
}