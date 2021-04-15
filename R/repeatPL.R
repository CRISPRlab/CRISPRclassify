## These functions are specific to the repeat-based pipeline ##
deriveBioFeatures <- function(repeat_df){
  repeat_gc_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getGC)
  repeat_palidx_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getPalIdx)
  repeat_length_vect <- sapply(repeat_df$sequence, CRISPRclassify:::getLength)

  repeat_expanded_df <- repeat_df %>% mutate(
    locus = 1,
    id = row_number(),
    gc = repeat_gc_vect,
    palIdx = repeat_palidx_vect,
    length = repeat_length_vect)

  return(repeat_expanded_df)
}

getXGBoostPredictionsForRPPL <- function(sequence_analytic_df, xg_matrix_format_file, cas_vect, wrkDir){

  print("...Getting xgboost predictions...")

  getXGModelCache <- CRISPRclassify:::getXGModel

  all_feat_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% dplyr::select(X3) %>% distinct() %>% pull()

  # add seq and unique id
  sequence_analytic_df <- sequence_analytic_df %>%
    rename(seq = sequence) %>%
    mutate(unique_id = row_number(), contig_name = 'manual', range = 'unknown', distinct_repeat_count = 1)
  token_df <- CRISPRclassify:::getTokenDF(sequence_analytic_df)
  bio_feat_df <- sequence_analytic_df %>% inner_join(token_df, by = "unique_id")

  missing_col_vect <- setdiff(all_feat_vect , colnames(bio_feat_df))
  missing_feature_df <- matrix( rep( 0, len=length(missing_col_vect)*nrow(bio_feat_df)), nrow = nrow(bio_feat_df)) %>% as.data.frame()
  names(missing_feature_df) <- missing_col_vect

  bio_feat_df <- bio_feat_df %>% bind_cols(missing_feature_df)

  result_df <- tibble(seq = as.character(), cas_type = as.character(), probability = as.double())

  # cas <- 'I-E'
  allCasN <- length(cas_vect)
  adjustmentInc <- (allCasN + (.1 * allCasN)) ## since we start on .1, account for that ##
  for (cas in cas_vect){

    print(str_c("running: ", cas))

    if (file.exists(CRISPRclassify:::getModelFilePath(cas, wrkDir))){
      cas_feature_vect <- read_csv(xg_matrix_format_file, col_names = F) %>% filter(X1 == cas) %>% dplyr::select(X3) %>% pull()
      # feature_df <- bio_feat_df %>% select_at(cas_feature_vect)
      feature_df <- bio_feat_df %>% select(one_of(cas_feature_vect))
      feature_mat <- feature_df %>%  as.matrix()
      bstSparse <- getXGModelCache(cas, wrkDir)
      feature_sparse_mat <- as(feature_mat, "dgCMatrix")

      if(nrow(feature_sparse_mat) > 0){

        pred <- predict(bstSparse, feature_sparse_mat)
        temp_result_df <- sequence_analytic_df %>% bind_cols(probability = pred) %>% dplyr::select(contig_name, locus, id, seq, probability, range, distinct_repeat_count) %>% distinct() %>% mutate(cas_type = cas)

        result_df <- result_df %>% bind_rows(temp_result_df)

      } else {
        print(str_c("predictions not possible for subtype ", cas))
      }

    }

    else {print("Subtype excluded due to low volume")}
  }

  max_result_df <- result_df %>% group_by(id) %>% filter(probability == max(probability)) %>% ungroup() %>%  select(seq,cas_type , probability)

  lkpMaster <- read_csv(str_c(wrkDir, '/csvRef/master_repeat_lkp.csv'))

  # Lookup closest matching strain by edit distance #
  max_result_df_lkp <- transform(max_result_df, closestStrain = CRISPRclassify:::strainLkp(seq, lkpMaster)) %>%
    mutate(edit_dist = CRISPRclassify:::minEditDist(seq, closestStrain.repeat_seq))

  print(max_result_df_lkp)

  return(max_result_df_lkp %>%
           mutate(cas_type = case_when(cas_type == 'V-J' ~ 'V-K', cas_type == "I-U" ~ "I-G", T ~ cas_type)) %>%
           rename( `Repeat` = seq, `Subtype` = cas_type, `Probability` = probability, `Closest Strain` = closestStrain.org_desc, `Closest Strain Accession` = closestStrain.accession, `Closest Strain Location` = closestStrain.location, `Closest Strain Subtype` = closestStrain.type, `Closest Strain Repeat` = closestStrain.repeat_seq, `Edit Dist` = edit_dist)
         )
}
