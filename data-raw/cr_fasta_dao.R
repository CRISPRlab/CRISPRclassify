library(stringi)
#library(tidyverse)
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

