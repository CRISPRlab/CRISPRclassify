processFastaFile <- function(fastaFile){
  ## Move fasta to tmp directory ##
  datapath <- tempdir()
  print(paste0('...copying ', fastaFile,' to ', datapath))
  file.copy(fastaFile, datapath)

  ## Process file & extract repeats/spacers ##
  fileArr <- str_split(fastaFile, '/')
  fileName <- fileArr[[1]][length(fileArr[[1]])]
  file.rename(paste0(datapath,'/',fileName), paste0(datapath,'/0.fasta'))
  fastaWithPath <- paste0(datapath,'/0.fasta')
  system(sprintf("%s/lib/crVSkeleton.sh -c -f %s",CRISPRclassify:::getWrkDir(), fastaWithPath, intern = TRUE))
  return(fastaWithPath)
}
