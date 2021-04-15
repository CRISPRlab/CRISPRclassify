
getWrkDir <- function(){
  return (fs::path_package("CRISPRclassify"))
}

getMatrixFormatFile <- function(wrkDir){
  return(str_c(wrkDir,"/csvRef/matrix_format_2020-12-04_nchar_5.csv"))
}

getFeatImpFile <- function(wrkDir){
  return (str_c(wrkDir,"/csvRef/xg_feature_importance_5_2020-12-04.csv"))
}

getCasVect <- function(){
  return(c('I-A','I-B','I-C','I-D','I-E','I-F','I-U','II-A','II-B','II-C','III-A','III-B','III-C','III-D','III-E','III-F','IV-A','IV-C','V-A','V-B1','V-B2','V-F','V-J','V-U1','V-U2','V-U4','VI-A','VI-B','VI-C','VI-D'))
}

getCharCount <- function(){
  return(5)
}
