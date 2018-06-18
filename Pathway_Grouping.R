Pathway_Grouping <- function(PathwayDir, Pattern){
 
 PathwayFiles <- list.files(PathwayDir, Pattern)

 if(grepl("rds", Pattern)){
  PathwayList <- list()
  for(FileIter in 1:length(PathwayFiles)){
   PathwayList[[FileIter]] <- readRDS(paste(PathwayDir, PathwayFiles[FileIter],sep = "", collapse = ""))
  }
  names(PathwayList) <- gsub(Pattern, "", PathwayFiles)
 }else{
  stop("Please use Pathways_Preprocess_MSigDB.R script and save the files as rds!")
 }

 return(PathwayList)
}
