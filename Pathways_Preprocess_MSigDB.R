
MainDir <- "~/Desktop/SIGN/Pathways_Dataset/"
TargetFiles <- list.files(MainDir, pattern = "gmt")
DSnames <- unlist(lapply(TargetFiles, function(X){paste(strsplit(X, "\\.")[[1]][
  c(1:(length(strsplit(X, "\\.")[[1]])-4), (length(strsplit(X, "\\.")[[1]])-1))], sep = "_", collapse = "_")}))

Pathway_Preprocess <- function(MainDir, File, Name){
  TargetDB_Pathways <- read.table(paste(MainDir, File, sep = "", collapse = ""), fill = T, sep = "\t")
  TargetDB_Pathways <- as.matrix(TargetDB_Pathways)
  
  TargetDB_Names <- as.character(TargetDB_Pathways[,1])
  
  TargetDB_Pathways <- as.matrix(TargetDB_Pathways)
  if(grepl("entrez", File)){
    TargetDB_Genes <- apply(TargetDB_Pathways, 1, function(X){as.numeric(X[3:length(X)])[which(!is.na(as.numeric(X[3:length(X)])))]})
  }else if(grepl("symbols", File)){
    TargetDB_Genes <- apply(TargetDB_Pathways, 1, function(X){as.character(X[3:length(X)])[which(!is.na(as.character(X[3:length(X)])) & 
                                                                                                   as.character(X[3:length(X)]) != "")]})
  }

  
  names(TargetDB_Genes) <- TargetDB_Names
  
  return(TargetDB_Genes)
}


invisible(lapply(seq(1:length(TargetFiles)), function(FileIter){saveRDS(Pathway_Preprocess(MainDir, TargetFiles[FileIter], DSnames[FileIter]),
                                                              paste(MainDir, DSnames[FileIter], "_ProcessPathways.rds", sep = "", collapse = ))}))

################
GO_ParentRemoval <- function(GOlist){

  ParentInd <- c()
  for(GOiter in 1:length(GOlist)){ #
    print(paste("GO", GOiter, sep = "_"))
    TargetPathGenes <- GOlist[[GOiter]]
    if(length(TargetPathGenes) >= 5){
      ParentInd_Tmp <- lapply(GOlist, function(X){(length(intersect(TargetPathGenes, X)) == length(TargetPathGenes)
                                                   & length(X) > length(TargetPathGenes))})
      ParentInd <- c(ParentInd, unique(which(as.character(unlist(ParentInd_Tmp)) == "TRUE")))
    }else{
      ParentInd <- c(ParentInd, GOiter)
    }
    
    print(length(unique(ParentInd)))
  }
  ParentInd <- unique(ParentInd)
  print(length(ParentInd))
  if(length(ParentInd) > 0){
    GOlist <- GOlist[-ParentInd]
  }
  return(GOlist)
}
#################
GO_Genes <- readRDS(paste(MainDir, "c5_all_entrez_ProcessPathways.rds", sep = "", collapse = ""))
GO_Genes <- GO_ParentRemoval(GO_Genes)

saveRDS(GO_Genes, file = paste(MainDir, "c5_all_entrez_ProcessPathways.rds", sep = "", collapse = ""))
##################
##################
GO_Genes <- readRDS(paste(MainDir, "c5_all_symbols_ProcessPathways.rds", sep = "", collapse = ""))
GO_Genes <- GO_ParentRemoval(GO_Genes)

saveRDS(GO_Genes, file = paste(MainDir, "c5_all_symbols_ProcessPathways.rds", sep = "", collapse = ""))
