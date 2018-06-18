Pathway_similarity <- function(ExpMat1, ExpMat2, GeneVec, PathwaySet, Name){
  
  TSC_Vec <- c()
  BubbleSort_Vec <- c()
  Pearson_Vec <- c()
  for(PathwayIter in 1:length(PathwaySet)){
   TargetGenes  <- PathwaySet[[PathwayIter]]
   MatchedInd <- which(GeneVec %in% TargetGenes)

   if(length(MatchedInd) > 1){
    PathExpMat1 <- matrix(ExpMat1[MatchedInd,], ncol = ncol(ExpMat1))
    PathExpMat2 <- matrix(ExpMat2[MatchedInd,], ncol = ncol(ExpMat2))
    TSC_Vec <- c(TSC_Vec, TSC(PathExpMat1, PathExpMat2))

    MedianVecTMP1 <- as.numeric(apply(PathExpMat1,1 ,function(X){median(na.omit(as.numeric(X)))}))
    MedianVecTMP2 <- as.numeric(apply(PathExpMat2,1 ,function(X){median(na.omit(as.numeric(X)))}))

    BubbleSort_Vec <- c(BubbleSort_Vec, as.numeric(BubbleSort(MedianVecTMP1, MedianVecTMP2)))
    Pearson_Vec <- c(Pearson_Vec, as.numeric(cor.test(MedianVecTMP1, MedianVecTMP2)$estimate))
   }
 }

  Similarity_Output <- cbind(TSC_Vec, BubbleSort_Vec, Pearson_Vec)
  colnames(Similarity_Output) <- paste(Name, c("TSC_Vec", "BubbleSort_Vec", "Pearson_Vec"), sep = "_")

  return(Similarity_Output)
}
