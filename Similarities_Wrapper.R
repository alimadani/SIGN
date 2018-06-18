Similarities_Wrapper <- function(ExpMat_Test, ExpMat_Ref, GeneVec, PathwaySet, RefID, TestClassIter, SampleIter){
 
  rownames(ExpMat_Ref) <- NULL

  if(TestClassIter == RefID){
   ExpMat_Ref <- ExpMat_Ref[,-SampleIter]
  }
  SimMat <- Pathway_similarity(matrix(ExpMat_Test[,SampleIter], ncol = 1), ExpMat_Ref, GeneVec, PathwaySet, RefID)

  GSVA_Out <- GSVA_Calculation(matrix(ExpMat_Test[,SampleIter], ncol = 1), ExpMat_Ref, GeneVec, PathwaySet, RefID)

  Similarities <- list(SimMat, GSVA_Out)
  names(Similarities) <- c("Pathway_SimMat", "GSVA_Out")
  return(Similarities)
}
