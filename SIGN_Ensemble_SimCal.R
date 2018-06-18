SIGN_Ensemble_SimCal <- function(ExpList, RefClassID, TestClassID, GeneID, PathwaySets){

  ExpList <- GeneMatching(ExpList)

  if(GeneID == "Symbol"){
    GeneVec <- as.character(rownames(ExpList[[1]]))
  }else if(GeneID == "EntrezID"){
    GeneVec <- as.numeric(rownames(ExpList[[1]]))
  }

   RefID <- 1
   ExpMat_Ref1 <- ExpList[[RefClassID[RefID]]]

   RefID <- 2
   ExpMat_Ref2 <- ExpList[[RefClassID[RefID]]]

   PathSet_PathSim_List <- list()
   PathSet_GeneSim_List <- list()
  for(PathSetIter in c(1:length(PathwaySets))){ 
   print(names(PathwaySets)[PathSetIter]) 
   PathwaySet <- PathwaySets[[PathSetIter]]
  
   PathSim_List <- list()
   GeneSim_List <- list()

   for(TestClassIter in TestClassID){
     print(TestClassIter)
     ExpMat_Test <- ExpList[[TestClassIter]]
     rownames(ExpMat_Test) <- NULL
     colnames(ExpMat_Test) <- NULL

    PathSim_Mat <- c()
    GeneSim_Mat <- c()
   
    for(SampleIter in 1:ncol(ExpMat_Test)){
      print(SampleIter)
      SimOut_Ref1 <- Similarities_Wrapper(ExpMat_Test, ExpMat_Ref1, GeneVec, PathwaySet, RefID=RefClassID[1], TestClassIter,  SampleIter)
     
      SimOut_Ref2 <- Similarities_Wrapper(ExpMat_Test, ExpMat_Ref2, GeneVec, PathwaySet, RefID=RefClassID[2], TestClassIter, SampleIter)
      
      PathSim_Mat <- rbind(PathSim_Mat, c(SimOut_Ref1[["GSVA_Out"]], SimOut_Ref2[["GSVA_Out"]],
                            SimSummary_2Class(SimOut_Ref1[["Pathway_SimMat"]], SimOut_Ref2[["Pathway_SimMat"]])))

      #ASSIGN_Mat <- rbind(ASSIGN_Mat, ASSIGN_Wrapper(ExpMat_Test, ExpMat_Ref1, ExpMat_Ref2, GeneVec, PathwaySet))
      if(PathSetIter == 1){
       GeneSim_Mat <- rbind(GeneSim_Mat, Genes_SimCal(ExpMat_Test, ExpMat_Ref1, ExpMat_Ref2, RefClassID, TestClassIter, SampleIter))
      }
    }
    PathSim_List[[TestClassIter]] <- PathSim_Mat
    if(PathSetIter == 1){
     GeneSim_List[[TestClassIter]] <- GeneSim_Mat
    }
   }
   PathSet_PathSim_List[[PathSetIter]] <- PathSim_List
   if(PathSetIter == 1){
    PathSet_GeneSim_List[[PathSetIter]] <- GeneSim_List
   }
  }

  names(PathSet_PathSim_List) <- names(PathwaySets)
  names(PathSet_GeneSim_List) <- names(PathwaySets)[1]

  PathSet_SimList <- list(PathSet_PathSim_List, PathSet_GeneSim_List)
  names(PathSet_SimList) <- c("PathSim", "GeneSim")

  return(PathSet_SimList)
}

