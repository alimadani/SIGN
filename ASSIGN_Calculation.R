ASSIGN_Calculation  <- function(ExpMat_Test, ExpMat_Ref, GeneVec, PathwaySet){

 ASSIGN_Vec <- c()
 for(PathwayIter in 1:length(PathwaySet)){
   TargetGenes  <- PathwaySet[[PathwayIter]]
   MatchedInd <- which(GeneVec %in% TargetGenes)

  if(length(MatchedInd) > 1){
    RandomNumber <- runif(1, 100,10000)
    system(paste("mkdir ", MainDir, "TMPDir_", RandomNumber, sep = "", collapse = ""))
    invisible(assign.wrapper(trainingData = data.frame(ExpMat_Ref), testData = data.frame(ExpMat_Test),
               trainingLabel = list(control = list(ctrl=1:ncol(ExpMat_Ref)), ctrl = 1:ncol(ExpMat_Ref)),
               testLabel=rep(c("test"), c(1:ncol(ExpMat_Test))),
               geneList = list(sig1 = rownames(ExpMat1)[MatchedInd]),
               outputDir=paste(MainDir, "TMPDir_", RandomNumber, sep = "", collapse = "")))

    Training_ASSIGN <- read.csv(paste(MainDir, "TMPDir_", RandomNumber, "/pathway_activity_trainingset.csv",
                                  sep = "", collapse = ""), stringsAsFactors = F, check.names = F)
    Training_ASSIGN <- as.numeric(Training_ASSIGN[,2])
    Test_ASSIGN <- read.csv(paste(MainDir, "TMPDir_", RandomNumber, "/pathway_activity_testset.csv",
                              sep = "", collapse = ""), stringsAsFactors = F, check.names = F)
    ASSIGN_Vec <- c(ASSIGN_Vec, median(as.numeric(Test_ASSIGN[,2])))
    system(paste("rm -r ", MainDir, "TMPDir_", RandomNumber, sep = "", collapse = ""))
   }
 }

 return(ASSIGN_Vec)

}

