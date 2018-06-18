
options(stringsAsFactors = F)
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 4){
  cat("Uasge: Rscript SIGN_Classification_L1out.R <Path> <Query_Samples> <Pathways>
      <GeneID>\n")
  stop("Provide at least first 5 arguments")
}
#### load required libraries
require(gdata)
require(e1071)
library(Biobase)
require(xtable)
library(devtools)
library(preprocessCore)
#library(rgl)
#library(qpcR)
library(data.table)
##################
Path <- args[1]
QueryMat <- readRDS(args[2])
load(args[3]) ### load R object for pathways
GeneID <- args[4]
########################
#######################
if(GeneID == "Symbol"){
  Query_GeneVec <- as.character(rownames(QueryMat))
}else if(GeneID == "EntrezGeneId"){
  Query_GeneVec <- as.numeric(rownames(QueryMat))
}
#####################
#################### Identifying target GO terms for classifcation
if(!is.na(args[5]) & args[5] != "NA"){
  MinGeneNum <- as.numeric(args[5])
}else{
  MinGeneNum <- 5
}

if(!is.na(args[6]) & args[6] != "NA"){
  MaxGeneNum <- as.numeric(args[6])
}else{
  MaxGeneNum <- 30
}
#################
GoTerms_checkedPathways <- GoTermsIDs_GeneSet[[1]]
GoIDs_checkedPathways <- GoTermsIDs_GeneSet[[2]]
GeneMatchedIDs_Go <- GoTermsIDs_GeneSet[[3]]
GeneMatchedSymbol_Go <- GoTermsIDs_GeneSet[[4]]

TargetGOind <- c()
for(PathwayIter in 1:length(GoTerms_checkedPathways)){
  PathwayGeneID <- GeneMatchedIDs_Go[[PathwayIter]]
  if(length(PathwayGeneID) >= MinGeneNum & 
     length(PathwayGeneID) <= MaxGeneNum){
    TargetGOind <- c(TargetGOind, PathwayIter)
  }
}
##################
################## function of TSC calculation
TSC_Calculation <- function(PathwayExp1, PathwayExp2){
  AA <- PathwayExp1%*%t(PathwayExp1)
  BB <- PathwayExp2%*%t(PathwayExp2)
  AA0 <- AA - diag(diag(AA))
  BB0 <- BB - diag(diag(BB))
  TSC <- sum(diag(AA0%*%BB0))/sum(AA0^2)^.5/sum(BB0^2)^.5
  return(TSC)
}
################# Determining TSC for the target pathways
TSC_Pathways <- function(RNAseqData1, RNAseqData2, GeneVec, GeneID, TargetGOind){
  
  TSC_PathwayVec <- c()
  for(PathwayIter in TargetGOind){
    if(GeneID == "Symbol"){
      PathwayGeneID <- GeneMatchedSymbol_Go[[PathwayIter]]
    }else if(GeneID == "EntrezGeneId"){
      PathwayGeneID <- GeneMatchedIDs_Go[[PathwayIter]]
    }
    
    PathwayGeneInd <- which(GeneVec %in% PathwayGeneID)
    DupInd <- which(duplicated(GeneVec[PathwayGeneInd]))
    
    if(length(DupInd) > 0){
      PathwayGeneInd <- PathwayGeneInd[-DupInd]
    }
    PathwayExp1 <- RNAseqData1[PathwayGeneInd,]
    PathwayExp2 <- RNAseqData2[PathwayGeneInd,]
    
    if(!is.null(nrow(PathwayExp1))){
      if(nrow(PathwayExp1) > 0){
        TSC <- TSC_Calculation(PathwayExp1, PathwayExp2)
        TSC_PathwayVec <- c(TSC_PathwayVec, TSC)
      }else{
        
        TSC_PathwayVec <- c(TSC_PathwayVec, 0)
      }
    }else{
      
      TSC_PathwayVec <- c(TSC_PathwayVec, 0)
    }
  }
  return(TSC_PathwayVec) 
}
########################## postprocessing TSC vectors to identify the most similar phenotype
Classifiction <- function(TSCMat){
  SimilarityVec <- apply(TSCMat, 2, function(x){median(x[which(!is.na(x) & x > 0)]
  )/mad(x[which(!is.na(x) & x > 0)])})
  MostSim <- colnames(TSCMat)[which(SimilarityVec == max(SimilarityVec))]
  return(MostSim)
}
########################### Bagging for classsification
SIGN_Bagging <- function(ExpMat_Query, RNAseqData1, GeneVec, TargetSamples,
                         UniquePhen, TarSam, BagSampleNum){
  
  PredBag <- c()
  for(BagIter in 1:BaggingNum){
    TSCMat <- c()
    for(PopulationIter in 1:length(UniquePhen)){
      RefSams <- which(TargetSamples == UniquePhen[PopulationIter] & 
                         c(1:ncol(RNAseqData1)) != TarSam)
      ExpMatrix <- matrix(as.numeric(RNAseqData1[,RefSams]),
                          ncol = length(RefSams))
      BagSamples <- sample(1:ncol(ExpMatrix),
                           BagSampleNum, replace = FALSE)
      ExpMatrix <- ExpMatrix[,BagSamples]
      
      if(ncol(ExpMatrix) == 1){
        RNAseqData2 <- matrix(rep(as.numeric(ExpMatrix),
                                  40)*runif(length(rep(as.numeric(ExpMatrix),
                                                       40)), min = 0.8, max = 1.2), ncol = 40)
      }else{
        RNAseqData2 <- ExpMatrix
      }
      TSC_PathwayVec <- TSC_Pathways(ExpMat_Query, RNAseqData2,
                                     GeneVec, GeneID, TargetGOind)
      
      TSCMat <- cbind(TSCMat, TSC_PathwayVec)
    }
    colnames(TSCMat) <- UniquePhen
    PredBag <- c(PredBag, Classifiction(TSCMat))
  }
  
  return(PredBag)
}


######################### main function to identify the most similar samples (bagging is used)
SIGN_class_1out <- function(RNAseqData1, TargetSamples,GeneVec, BaggingNum, BagSampleNum){
  
  UniquePhen <- unique(TargetSamples)
  PredAnnot <- c()
  for(SampleIter in 1:ncol(RNAseqData1)){
    print(paste("Sample_", SampleIter, sep = "", collapse = ""))
    
    ExpMat_Query <- matrix(rep(as.numeric(RNAseqData1[,SampleIter]),40)*runif(length(
      rep(as.numeric(RNAseqData1[,SampleIter]),40)), min = 0.8, max = 1.2), ncol = 40)
    
    PredBag <- SIGN_Bagging(ExpMat_Query, RNAseqData1, GeneVec,
                            TargetSamples, UniquePhen, SampleIter, BagSampleNum)

    PredBag_table <- table(PredBag)
    IdentifiedAnnot <- names(PredBag_table)[which(
      as.numeric(PredBag_table) == max(PredBag_table))]
    if(length(IdentifiedAnnot) > 1){
      print("More than 1 match")
    }
    PredAnnot <- c(PredAnnot, IdentifiedAnnot)
    
    print(PredAnnot[length(PredAnnot)])
  }
  
  return(PredAnnot)
}
####################### Determining number of samples to be used in bagging
RefSamNum <- c()
RefSamNum <- as.numeric(table(colnames(QueryMat)))
BagSamNum <- min(RefSamNum)

if(!is.na(args[7]) & args[7] != "NA"){
  BaggingNum <- as.numeric(args[7])
}else{
  BaggingNum <- 10
}
######################## 
RNAseqData1 <- as.matrix(QueryMat)
TargetSamples <- as.character(colnames(QueryMat))
colnames(RNAseqData1) <- NULL
rownames(RNAseqData1) <- NULL
RNAseqData1 <- matrix(as.numeric(RNAseqData1), ncol = ncol(RNAseqData1))

PredAnnot <- SIGN_class_1out(RNAseqData1, TargetSamples,
                             Query_GeneVec, BaggingNum, BagSamNum)

PredictionMat <- cbind(TargetSamples, PredAnnot)
colnames(PredictionMat) <- c("Actual", "Prediction")
##################### statistics of predictions
TotalTPR <- length(which(as.character(PredictionMat[,"Actual"]) ==
                           as.character(PredictionMat[,"Prediction"])))/nrow(PredictionMat)
UniquePhen <- unique(as.character(PredictionMat[,"Actual"]))
TPR_Phen <- lapply(UniquePhen,function(x){length(which(as.character(PredictionMat[,"Actual"]) ==
                                                         as.character(PredictionMat[,"Prediction"]) & 
                                                         as.character(PredictionMat[,"Actual"]) ==
                                                         x))/length(which(as.character(PredictionMat[,"Actual"])
                                                                          == x))})
names(TPR_Phen) <- UniquePhen
####################
Predictions <- list(PredictionMat, TPR_Phen, TotalTPR)
names(Predictions) <- c("PredictionMat", "TPR_Phenotypes", "TRP_Total")
saveRDS(Predictions, file = paste(Path, "/", "SIGN_ClassifiedSamples_GO", MinGeneNum,
                                  "to", MaxGeneNum, "_",BaggingNum, "Bagging", "_L1out.rds",
                                  sep = "", collapse = ""))

