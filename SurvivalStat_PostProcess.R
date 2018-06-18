SurvivalStat_PostProcess <- function(StatList){
 
 Cindex <- c()
 Cindex_std <- c()
 LogTest_pval <- c()
 NamesVec <- c()
 for(Iter in 1:length(StatList)){
  CoxList <- StatList[[Iter]]
  if(class(CoxList) == "summary.coxph"){
   Cindex <- c(Cindex, as.numeric(CoxList$concordance[1]))
   Cindex_std <- c(Cindex_std, as.numeric(CoxList$concordance[2]))
   LogTest_pval <- c(LogTest_pval, as.numeric(CoxList$logtest[3]))
   NamesVec <- c(NamesVec, names(StatList)[Iter])
  }else if(class(CoxList) == "list"){
   for(ModelIter in 1:length(CoxList)){
    CoxModel <- CoxList[[ModelIter]]
    Cindex <- c(Cindex, as.numeric(CoxModel$concordance[1]))
    Cindex_std <- c(Cindex_std, as.numeric(CoxModel$concordance[2]))
    LogTest_pval <- c(LogTest_pval, as.numeric(CoxModel$logtest[3]))
    NamesVec <- c(NamesVec, names(StatList[[Iter]])[ModelIter])
   }
  }
 }

 names(Cindex) <- NamesVec
 names(Cindex_std) <- NamesVec
 names(LogTest_pval) <- NamesVec
 
 StatList_Summary <- list(Cindex, Cindex_std, LogTest_pval)
 names(StatList_Summary) <- c("Cindex", "Cindex_std", "LogTest_pval")
 
 return(StatList_Summary)
}
