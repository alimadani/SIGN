SIGN_Aggregate <- function(ScoreList, TimeList, EventList){

 MethodNum <- length(ScoreList)
 ScoreMat_Agg <- c()
 for(MethodIter in 1:MethodNum){
  ScoreList_Path <- ScoreList[[MethodIter]]
  PathwayNum <- length(ScoreList_Path)
  print(PathwayNum)
  if(PathwayNum > 0){
   for(PathwayIter in 1:PathwayNum){
    ScoreList_TMP <- ScoreList_Path[[PathwayIter]]
    print(length(ScoreList_TMP))
    ScoreMat_TMP <- c()
    for(ClassIter in 1:length(ScoreList_TMP)){
     ScoreMat_TMP <- rbind(ScoreMat_TMP, ScoreList_TMP[[ClassIter]])
     print(dim(ScoreList_TMP[[ClassIter]]))
    }
    print(dim(ScoreMat_TMP))
    colnames(ScoreMat_TMP) <- paste(rep(names(ScoreList)[MethodIter], length(colnames(ScoreMat_TMP))),
                                    rep(names(ScoreList_Path)[PathwayIter], length(colnames(ScoreMat_TMP))),
                                    colnames(ScoreMat_TMP), sep = "_")
    ScoreMat_Agg <- cbind(ScoreMat_Agg, ScoreMat_TMP)
   }
  }
 }

 TimeVec <- c()
 EventVec <- c()
 for(ClassIter in 1:length(TimeList)){
  TimeVec <- c(TimeVec, TimeList[[ClassIter]])
  EventVec <- c(EventVec, EventList[[ClassIter]])
 }

 SIGN_Agg_Out <- list(ScoreMat_Agg, TimeVec, EventVec)
 names(SIGN_Agg_Out) <- c("Scores", "Time", "Event")

 return(SIGN_Agg_Out)
}
