EventRenaming <- function(EventVec, Censored_Annot){

 CensoredInd <- which(EventVec == Censored_Annot)

 EventVec_Renamed <- EventVec
 EventVec_Renamed[CensoredInd] <- 0
 EventVec_Renamed[-CensoredInd] <- 1 

 return(EventVec_Renamed)
}
