ExpPhen_Matching <- function(ExpMat, MetaMat, SamID_Meta){
 
 SamIntersect <- intersect(as.character(colnames(ExpMat)), as.character(MetaMat[,SamID_Meta]))

 ExpMat_Matched <- ExpMat[,which(as.character(colnames(ExpMat)) %in% SamIntersect)]
 MetaMat_Matched <- MetaMat[which(as.character(MetaMat[,SamID_Meta]) %in% SamIntersect),]

 ExpMat_Matched <- ExpMat_Matched[,sort(as.character(colnames(ExpMat_Matched)), index.return = T)[[2]]]
 MetaMat_Matched <- MetaMat_Matched[sort(as.character(MetaMat_Matched[,SamID_Meta]), index.return = T)[[2]],]

 ExpMeta_List <- list(ExpMat_Matched, MetaMat_Matched)
 names(ExpMeta_List) <- c("Expression", "Meta")

 return(ExpMeta_List)

}
