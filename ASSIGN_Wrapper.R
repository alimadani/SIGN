ASSIGN_Wrapper <- function(ExpMat_Test, ExpMat_Ref1, ExpMat_Ref2, GeneVec, PathwaySet){

 ASSIGN_Vec1 <- ASSIGN_Calculation(ExpMat_Test, ExpMat_Ref1, GeneVec, PathwaySet)
 ASSIGN_Vec2 <- ASSIGN_Calculation(ExpMat_Test, ExpMat_Ref2, GeneVec, PathwaySet)

 ASSIGN_SimOut <- c(cor.test(ASSIGN_Vec1, ASSIGN_Vec2)$estimate,
                   wilcox.test(ASSIGN_Vec1, ASSIGN_Vec2, paired = TRUE)$p.value,
                   (median(ASSIGN_Vec1)/mad(ASSIGN_Vec1)-median(ASSIGN_Vec2)/mad(ASSIGN_Vec2))/(median(ASSIGN_Vec1)/mad(ASSIGN_Vec1)+median(ASSIGN_Vec2)/mad(ASSIGN_Vec2)))

 names(ASSIGN_SimOut) <- c("ASSIGN_Pearson", "ASSIGN_WilcoxPaired", "ASSIGN_MedianMad")
 return(ASSIGN_SimOut)
}
