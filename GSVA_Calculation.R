
GSVA_Calculation <- function(ExpMat1, ExpMat2, GeneVec, GeneSets, Name){

 rownames(ExpMat1) <- GeneVec
 rownames(ExpMat2) <- GeneVec
 
 ExpMat_Combined <- cbind(ExpMat1, ExpMat2) 
 GSVA_TMP <- invisible(gsva(ExpMat_Combined, GeneSets, method='gsva'))

 GSVAMat1 <- matrix(GSVA_TMP$es.obs[,1], ncol = 1)
 GSVAMat2 <- matrix(GSVA_TMP$es.obs[,2:ncol(ExpMat_Combined)], ncol = (ncol(ExpMat_Combined)-1)) 

 Pearson_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(cor.test(as.numeric(X),as.numeric(Y))$estimate)})})
 BubbleSort_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(BubbleSort(as.numeric(X),as.numeric(Y)))})})
 WilcoxonPaired_Vec <- apply(GSVAMat1, 2, function(X){apply(GSVAMat2, 2, function(Y){as.numeric(wilcox.test(as.numeric(X),as.numeric(Y), paired=TRUE)$p.value)})})

 GSVA_output <- c(median(Pearson_Vec), median(BubbleSort_Vec), (1 - pchisq( -2*sum(log(na.omit(WilcoxonPaired_Vec))), 2*length(na.omit(WilcoxonPaired_Vec)))))
 names(GSVA_output) <- paste(Name, c("GSVA_Pearson", "GSVA_BubbleSort", "GSVA_WilcoxonPaired"), sep = "_")

 return(GSVA_output)
}
