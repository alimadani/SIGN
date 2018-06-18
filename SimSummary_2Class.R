SimSummary_2Class <- function(SimMat1, SimMat2){
 
  MedMad1 <- unlist(lapply(seq(1:ncol(SimMat1)), function(SimMeasure){median(as.numeric(SimMat1[,SimMeasure]))/mad(as.numeric(SimMat1[,SimMeasure]))}))
  MedMad2 <- unlist(lapply(seq(1:ncol(SimMat2)), function(SimMeasure){median(as.numeric(SimMat2[,SimMeasure]))/mad(as.numeric(SimMat2[,SimMeasure]))}))
  MedMad_Out <- ((MedMad2-MedMad1)/(MedMad2+MedMad1))

  WilcoxPaired_Out <- unlist(lapply(seq(1:ncol(SimMat1)), function(SimMeasure){wilcox.test(as.numeric(SimMat1[,SimMeasure]), as.numeric(SimMat2[,SimMeasure]), paired = TRUE)$p.value}))
  
  BubbleSort_Out <- unlist(lapply(seq(1:ncol(SimMat1)), function(SimMeasure){BubbleSort(as.numeric(SimMat1[,SimMeasure]), as.numeric(SimMat2[,SimMeasure]))}))

  Similarity_Out <- c(MedMad_Out, WilcoxPaired_Out, BubbleSort_Out)
  names(Similarity_Out) <- c(paste(rep("MedMad", ncol(SimMat1)), colnames(SimMat1), sep = "_"), paste(rep("WilcoxPaired", ncol(SimMat1)), colnames(SimMat1), sep = "_"), paste(rep("BubbleSort", ncol(SimMat1)), colnames(SimMat1), sep = "_"))

  return(Similarity_Out)
}
