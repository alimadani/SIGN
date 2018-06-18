############  Calculating summary metric for the prediction of survival of samples
############  Plot the survival curve based on the predicted classes
Survival_KMplot <- function(TimeVec, EventVec, Class, MaxTime, OutDir, Name){

 surv.obj <- survfit(Surv(TimeVec,EventVec) ~ Class)
 dind <- D.index(x=Class, surv.time=TimeVec, surv.event=EventVec,na.rm = TRUE)
 print(dind$p.value)

 pdf(paste(OutDir, Name, ".pdf",
           sep = "", collapse = ""), width=5,height=5)
 par(mar=c(5,5,2,3)+0.1,mgp=c(3,1,0))
 plot(main = "", surv.obj,col =c("blue", "red"),
      lty = 1,lwd = 3, xlim = c(0,MaxTime), cex.lab = 2, cex.axis = 2,
      xlab = "Time (years)",ylab = "Probability of Overall Survival")
 legend("bottomleft",fill = c("red", "blue"), cex = 2,
        legend = c("Positive", "Negative"),bty = "n")
 dev.off()

}
