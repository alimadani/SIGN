TSC <- function(PathwayExp1, PathwayExp2){
  AA <- PathwayExp1%*%t(PathwayExp1)
  BB <- PathwayExp2%*%t(PathwayExp2)
  AA0 <- AA - diag(diag(AA))
  BB0 <- BB - diag(diag(BB))
  TSC <- sum(diag(AA0%*%BB0))/sum(AA0^2)^.5/sum(BB0^2)^.5
  return(TSC)
}
