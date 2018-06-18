BubbleSort <- function(Vec1,Vec2){
  if(length(Vec1)!=length(Vec2)){
    stop("permutations need to be of same length!")
  }
  
  Vec1_rank <- rank(Vec1)
  Vec2_rank <- rank(Vec2)

  return(((choose(length(Vec1_rank),2) - cov(Vec1_rank,Vec2_rank,method='kendall')/2)/2)/choose(length(Vec1_rank),2))
}
