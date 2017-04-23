
#----- private functions ------

#' This function maps a vector of p-values to the open unit interval 
#' (that is, it moves them away from 0 and 1).  
#' Needed for input to p-value combining function.
#' @keywords internal
#' @example 
open01<-function (p, B) 
{
  (p + 1/(2 * B))/(1 + 1/B)
}

#' combine p-vals using Fisher's method
#' @keywords internal
#' @example 
fisherComb<-function(p, B) 
{
  if(missing(B)){
    comb=-2 * sum(log(p))
  }else{
    comb=-2 * sum(log(open01(p, B)))
  }
  1-pchisq(comb,df=2*length(p))
}


