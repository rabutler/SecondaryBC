"calcdf" <- function(x)
{
 x <- sort(x)
 n <- length(x)
 y <- rep(0,n)
 for (i in 1:n) {
  y[i]=i/(n+1)
 }
 cdf=list("x"=x,"y"=y)
 return(cdf)
} 
