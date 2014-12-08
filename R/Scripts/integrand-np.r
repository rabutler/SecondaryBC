integrand=function(xeval,x){

#fuction to be evaluated for integration..

# the NP estimator using Normal kernel..
library(sm)
band=hnorm(x)

n=length(x)

tval=(xeval-x)/band
xdens=sum(dnorm(tval,0,1))

xdens=xdens/(n*band)

xdens

}
