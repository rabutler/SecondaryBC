estimatesmoothcdf <- function(data)

{

library(sm)

source("qgaus.r")
source("integrand-np.r")
source("cdf2.r")

#source("calcdf.R")

#bandd=hsj(data)		#Sheather Jones approach for bandwidth, SG_20110323
#bandd=hnorm(data)		#for the Normal bandwidth..

#go one bandwidth from the min and max of the data..
#xlow=min(data)-bandd
xlow=min(data)        #SG_20110323


#xlow=max(0,xlow)
#xhigh=max(data)+bandd
xhigh=max(data)      #SG_20110323


#create 50 points equally spaced between xlow and xhigh..
xeval=seq(xlow,xhigh,length=50)

neval=length(xeval)

#smoothed cdf of the data vector
xcdf=cdf(data,xeval)

#empirical (Weibull) cdf can also be calculated using the function calcdf
#z=calcdf(data)

#plot(ecdf(data),xlim=range(data,xeval),main="")
#lines(xeval,xcdf,col="red")

#title(main="Empirical and Smoothed CDF")

return(list(x=xeval,y=xcdf))

} #end function