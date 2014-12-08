cdf = function(x,xx1){
#function to calulate the CDF of a vector..
#source("code/cdf.r")

#Estimate the CDF at N points on a grid
#N=500
#x=rnorm(N)
#x=flowData[,2]
#xx1=seq(min(x)-sd(x), max(x)+sd(x),length=N) #replace with eval.points
N=length(xx1)

#for data with 0 as lower bound..
xx=c(0,xx1) 
#for data with -infinity as lower bound..
#xx=c(-999,xx1)

#source("code/qgaus.r")
#source("code/integrand-np.r")

zz=1:N
N1=N+1
for(i in 2:N1){
	x1=xx[i-1]
	x2=xx[i]
	zz[i-1]=qgaus(x1,x2,x)
}
zz=cumsum(zz/sum(zz))

#plot(xx[2:(N+1)],zz,xlab="X",ylab="CDF")

return(zz)
}
