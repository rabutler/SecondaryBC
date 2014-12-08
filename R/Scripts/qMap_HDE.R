"qMap" <- function(obsfn,histfn,ccsfn,nproj,nsimyrs,nhistyrs, station, histstart, histend,plotdir){
#This code is different from the standard qMAP file because it reads in separte files for the historical
# simulation period and the future simulation period.  Files should be preprocessed so tha the historical simulation
# and observed date ranges are the same. 
# Note that the bias corrected flow outputs (bc1mon and bc2mon) still contain both the historical and future
# simulation data appended together.  However, because the historical simulation is the same for all of the projections
# The outputs are repeated nprojection times for the historical years. As such all plotting etc. of historical bias
# corrected data just grabs the first projection since the results will be the same for all. 

#defaults
source("estimatesmoothcdf.R")
source("calcdf.R")
source("zfill.R")
library(fields)
flush.console()


##### GHANGE PLOTFLAGS
plotflag_month_hydro=FALSE  # seprate hydrographs by month and projection
plotflag_ann_hydro=FALSE	# Annual hyrdrograph with lines for observed and simulated
plotflag_cdf=FALSE			# separate cdf plots by month
plotflag_hydro_barplot= TRUE # TRUE Annual hydrograph with lines for observed, simulated, and both bias correction steps and barplot of annual average flow for overlap time period with bars for observed, simulated and bias corrected
plotflag_timeseries = TRUE   #  TRUE timeseries plot split into 4 sections with lies for observed, simulated and bias corrected
daysinmth=matrix(scan("../daysinmonth.txt"),ncol=3,byrow=T) 

#I. Assumes Stationarity - no change in distribution parameters between historical and future simulations
#read flow time-series
# Assumes all flows are in cfs
obsdat=matrix(scan(obsfn),ncol=3,byrow=T)    #year,month,observedflow
histdat=matrix(scan(histfn),ncol=3,byrow=T)    #year,month,historicalSimulatedflow
simdat=matrix(scan(ccsfn, skip=1),ncol=(nproj+2),byrow=T)  #year,month,proj1 ... nproj

#convert flows to TAF (thousand acre-ft)
obsflow=obsdat[,3]
histflow=histdat[,3]
for(i in 1:(nhistyrs*12)){
    year=obsdat[i,1]
	month=obsdat[i,2]
	if(year%%4 ==0){
		days=daysinmth[month,3]} else {days=daysinmth[month,2]}
	    obsflow[i]=(obsflow[i]*86400*days)/(43560*1000) #CFS to TAF
	    histflow[i]=(histflow[i]*86400*days)/(43560*1000) #CFS to TAF
} # end for i
obsmatcy=matrix(obsflow,ncol=12,byrow=T)                 #50 x 12 matrix of monthly values, calendar year basis
histmatcy=matrix(histflow,ncol=12,byrow=T)                 #50 x 12 matrix of monthly values, calendar year basis

simflow=simdat[,3:(nproj+2)]
for(i in 1:(nsimyrs*12)){
    year=simdat[i,1]
	month=simdat[i,2]
	if(year%%4 ==0){
		days=daysinmth[month,2]} else{days=daysinmth[month,3]}
		simflow[i,]=(simflow[i,]*86400*days)/(43560*1000) #CFS to TAF
} # end for i

simmatcy=array(NA,dim=c(nproj,nsimyrs,12))
for (iproj in 1:nproj){
  simmatcy[iproj,,]=matrix(simflow[,iproj],ncol=12,byrow=T)  #112 x 150 x 12 matrix of monthly values, calendar year basis
} #iproj

#water year basis
obsmatwy=matrix(NA,nrow=(nhistyrs-1),ncol=12)             #49 x 12 values, oct...sep values
obsmatwy[1:(nhistyrs-1),1:3]=obsmatcy[1:(nhistyrs-1),10:12]
obsmatwy[1:(nhistyrs-1),4:12]=obsmatcy[2:nhistyrs,1:9]

histmatwy=matrix(NA,nrow=(nhistyrs-1),ncol=12)             #49 x 12 values, oct...sep values
histmatwy[1:(nhistyrs-1),1:3]=histmatcy[1:(nhistyrs-1),10:12]
histmatwy[1:(nhistyrs-1),4:12]=histmatcy[2:nhistyrs,1:9]

simmatwy=array(NA,dim=c(nproj,(nsimyrs-1),12))
simmatwy[,1:(nsimyrs-1),1:3]=simmatcy[,1:(nsimyrs-1),10:12]
simmatwy[,1:(nsimyrs-1),4:12]=simmatcy[,2:nsimyrs,1:9]

obsannwy=apply(obsmatwy,1,"sum")
histannwy=apply(histmatwy,1,"sum")
simannwy=matrix(NA,nproj,(nsimyrs-1))
for (iproj in 1:nproj){
  simannwy[iproj,]=apply(simmatwy[iproj,,],1,"sum") #112 x 149, each projection along row and each wy along columns
} #

#observed time-series, water years 1951-1999
#-------------------------------------------
#develop monthly CDFs using the 49 water year values, the column values are from oct-sep
monobscdf=array(NA,dim=c(12,50,2)) #12, months; 50, number of points used in smoothing the cdfs, see R code, estimatesmoothcdf.R
annobscdf=matrix(NA,nrow=50,ncol=2) #50, number of points used in smoothing the cdfs, see R code, estimatesmoothcdf.R
month=c("oct","nov","dec","jan","feb","mar","apr","may","jun","jul","aug","sep")

for(im in 1:12){
  mondata=obsmatwy[,im] #1=oct, 2=nov, ...
  if(max(mondata)==0){
  print(paste("month", im, "has all zeros for station", station))
  monobscdf[im,,1]= rep(0,50)
  monobscdf[im,,2]=seq(from=0, to=1, length.out=50)
  } else {
  smcdf=estimatesmoothcdf(mondata) #smcdf -- variable used to temporary store the estimated smooth cdf
  monobscdf[im,,1]=smcdf$x
  monobscdf[im,,2]=smcdf$y
  }
  #ecdf=calcdf(mondata)
  #windows(6,6)
  #plot(ecdf$x,ecdf$y,type="l",xlab="flow (taf)",ylab="cumulative probability")
  #lines(smcdf$x,smcdf$y,col="red")
  #title(main=month[im])
  #fout=paste(month[im],".pdf",sep="")
  #dev.print(pdf,file=fout)

  #fout=paste(im,".emf",sep="")
  #savePlot(fout,type="emf")
} #im

smcdf=estimatesmoothcdf(obsannwy)
annobscdf[,1]=smcdf$x
annobscdf[,2]=smcdf$y


#simulated time-series for the observed time period
#---------------------------------------------------------------
monsimcdf=array(NA,dim=c(12,50,2)) ###CHANGE

for (im in 1:12){
    #print(im)
    mondata=histmatwy[1:(nhistyrs-1),im] #1:(nhistyrs-1), water years, 1951-1999
	if(max(mondata)==0){
		print(paste("month", im, "has all zeros for station", station))
		monsimcdf[im,,1]= rep(0,50)
		monsimcdf[im,,2]=seq(from=0, to=1, length.out=50)
		} else {
    smcdf=estimatesmoothcdf(mondata)
    monsimcdf[im,,1]=smcdf$x
    monsimcdf[im,,2]=smcdf$y
	}
} #im
  
  #anndata=simannwy[iproj,1:(nhistyrs-1)] #1:(nhistyrs-1), water years, 1951-1999
  #smcdf=estimatesmoothcdf(anndata)
  #annsimcdf[iproj,,1]=smcdf$x
  #annsimcdf[iproj,,2]=smcdf$y


## test cdf plot
#for(im in 1:12){
#  xrange=range(c(monobscdf[im,,1],monsimcdf[iproj,im,,1])) #im=1, oct, etc.
#  plot(monobscdf[1,,],xlim=xrange,ylim=c(0,1),col="red",type="l",lwd=2,
#       xlab="flow (taf)",ylab="cumulative probability")
#  title(main=month[im])
#  lines(monsimcdf[iproj,5,,],col="light blue")
#  fout=paste(station, month[im],".pdf",sep="")
#  dev.print(pdf,file=fout)
#  print(max(monsimcdf[iproj,im,,1]))
#} #im
###


if (plotflag_cdf){
 par(mfrow=c(3,4))
for(im in 1:12){
 
  xrange=range(c(monobscdf[im,,1],monsimcdf[im,,1])) #im=1, oct, etc.
  plot(monobscdf[1,,],xlim=xrange,ylim=c(0,1),col="red",type="l",lwd=2,
       xlab="flow (taf)",ylab="cumulative probability")
  title(main=month[im])
  lines(monsimcdf[im,,],col="light blue")
  #fout=paste(station, month[im],".pdf",sep="")
  #dev.print(pdf,file=fout)

} #im

} #plotflag_cdf

#apply monthly quantile maps to perform first-pass adjustment
#------------------------------------------------------------

bc1mon=array(NA,dim=c(nproj,(nhistyrs+nsimyrs-2),12)) #bias corrected flows for the historical and future time periods
bc1ann=array(NA,dim=c(nproj,(nhistyrs+nsimyrs-2)))
errorL_matrix=errorH_matrix=matrix(0, nrow=nproj, ncol=12)

for (iproj in 1:nproj){
  for (im in 1:12){            #Oct-Sep
    mondata=c(histmatwy[,im], simmatwy[iproj,,im])  #Combine the historical and simulated flows
    z=approx(monsimcdf[im,,1],monsimcdf[im,,2],mondata,yleft=-Inf,yright=Inf)
    qsimcdf=z$y #quantile from simulated cdf
    minobs=monobscdf[im,1,1]; maxobs=monobscdf[im,50,1]  #min index =1; max index = 50, since 50 points are used to smooth the cdfs
    minsim=monsimcdf[im,1,1]; maxsim=monsimcdf[im,50,1]
    z=approx(monobscdf[im,,2],monobscdf[im,,1],qsimcdf,yleft=-Inf,yright=Inf) #note that the input is quantile and the output is flow
                                                                       #for some cases in the lower bound, the cumulative probability in qsimcdf
	  																  #could be smaller than the cumulative probability value in monobscdf.
																	  #The upper bound of the cumulative probability in the smooth cdf estimation
																	  #is always equal to one.
									      
    bc1flow=z$y  #there will be -Inf and Inf here

    maxindexlist=which(bc1flow == Inf)   #magnitude greater than overlap period simulated flows
    minindexlist=which(bc1flow == -Inf)  #magnitude less than overlap period simulated flows
    
	
	##### START CHANGES TO UNCONSTRAINED TAILS LEC#######
	# Find the 5%  and 95% value for the bc flows that were mapped (i.e. excluding the -Inf and Inf)
	outliers=c(maxindexlist, minindexlist)
	if(length(outliers)>0){threshold_l=quantile(bc1flow[-outliers], 0.05)}else{threshold_l=9999}
    if(length(outliers)>0){threshold_h=quantile(bc1flow[-outliers], 0.95)}else{threshold_h=9999}
	
	# Find the max and min tail values that would result from unconstrained mapping
	if(length(minindexlist)>0){test_low= max((minobs/minsim)*mondata[minindexlist])} else {test_low=-9999}
	if(length(maxindexlist)>0){test_hi= min((maxobs/maxsim)*mondata[maxindexlist])} else {test_hi=9999999}
	
	# Now deal with the tails (i.e. the values that got mapped to -INF or +INF)
	#	- If unconstrained mapping values don't cross the 5th or 95th percentile values then proceed with the 
	#	- the standard unconstrained approach.  If they do cross then the new mapping will be constrained such
	# 	- that the new ratio is the max(min) of the mapped values (excluding the +/- Inf vales) divided by the 
	#   - min(max) of the original tail values. 
	
	# First look at the top values
	#Test if the scaled tails cross the 95% threshold of the mapped values
	if(test_hi>= threshold_l){
		# if the test is greater than the threshold calculate the tail the normal way	
		bc1flow[maxindexlist]=(maxobs/maxsim)*mondata[maxindexlist]
    } else{
		# Calcualte the new scaling factor (ratio of the minimum mapped flow to the maximum tail flow)
		ratio=max(bc1flow[-outliers])/min(mondata[maxindexlist])
		bc1flow[maxindexlist]=ratio*mondata[maxindexlist] 

		#Calculate the percentage of mapped values that would fall above the most extreme unconstrained
		# tail value (this is just a diagnostic to show how far off the unconstrained option is)
		test=min((maxobs/maxsim)*mondata[maxindexlist])
		position=length(which(bc1flow[-outliers]>test))/length(bc1flow[-outliers])*100
		error_matrixH[iproj, im]=position
	}

	# Second look at the bottom values
	#Test if the scaled tails cross the 5% threshold of the mapped values
	if(test_low <= threshold_l){
		# if the test is less than the threshold calculate the tail the normal way
		bc1flow[minindexlist]=(minobs/minsim)*mondata[minindexlist]
	} else{
		# If not then the scaling factor must be adjusted
		# The new scaling factor is the ratio of the minimum mapped flow to the maximum tail flow
		ratio=min(bc1flow[-outliers])/max(mondata[minindexlist])
		bc1flow[minindexlist]=ratio*mondata[minindexlist] 

		#Calculate the percentage of mapped values that would fall below the most extreme unconstrained
		# tail value (this is just a diagnostic to show how far off the unconstrained option is)
		test=max((minobs/minsim)*mondata[minindexlist])
		position=length(which(bc1flow[-outliers]<test))/length(bc1flow[-outliers])*100
		errorL_matrix[iproj, im]=position	
		##### END CHANGES TO UNCONSTRAINED TAILS LEC#######	
	}
   
   bc1mon[iproj,,im]=bc1flow

    #print(c(iproj,range(bc1flow)))
   
    if (plotflag_month_hydro){
    #test plots
    z1=estimatesmoothcdf(obsmatwy[,im])
    z2=estimatesmoothcdf(mondata[1:(nhistyrs-1)])
    z3=estimatesmoothcdf(bc1flow[1:(nhistyrs-1)])
    plot(z1$x,z1$y,xlim=range(c(z1$x,z2$x,z3$x)),type="l",xlab="flow (taf)",ylab="cumulative probability")
    lines(z2,col="red")
    lines(z3,col="blue")
    title(main=paste(month[im],"wy 1959-2002",sep="\n"))
    legend("bottomright",legend=c("obs","sim","bcf"),lty=c(1,1,1),col=c("black","red","blue"))

    ##fout=paste(month[im],".pdf",sep="")
    ##dev.print(pdf,file=fout)

    #fout=paste(station, month[im],"_",zfill(iproj),".emf",sep="")
    #savePlot(fout,type="emf")

    #dev.off()
    } #plotflag
    
  } #im
  #get annual total from first pass
  bc1ann[iproj,]=apply(bc1mon[iproj,,],1,"sum")

} #iproj
	print("Average % of mapped values which are more than the min tail with origial scale")
	print(round(apply(errorH_matrix,2,mean),1))
	print("Average % of mapped values which are less than the max tail with origial scale")
	print(round(apply(errorL_matrix,2,mean),1))

if (plotflag_ann_hydro){
#DEBUG, plot annual hydrographs for water years 1951-1999
xobs=apply(obsmatwy,2,"mean"); plot(xobs,type="l",xlab="month (oct-sep)",ylab="flow (taf)")
title(main=paste("Mean Monthly Flows","Water Years 1951-1999",sep="\n"))
X=array(NA,dim=c(nproj,12))
for (iproj in 1:nproj){
  xx=bc1mon[iproj,1:(nhistyrs-1),]
  X[iproj,]=apply(xx,2,"mean")
}
xsim=apply(X,2,"mean")
lines(xsim,col="red")
} #plotflag


#bias correct annual flows after the first pass
#----------------------------------------------

bc1anncdf=array(NA,dim=c(50,2))

#Get annual cdfs from bias corrected historical flow
#The hitorical data is the same for all projections so just
#use one to get the cdf no need to loop through projections
anndata=bc1ann[1,1:(nhistyrs-1)] 
smcdf=estimatesmoothcdf(anndata)
bc1anncdf[,1]=smcdf$x
bc1anncdf[,2]=smcdf$y


bc2ann=array(NA,dim=c(nproj,(nhistyrs+nsimyrs-2)))

for (iproj in 1:nproj){
  anndata=bc1ann[iproj,]
  z=approx(bc1anncdf[,1],bc1anncdf[,2],anndata,yleft=-Inf,yright=Inf)
  qsimcdf=z$y #quantile from simulated cdf
  minobs=annobscdf[1,1]; maxobs=annobscdf[50,1]  #min index =1; max index = 50, since 50 points are used to smooth the cdfs
  minsim=bc1anncdf[1,1]; maxsim=bc1anncdf[50,1]
  z=approx(annobscdf[,2],annobscdf[,1],qsimcdf,yleft=-Inf,yright=Inf) #note that the input is quantile and the output is flow
                                                                      #for some cases in the lower bound, the cumulative probability in qsimcdf
								      #could be smaller than the cumulative probability value in observed cdf.
  bc2flow=z$y  #there will be -Inf and +Inf here
  maxindexlist=which(bc2flow == Inf)   #magnitude greater than overlap period simulated flows
  minindexlist=which(bc2flow == -Inf)  #magnitude less than overlap period simulated flows

  bc2flow[maxindexlist]=(maxobs/maxsim)*anndata[maxindexlist]
  bc2flow[minindexlist]=(minobs/minsim)*anndata[minindexlist]
    
  bc2ann[iproj,]=bc2flow

} #iproj

#final monthly values
#--------------------
bc2mon=array(NA,dim=c(nproj,(nhistyrs+nsimyrs-2),12))
for (iproj in 1:nproj){
  for (iyr in 1:(nhistyrs+nsimyrs-2)){
	if(bc1ann[iproj,iyr]> 0){    #CHANGED LEC 9/20/11: Added the if else loop to avoid dividing by zero and getting a NAN when bc1ann is zero
		bc2mon[iproj,iyr,]=bc1mon[iproj,iyr,]*(bc2ann[iproj,iyr]/bc1ann[iproj,iyr])
	} else {
		bc2mon[iproj,iyr,]=bc1mon[iproj,iyr,]
	}
  } #iyr
} #iproj

#DEBUG, diagnostic plots -

#Once again all of the historical data is the same so just
# use the historical data in with the first projection to calculate no need to loop through projections
X1 <- X2 <- X0 <- array(NA,dim=c(1,12))
xx1=bc1mon[1, 1:(nhistyrs-1),]
X1[1,]=apply(xx1,2,"mean")
xx2=bc2mon[1,1:(nhistyrs-1),]
X2[1,]=apply(xx2,2,"mean")
xx0=simmatwy[1,1:(nhistyrs-1),]
X0[1,]=apply(xx0,2,"mean")


xobs=apply(obsmatwy,2,"mean");
xsim0=apply(X0,2,"mean");
xsim1=apply(X1,2,"mean");
xsim2=apply(X2,2,"mean");



##  plot annual hydrographs for water years 1950-1999
if (plotflag_hydro_barplot){
par(mfrow=c(1,2))
# plot hydrograph
par(xaxt="n")
plot(xobs,type="l", lwd=2, ylim=range(c(xobs,xsim0,xsim2)),
     xlab="",ylab="TAF")
title(main=paste("Mean Monthly Volumes", "\n", "Water Years", histstart+1, "-", histend,sep=""))
lines(xsim0,lwd=2, col=2)
lines(xsim2,col=13, lty=4, lwd=2)
par(xaxt="s"); axis(1,at=1:12,labels=month)
legend("topleft",c("OBS","SIM","BCF"),lwd=c(2,2,2),lty=c(1,1,4), col=c(1,2,13))

# Make barplot
sim_ann_avg=mean(simannwy[,1:(nhistyrs-1)])
obs_ann_avg=mean(obsannwy)
bcf_ann=apply(bc2mon[1,1:(nhistyrs-1),],1,"sum")
bcf_ann_avg = mean(bcf_ann)
plot_data=c(obs_ann_avg, sim_ann_avg, bcf_ann_avg)

barplot(plot_data, names.arg=c("OBS", "SIM", "BCF"), ylim=c(0,max(plot_data)*1.2), col=c(1,2,13), ylab="TAF")
title(main=paste("Annual Mean Volumes", "\n", "Water Years", histstart+1, "-", histend,sep=""))
box()

# Write to file
fout=paste(plotdir, station, "_Hydrograph_Barplot.emf",sep="")
    savePlot(fout,type="emf")
}

#Make timeseries plot
if (plotflag_timeseries){
obs_ts = sim_ts = bc2mon_ts = time_ts =rep(NA, ((nhistyrs-1)*12)) 
k=1
for(i in 1:(nhistyrs-1)){
	for(j in 1:12){
		bc2mon_ts[k] = bc2mon[1,i,j]
		#sim_ts[k] = simmatwy[1,i,j]
		sim_ts[k] = histmatwy[i,j]
		obs_ts[k] = obsmatwy[i,j]
		time_ts[k]=i+(histstart-1)+j/12
		k=k+1
		} # end for j
	} # end for i

#Divide the timeseris into four sections
nplot=round((nhistyrs -1)/4,0)*12	
nlast=(nhistyrs-1)*12 - (nplot*3)

# plot in sections
###START HERE --- NEED TO FIGURE OUT HOW TO MAKE THIS PLOT FOR DIFFERENT YEAR RANGES
par(mfrow=c(5,1))
par(mar = c(1, 4, 2.5, 1)) 
tstitle=paste("Monthly Volumes, WY ", histstart+1, "-", histend, sep="")
t1 = 1
t2 = nplot
range_y=range(c(obs_ts, sim_ts, bc2mon_ts))
plot(time_ts[t1:t2], obs_ts[t1:t2], ylim=range_y, col=1, lwd=2, type="l", xlab="Water Year", ylab="TAF", main=tstitle)
lines(time_ts[t1:t2],sim_ts[t1:t2], col=2)
lines(time_ts[t1:t2], bc2mon_ts[t1:t2], col=13, lty=4)

t1 = nplot+1
t2 = 2*nplot
plot(time_ts[t1:t2], obs_ts[t1:t2], ylim=range_y, col=1, lwd=2, type="l", xlab="Water Year", ylab="TAF")
lines(time_ts[t1:t2],sim_ts[t1:t2], col=2)
lines(time_ts[t1:t2], bc2mon_ts[t1:t2], col=13, lty=4)

t1 = 2*nplot+1
t2 = 3*nplot
plot(time_ts[t1:t2], obs_ts[t1:t2], ylim=range_y, col=1, lwd=2, type="l", xlab="Water Year", ylab="TAF")
lines(time_ts[t1:t2],sim_ts[t1:t2], col=2)
lines(time_ts[t1:t2], bc2mon_ts[t1:t2], col=13, lty=4)

t1 = 3*nplot+1
t2 = (nhistyrs-1)*12
plot(time_ts[t1:t2], obs_ts[t1:t2], ylim=range_y, col=1, lwd=2, type="l", xlab="Water Year", ylab="TAF")
lines(time_ts[t1:t2],sim_ts[t1:t2], col=2)
lines(time_ts[t1:t2], bc2mon_ts[t1:t2], col=13, lty=4)	

plot(c(0,5),c(0,5),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
legend("top",c("OBS","SIM","BCF"),lwd=c(2,1,1),lty=c(1,1,2), col=c(1,2,13), ncol=3, bty='n')

fout=paste(plotdir,station, "_Timeseries.emf",sep="")
    savePlot(fout,type="emf")
}

	
return (bc2mon)
#return(sim_ts)
} #end function

