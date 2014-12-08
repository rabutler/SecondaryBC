#driver program for quantile mapping of HDE data
# This code is different from the standard driver_qMap because it expects there to be
# A separate historical simulation file for each station (with the other code the historical and future were all 
# in the same matrix).  With HDE there is just one historical time series per station and not a separate historical
# time sereies for each climate projection
# Data should already be preproccesed to that the period of record in the historical simulation file matches the observed.

#Analysis steps
#--------------
#Step  0. Read observed and simulated flow time series.
#Step  1. Estimate smooth cdf from observed monthly flow time-series.
#Step  2. Estimate smooth cdf from historical simulated monthly flow time series.
#Step  3. Do quantile mapping of monthly simulated flows (historical and future). Quantile map each month for all years.
#Step  4. Add the monthly quantile mapped simulated flows for each month (Step 3) to calculate annual flows.
#Step  5. Estimate cdf from observed annual flow time-series.
#Step  6. Estimate cdf from simulated annual flow time-sereris (historical and future)
#Step  7. Lookup quantile from Step 4 values in CDF from Step 6.
#Step  8. Use quantile from Step 7 to get flow values from CDF in Step 5.
#Step  9. Get flow correction fraction from values in Steps 4 and 8.
#Step 10. Adjust monthly quantile mapped flows (Step 3) with flow correction factors from Step 9.

#Assumptions
#-----------
#Inflow units in cfs, Outflow units TAF/month
#Stationarity

#defaults
#--------
rm(list=ls())
setwd("C:/BOR_Work/UC/AAO/HDE_BiasCorrection/R_codes")
ptm <- proc.time()
source("qMap_HDE.R")

# Set the simulation information
scenarios=c("central", "hotdry", "hotwet", "warmdry", "warmwet")
periods=c(2025, 2055, 2085)
nscen=length(scenarios)
nper=length(periods)
nproj = nscen*nper+1 # The number of climate change projections (the +1 is for Maurer)
nsimyrs= 51 # The total number of years in the future simulation period (i.e. the length of the HDE time series)
simstart=1949 #Starting year for the simulation time period
outdir="../Bias_Corrected_Results/URGSIM/"
proj_list=c("2025_central", "2025_hotdry", "2025_hotwet", "2025_warmdry", "2025_warmwet",
			"2055_central", "2055_hotdry", "2055_hotwet", "2055_warmdry", "2055_warmwet",
			"2085_central", "2085_hotdry", "2085_hotwet", "2085_warmdry", "2085_warmwet", "Maurer")

# Read in the station information
station_input=matrix(scan("../Station_List.txt", what=character(0),skip=1), ncol=7, byrow=T)
station_list=station_input[,2]  #Station numbers used with VIC routing
folder=station_input[,3]    #The folder name containing the formatted observed data
ystart=as.numeric(station_input[,4])    #Starting year for observed data (note that these dates shoul match the range in the historical simualtion file)
mstart=as.numeric(station_input[,5])    #Starting month for observed data
yend=as.numeric(station_input[,6])	   #Ending year for observed data
mend=as.numeric(station_input[,7])	   #Ending month for observed data
sim_list = obs_list = station_list
nstat=length(sim_list)
print(nstat)

#Set directory locations and read in supporting files
obsdir="C:/BOR_Work/UC/AAO/HDE_BiasCorrection/Formatted_Observed/"
histdir="C:/BOR_Work/UC/AAO/HDE_BiasCorrection/Formatted_Historical_Sim/"
simdir="C:/BOR_Work/UC/AAO/HDE_BiasCorrection/Formatted_Future_Sim/"
plotdir="../Plots/URGSIM/"
daysinmth=matrix(scan("../daysinmonth.txt"),ncol=3,byrow=T) 
wy_cy=matrix(scan("../WY_CY_Conversion.txt", skip=1), ncol=3, byrow=T)


# Loop through stations
for(stat in 1:nstat){
#for(stat in 23:24){
	station=station_list[stat]
	obsfn = paste(obsdir,"Riog_ObsFlow_cfs_", station_list[stat], ".txt", sep="")
	histfn=paste(histdir,station_list[stat], "_MaurerSimulation.txt", sep="")
	ccsfn = paste(simdir, station_list[stat], "_HDESimulation.txt", sep="")
	nhistyrs=yend[stat]-ystart[stat]+1
	bc2mon=qMap(obsfn,histfn,ccsfn,nproj, nsimyrs, nhistyrs, station, ystart[stat], yend[stat],plotdir) #bc2mon has the historical with the future right now...
	print(stat)
	
	## Format monthly bias corrected output and convert from TAF to cfs
	#First for historical (note that the historical is the same for all projections so it only needs to be written once)
	temp_hist=bc2mon[1,1:(nhistyrs-1),]
	output_hist=matrix(NA, nrow=((nhistyrs-1)*12), ncol=4)
	row=1
	for(i in 1:(nhistyrs-1)){
		for(j in 1:12){
			output_hist[row,1]=station
			output_hist[row,2]= 'Flow'
			year=i+ystart[stat] +wy_cy[j,3]
			month=wy_cy[j,2]
			if(year%%4 ==0){
				day=daysinmth[month,3]} else{
				day=daysinmth[month,2]} # end if
			output_hist[row,3]=paste(year,month,day, sep="-")
			output_hist[row,4]=(temp_hist[i,j]*43560*1000)/(86400*day)
			row=row+1
		} # end for j
	} # end for i
	if(stat==1){
	write.table(output_hist, paste(outdir, "Bias_Corrected_CY_cfs.historical.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=F)
	} else {
	write.table(output_hist, paste(outdir,  "Bias_Corrected_CY_cfs.historical.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=T)
	}

	#Next for the future simulations
	for(iproj in 1:nproj){
		 temp_sim=bc2mon[iproj,(nhistyrs):(nhistyrs+nsimyrs-2),] 
		 output_sim=matrix(NA, nrow=((nsimyrs-1)*12), ncol=4) 

		#Simulation years
		row=1
		for(i in 1:(nsimyrs-1)){
			for(j in 1:12){
				output_sim[row,1]=station
				output_sim[row,2]= 'Flow'
				year=i+simstart +wy_cy[j,3]
				month=wy_cy[j,2]
				if(year%%4 ==0){
					day=daysinmth[month,3]} else{
					day=daysinmth[month,2]} # end if
				output_sim[row,3]=paste(year,month,day, sep="-")
				output_sim[row,4]=(temp_sim[i,j]*43560*1000)/(86400*day)
				row=row+1
			} # end for j
		} # end for i
	
		if(stat==1){
		write.table(output_sim, paste(outdir, proj_list[iproj], "_Bias_Corrected_CY_cfs.simulated.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=F)
		} else {
		write.table(output_sim, paste(outdir, proj_list[iproj], "_Bias_Corrected_CY_cfs.simulated.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=F, append=T)
		}
	} # end for iproj
} # end station loop

ptm=proc.time() - ptm
print(ptm)