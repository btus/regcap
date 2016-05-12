#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# ##For daily averages
# Rscript --vanilla REGCAP_rco_viewer.R <directory_containing_data_files> days

# #For hourly averages
# Rscript --vanilla REGCAP_rco_viewer.R hours

# #For monthly averages
# Rscript --vanilla REGCAP_rco_viewer.R months

#Or combine them:
# Rscript --vanilla REGCAP_rco_viewer.R months days hours

#follow with:
# mv *.months [new_directory]
# mv *.hours [new_directory]

library(lubridate); library(xts); library(data.table); library(zoo)

setwd(args[1])

#Set working directory
#setwd(args[1])

endCycleReassign=function(x){
	end_ind<-which(x==102)
	for(i in 1:length(end_ind)){
		x[end_ind[i]]<-x[(end_ind[i]-1)]
	}
	return(x)
}

AHUflagIndex=function(x){
	AH_none<-x[which(x==0)]
	AH_none[AH_none==0]<-1
	AH_heat<-x[which(x==1)]
	AH_cool<-x[which(x==2)]
	AH_cool[AH_cool==2]<-1
	AH_vent<-x[which(x==100)]
	AH_vent[AH_vent==100]<-1
	all<-merge.xts(x, AH_none, AH_heat, AH_cool, AH_vent, fill=0)
	colnames(all)<-c("AHflag", "None", "Heat", "Cool", "Vent")
	return(all)
}

#colSelct takes a data object and a flexible number of column indexes, and returns a new data object including just those columns selected by the user. The dygraph plotting utility requires an entire data object, rather than selecting columns from existing object.

colSelect=function(dat_obj, col1, ...){
	return(dat_obj[,c(col1, ...)])
}


rcoCondense=function(filename, time_period="hours"){
	
	#Create date-time index
	date_mins<-seq(ymd_hm("2015-01-01 00:00", tz='GMT'), ymd_hm("2015-12-31 23:59", tz='GMT'), by="min")

	#Read in file as data.table object
	dat_rco<-fread(filename, sep="\t" , header=TRUE)

	#Convert to xts time series object
	dat_rco.xts<-xts(dat_rco, order.by=date_mins)

	temps<-4:9
	for(i in temps){
		dat_rco.xts[,i]<-dat_rco.xts[,i]-273.15
	}

	#Replace any end-of-cycle 102 values with the prior index value
	dat_rco.xts$AHflag<-endCycleReassign(dat_rco.xts$AHflag)
	
	#Splits the AHflag column into separate columns for each AH state (Heat, Cool, Vent, None)
	dat_rco.xts<-merge(dat_rco.xts, AHUflagIndex(dat_rco.xts$AHflag))
	
	#Endpoint index for hours of the year; can also do “us” (microseconds), “microseconds”, “ms” (milliseconds), “milliseconds”, “secs” (seconds), “seconds”, “mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
	time_period<-endpoints(dat_rco.xts, time_period)
	
	#Calculate the mean hourly value for each column of the file
	dat_hourly<-period.apply(dat_rco.xts, INDEX=time_period, FUN=mean)
	
	return(dat_hourly[,3:ncol(dat_hourly)])
}

files<-list.files()[grep(".rco", list.files(), fixed=TRUE)]

for(i in files[1:length(files)]){
	for(j in 2:length(args)){
		sim_name<-substr(i, 1 ,(nchar(i)-4))
		extname<-paste(".", args[j], sep="")
		textname<-paste(sim_name, extname, sep="")
		dat<-rcoCondense(i, args[j]) #args[1] should be set to "hours"
		write.zoo(dat, textname, sep=",")
	}
}


















# # 
# #Create date-time index
# date_mins<-seq(ymd_hm("2015-01-01 00:00"), ymd_hm("2015-12-31 23:59"), by="min")

# #Read in file as data.table object
# dat_rco<-fread("test01_11.rco", sep="\t" , header=TRUE)

# #Convert to xts time series object
# dat_rco.xts<-xts(dat_rco, order.by=date_mins)

# temps<-4:9
# for(i in temps){
	# dat_rco.xts[,i]<-dat_rco.xts[,i]-273.15
# }

# #Replace any end-of-cycle 102 values with the prior index value
# dat_rco.xts$AHflag<-endCycleReassign(dat_rco.xts$AHflag)

# #Splits the AHflag column into separate columns for each AH state (Heat, Cool, Vent, None)
# dat_rco.xts<-merge(dat_rco.xts, AHUflagIndex(dat_rco.xts$AHflag))

# #Endpoint index for hours of the year; can also do “us” (microseconds), “microseconds”, “ms” (milliseconds), “milliseconds”, “secs” (seconds), “seconds”, “mins” (minutes), “minutes”, “hours”, “days”, “weeks”, “months”, “quarters”, and “years”
# hours<-endpoints(dat_rco.xts, "hours")

# #Calculate the mean hourly value for each column of the file
# dat_hourly<-period.apply(dat_rco.xts, INDEX=hours, FUN=mean)

 # # [1] "Time"             "Min"              "windSpeed"        "tempOut"         
 # # [5] "tempHouse"        "setpoint"         "tempAttic"        "tempSupply"      
 # # [9] "tempReturn"       "AHflag"           "AHpower"          "Hcap"            
# # [13] "compressPower"    "Ccap"             "mechVentPower"    "HR"              
# # [17] "SHR"              "Mcoil"            "housePress"       "Qhouse"          
# # [21] "ACH"              "ACHflue"          "ventSum"          "nonRivecVentSum" 
# # [25] "fan1"             "fan2"             "fan3"             "fan4"            
# # [29] "fan5"             "fan6"             "fan7"             "rivecOn"         
# # [33] "turnover"         "relExpRIVEC"      "relDoseRIVEC"     "occupiedExpReal" 
# # [37] "occupiedDoseReal" "occupied"         "occupiedExp"      "occupiedDose"    
# # [41] "DAventLoad"       "MAventLoad"       "HROUT"            "HRhouse"         
# # [45] "RH.house"         "RHind60"          "RHind70"          "AHflag.1"        
# # [49] "None"             "Heat"             "Cool"             "Vent"            


# sub_dat<-colSelect(dat_hourly,4,5,6,7,8,9)

# dygraph(sub_dat) %>%
	# dyRangeSelector() %>%
	# dyRoller(rollPeriod=24)

