
#Read in schedule file with a row for each hour of the day and a column for each of 5 fans
#Dryer, kitchen exhaust, bathOne, bathTwo and bathThree (cols 1:5).

setwd('/Users/brennanless/GoogleDrive/regcap/schedules/')

dat <- read.table('/Users/brennanless/GoogleDrive/SmartVentilation_OccupancyControl/TEST_Fan.txt')

#schedule matrix with a row for each minute of the year, and a column for each fan. 

dryerFan<-matrix(NA, ncol=ncol(dat), nrow=525600)

for(fan in 1:ncol(dat)){ #Cycle through each fan.
  
  min_hour <- 0
  hour_day <- 1
  dryer_Incr<-0
  dryer_Target<-0
  dryer_Rand <- 0
  
  for(min in 1:525600){ #Cycle through each minute of the year.
    #Minute and hour counter accounting.
    min_hour<-min_hour + 1
    if(min_hour == 60){
      min_hour <- 0
      hour_day <- hour_day + 1
    }
    if(hour_day == 25){
      hour_day <- 1
    }
    
    #dryerFan[min, 6]<-min_hour
    #dryerFan[min, 7]<-hour_day
  
    #Fan operation. If an hour has fan operation, then at the first minute of that hour a random
    #integer is assigned. When the min_hour = that random integer, the fan turns on, and it turns off
    #When its counter reaches the specficied number of minutes of operation (Target).
    
    if(min_hour == 0 & dat[hour_day, fan] != 0){
      dryer_Rand <- sample(0:59, 1, replace=TRUE)
    }
    if(min_hour == dryer_Rand){
      dryerFan[min, fan] <- 1
      dryer_Incr <- dryer_Incr + 1
      dryer_Target <- dat[hour_day, fan]
    }
    if(dryer_Incr != 0 & dryer_Incr <= dryer_Target){
      dryerFan[min, fan] <- 1
      dryer_Incr <- dryer_Incr + 1
    } else{
      dryerFan[min, fan] <- 0
      dryer_Incr <- 0
      dryer_Target <- 0
    }
  }
}

write.table(dryerFan, 'REGCAP_FanSchedule.txt', quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
