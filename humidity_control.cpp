// Humidity Control Logic
// This class has NOT been compiled or tested - just moved here to get it out of main

humidity_control_constructor()
// The following need to be initialized:
		//double wCutoff;			// Humidity Ratio cut-off calculated as some percentile value for the climate zone. Brennan.
		//double wDiffMaxNeg;		// Maximum average indoor-outdoor humidity differene, when wIn < wOut. Climate zone average.
		//double wDiffMaxPos;		// Maximum average indoor-outdoor humidity differene, when wIn > wOut. Climate zone average.
		//double W25[12]; 			// 25th percentiles for each month of the year, per TMY3
		//double W75[12]; 			// 75th percentiles for each month of the year, per TMY3
		//int FirstCut;				// Monthly Indexes assigned based on climate zone
		//int SecondCut; 			// Monthly Indexes assigned based on climate zone
		//double doseTarget;		// Targeted dose value
		//double HiDose;				// Variable high dose value for real-time humidity control, based on worst-case large, low-occupancy home. 
		//int HiMonths[3];
		//int LowMonths[3];
		//double HiMonthDose;
		//double LowMonthDose;

// This was in the hour loop. Better to do the hcFlag decision once per day
				if (hour == peakEnd)
					peakFlag = 1;			// Prevents two peak periods in the same day when there is heating and cooling


int humidity_control(double HROut, double HRHouse, double RHhouse, double relExp, double relDose, int hour, int month, int AHflag, double& occupied) {
	int rivecOn;
	double relExpTarget;		//This is a relative expsoure value target. It varies between 0 and 2.5, depending on magnitude of indoor-outdoor humidity ratio difference.

//	Cooling system tie-in.
	if(HumContType == 1){			   				
		if(hcFlag == 2){ //test if we're in cooling season
			if(AHflag == 2){
				rivecOn = 1;
			} else
				if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;					
				}
		} else //if NOT in cooling season
			if(relExp >= 0.95 || relDose > 1.0){ //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;					
				}
		}

//	Fixed control. 	   				
//	Indoor and Outdoor sensor based control. 
	if(HumContType == 2){
		if(RHhouse >= 55){ //Engage increased or decreased ventilation only if house RH is >60 (or 55% maybe?). 
			if(HROut > HRHouse){ //do not want to vent. Add some "by what amount" deadband value. 
				if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relExp >= 0.50 || relDose > 1.0){ //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(relExp >= 0.95 || relDose > 1.0){ 
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Fixed control + cooling system tie-in.		   				
//	Indoor and Outdoor sensor based control. 
	if(HumContType == 3){							
		if(RHhouse >= 55){ //Engage increased or decreased ventilation only if house RH is >60 (or 55% maybe?). 
			if(HROut > HRHouse){ //do not want to vent. Add some "by what amount" deadband value. 
				if(AHflag == 2){
					rivecOn = 1;
				} else
					if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
						rivecOn = 1;
					} else { //otherwise off
						rivecOn = 0;
					}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relExp >= 0.50 || relDose > 1.0){ //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(relExp >= 0.95 || relDose > 1.0){ 
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Proportional control.  				   				
//	Indoor and Outdoor sensor based control. 
	if(HumContType == 4) {
		if(RHhouse >= 55) {
			if(HROut > HRHouse){ //More humid outside than inside, want to under-vent.  
				relExpTarget = 1 + (2.5-1) * abs((HRHouse-HROut) / (wDiffMaxNeg)); //wDiffMax has to be an avergaed value, because in a real-world controller you would not know this. 
				if(relExpTarget > 2.5){
					relExpTarget = 2.5;
				}
				if(relExp >= relExpTarget || relDose > 1){ //relDose may be over a 1-week or 2-week time span...
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { // More humid inside than outside, want to over-vent
				relExpTarget = 1 - abs((HRHouse- HROut) / (wDiffMaxPos));
				if(relExpTarget < 0){
					relExpTarget = 0;
				}
				if(relExp >= relExpTarget || relDose > 1){ //
					rivecOn = 1;
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(relExp >= 0.95 || relDose > 1){ 
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Proportional control + cooling system tie-in. 				   				
//	Indoor and Outdoor sensor based control.
	if(HumContType == 5) { 
		if(RHhouse >= 55) {
			if(HROut > HRHouse) { //More humid outside than inside, want to under-vent.  
				relExpTarget = 1 + (2.5-1) * abs((HRHouse-HROut) / (wDiffMaxNeg)); //wDiffMax has to be an avergaed value, because in a real-world controller you would not know this. 
				if(relExpTarget > 2.5) {
					relExpTarget = 2.5;
				}
				if(AHflag == 2) {
					rivecOn = 1;
				} else if(relExp >= relExpTarget || relDose > 1) { //relDose may be over a 1-week or 2-week time span...
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { // More humid inside than outside, want to over-vent
				relExpTarget = 1 - abs((HRHouse- HROut) / (wDiffMaxPos));
				if(relExpTarget < 0) {
					relExpTarget = 0;
				}
				if(relExp >= relExpTarget || relDose > 1) {
					rivecOn = 1;
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(relExp >= 0.95 || relDose > 1) { 
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Monthly Seasonal Control
//	Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
	if(HumContType == 6) {
		if(month <= FirstCut || month >= SecondCut) { //High ventilation months with net-humidity transport from inside to outside.
			if(relDose > doseTarget) {
				rivecOn = 1;
			} else {						
				rivecOn = 0;
			} 
		} else { //Low ventilation months with net-humidity transport from outside to inside.
			if(relExp >= 2.5 || relDose > HiDose) { //Brennan changed from fixed 1.5 to variable HiDose value.
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Fixed control + cooling system tie-in + Monthly Seasonal Control.		   				
//	Indoor and Outdoor sensor based control.
	if(HumContType == 7) { 
		if(HROut > HRHouse) { //do not want to vent. 
			if(AHflag == 2) {
				rivecOn = 1;
			} else if(relExp >= 2.5 || relDose > HiDose) { 
				rivecOn = 1;
			} else { //otherwise off
				rivecOn = 0;
			}
		} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
			if(relDose > doseTarget) { 
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		}
	}

//	Fixed control + cooling system tie-in + Monthly Seasonal Control.		   				
//	Indoor and Outdoor sensor based control.
	if(HumContType == 7) { 
		if(HROut > HRHouse) { //do not want to vent. 
			if(RHhouse >= 55){
			if(AHflag == 2) {
				rivecOn = 1;
			} else if(relExp >= 2.5 || relDose > HiDose) { 
				rivecOn = 1;
			} else { //otherwise off
				rivecOn = 0;
			}
			}
			else if(relExp >= 0.95 || relDose > 1){ 
				rivecOn = 1;
			} 
			else{
				rivecOn = 0;
			}		
		}		
		else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
			if(relDose > doseTarget) { 
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		}
	}

//	Monthly Seasonal controller + time of day
//	Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average
//	targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
	if(HumContType == 8) {
		if(month <= FirstCut || month >= SecondCut) { //High ventilation months with net-humidity transport from inside to outside.
			if(hour >= 3 && hour <= 7) { //Was 3 and 7. Maybe change this to the warmest hours of the day.
				rivecOn = 1; //Relatively dry time of day, vent more.
			} else {
				if(relDose > doseTarget) {
					rivecOn = 1;
				} else {						
					rivecOn = 0;
				} 
			}
		} else { //Low ventilation months with net-humidity transport from outside to inside.
			if(hour >= 14 && hour <= 18) { //Brennan changed from 12 and 6, to 14 to 18. Always vent during peak cooling period.
				rivecOn = 1;
			} else if(hour >= 7 && hour <= 11){
				rivecOn = 0;
			} else {
				if(relExp >= 2.5 || relDose > HiDose) {
					rivecOn = 1;
				} else {
					rivecOn = 0;
				}
			}
		}
	}

//	Fixed control + Monthly Seasonal Control.		   				
//	Indoor and Outdoor sensor based control.
	if(HumContType == 9) { 
		if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
			if(relExp >= 2.5 || relDose > HiDose) { //have to with high exp 0.61
				rivecOn = 1;
			} else { //otherwise off
				rivecOn = 0;
			}
		} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
			if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
				rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
			} else {
				rivecOn = 0;
			}
		}
	}

//	Fixed control + Monthly Seasonal Control.		   				
//	Indoor and Outdoor sensor based control.
	if(HumContType == 9) { 
		if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
			if(RHhouse >= 55) {
			if(relExp >= 2.5 || relDose > HiDose) { //have to with high exp 0.61
				rivecOn = 1;
			} 
			else { //otherwise off
				rivecOn = 0;
			}
			}	
			else if(relExp >= 0.95 || relDose > 1){ 
				rivecOn = 1;
			} 
			else{
				rivecOn = 0;
			}	
		}			
		else if(relDose > doseTarget){ //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
			rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
		} 
		else {
			rivecOn = 0;
		}
	}						


//	The real opportunities for control based on outside are when the outside value is changing rapidly. Sharp increases, decrease ventilation. Sharp decreases, increase ventilation. 
//	Need to undervent at above the 75th percentile and overvent below the 25th percentile based on a per month basis. 
//	Outdoor-only sensor based control
	if(HumContType == 10) {
		if(HROut > 0.012) { //If humid outside, reduce ventilation. 
			if(relExp >= 2.5 || relDose > 1) {
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		} else { //If dry outside, increase vnetilation.
			if(relExp >= 0.95 || relDose > 1) { //but we do if exp is high
				rivecOn = 1;
			} else {
				rivecOn = 0; //otherwise don't vent under high humidity condition.
			}
		}
	}

//	Monthly Advanced Seasonal Control
//	Monthly timer-based control, based on mean HRdiff by month. 
//	doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
//	Setting the appropriate Dose Target based on the month						
	if(HumContType == 11) {
		double doseTargetTmp;
		if(month == HiMonths[0] || month == HiMonths[1] || month == HiMonths[2]){
			doseTargetTmp = HiMonthDose;
		} else if(month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
			doseTargetTmp = LowMonthDose;
		} else if(month > FirstCut && month < SecondCut) {
			doseTargetTmp = HiDose; //Brennan changed from fixed 1.5 to HiDose.
		} else {
			doseTargetTmp = doseTarget;
		}
		if(doseTargetTmp >= HiDose) {
			if(relExp >= 2.5 || relDose > doseTargetTmp) {
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		}
		else {
			if(relDose > doseTargetTmp) {
				rivecOn = 1;
			} else {						
				rivecOn = 0;
			} 
		}
	}

//	Consider combining 8 and 12 - Monthly + Time of day + Cooling tie-in.
//	Monthly Seasonal Control + Cooling system tie-in.		   				
//	Indoor and Outdoor sensor based control. 
	if(HumContType == 12) {
		double doseTargetTmp;
		if(month <= FirstCut || month >= SecondCut) {
			doseTargetTmp = doseTarget;
		} else {
			doseTargetTmp = HiDose;
		}
		if(AHflag == 2) {
				rivecOn = 1;
		} else if(doseTargetTmp >= HiDose) {
			if(relExp >= 2.5 || relDose > doseTargetTmp) {
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		} else {
			if(relDose > doseTargetTmp) {
				rivecOn = 1;
			} else {						
				rivecOn = 0;
			} 
		}
	}

//	Outdoor-only sensor based control, with variable dose targets
	if(HumContType == 13) { 
		if(HROut >= wCutoff) { //If humid outside, reduce ventilation. 
			if(relExp >= 2.5 || relDose > 1.45) {
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		} else { //If dry outside, increase ventilation.
			if(relDose > 0.45) { //but we do if exp is high
				rivecOn = 1;
			} else {
				rivecOn = 0; //otherwise don't vent under high humidity condition.
			}
		}
	}

//	Outdoor-only sensor based control, control based on 25th and 75th percentile monthly values for each month and climate zone.
	if(HumContType == 14) { 
		if(HROut >= W75[month-1]) { //If humid outside, reduce ventilation. 
			if(relExp >= 2.5 || relDose > 1.5) {
				rivecOn = 1; 
			} else {
				rivecOn = 0;
			}
		} else if (HROut <= W25[month-1]) { //If dry outside, increase vnetilation.
			if(relDose > 0.5) { //but we do if exp is high
				rivecOn = 1;
			} else {
				rivecOn = 0; //otherwise don't vent under high humidity condition.
			}
		} else {
			if(relExp >= 0.95 || relDose > 1.0){ //Need to reduce this target to make equivalence work out...
				rivecOn = 1;
			} else {
				rivecOn = 0;
			}
		}
	}

//	Monthly Advanced Seasonal Control + Fixed Control + Cooling System Tie_in
//	Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
//	Setting the appropriate Dose Target based on the month
	if(HumContType == 15) {
		if(month == HiMonths[0] || month == HiMonths[1]  || month == HiMonths[2] ||
			month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
			doseTargetTmp = HiMonthDose;
			if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
				if(AHflag == 2) {
					rivecOn = 1;
				} else if(relExp >= 2.5 || relDose > LowMonthDose) { //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relDose > HiMonthDose) { //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
				if(AHflag == 2) {
					rivecOn = 1;
				} else if(relExp >= 2.5 || relDose > 1.5) { //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		}
	}

//	Monthly Advanced Seasonal Control + Fixed Control
//	Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
//	Setting the appropriate Dose Target based on the month
	if(HumContType == 16) {
		if(month == HiMonths[0] || month == HiMonths[1]  || month == HiMonths[2] || 
			month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
			doseTargetTmp = HiMonthDose;
			if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
				if(relExp >= 2.5 || relDose > LowMonthDose) { //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relDose > HiMonthDose) { //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		} else {
			if(HROut > HRHouse) { //do not want to vent. Add some "by what amount" deadband value. 
				if(relExp >= 2.5 || relDose > 1.5) { //have to with high exp
					rivecOn = 1;
				} else { //otherwise off
					rivecOn = 0;
				}
			} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
				if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
					rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
				} else {
					rivecOn = 0;
				}
			}
		}
	}

//	Rivec control of number 50 fan type, algorithm v6					     
	if(HumContType == 0) {
		if(occupied[weekend][hour]) {				            	// Base occupied
			if(relExp >= 0.95 || relDose >= 1.0)
				rivecOn = 1;
		} else {						                	// Base unoccupied
			if(relExp >= expLimit)
				rivecOn = 1;
		}
		if(hour >= peakStart && hour < peakEnd && peakFlag == 0) {		// PEAK Time Period
			rivecOn = 0;												// Always off
			if(relExp >= expLimit)
				rivecOn = 1;
		}
	}

	return rivecOn
}
				
