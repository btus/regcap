#include <iostream>
#include "equip.h"
#include "psychro.h"
#include "constants.h"
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

Compressor::Compressor(double EER, double tons, double charge) {
	// The following corrections are for TXV only
	chargeCapDry = (charge < .725) ? 1.2 + (charge - 1) : .925;
	chargeEffDry = (charge < 1) ? 1.04 + (charge - 1) * .65 : 1.04 - (charge - 1) * .35;
	chargeCapWet = (charge < .85) ? 1 + (charge - .85) : 1;
	if(charge < .85) {
		chargeEffWet = 1 + (charge - .85) * .9;
	} else if(charge > 1) {
		chargeEffWet = 1 - (charge - 1) * .35;
	} else {
		chargeEffWet = 1;
	}

	capacityRated = tons;
	efficiencyRated = EER;
	maxMoisture = 0.3 * tons;		// maximum mass on coil is 0.3 kg per ton of cooling
	onTime = 0;
	coilMoisture = 0;
}							
							
double Compressor::run(bool compOn, double hrReturn, double tReturn, double tOut, double fanFlow, double fanHeat, double mAH) {
	const int rampTime = 3;	// time to reach full SHR (minutes)
	double hReturn;			// return air enthalpy (btu/lb)
	double capacityTotal;	// total capacity at current conditions (btuh)
	double EER;					// efficiency at current conditions (kbtuh/kWh)
	double compressorPower;	// compressor power (Watts)

	if(compOn) {	// compressor running
		onTime++;
		tOut = KtoF(tOut);			// cap and eer equations use Tout in deg F
		tReturn = KtoF(tReturn);	// enthalpy calc uses temps in deg F
	
		hReturn = 0.24 * tReturn + hrReturn * (1061 + 0.444 * tReturn);
		hReturn += fanHeat / mAH / 2326;	// add fan heat (in btu/lb)
			
		// Cooling capacity air flow correction term (in CFM)
		double fanCFM = fanFlow / .0004719;
		double fanCorr = 1.62 - .62 * fanCFM / (400 * capacityRated) + .647 * log(fanCFM / (400 * capacityRated));

		// Set sensible heat ratio using humidity ratio - note only because we have small indoor dry bulb range			
		SHR = 1 - 50 * (hrReturn - .005);
		if(SHR < .25)
			SHR = .25;
		if(SHR > 1)
			SHR = 1;
		// SHR ramps up for first rampTime minutes of operation
		if(onTime < rampTime) {
			SHR += (rampTime - onTime) / rampTime * (1 - SHR);
		}

		if(SHR == 1) { // dry coil
			capacityTotal = (capacityRated * 12000) * .91 * fanCorr * chargeCapDry * (-.00007 * pow((tOut - 82),2) - .0067 * (tOut - 82) + 1);
			EER = efficiencyRated * fanCorr * chargeEffDry * ((-.00007) * pow((tOut - 82),2) - .0085 * (tOut - 82) + 1);
		} else {			// wet coil
			capacityTotal = capacityRated * 12000 * (1 + (hReturn - 30) * .025) * fanCorr * chargeCapWet * (-.00007 * pow((tOut - 95),2) - .0067 * (tOut - 95) + 1);
			EER = efficiencyRated * fanCorr * chargeEffWet * ((-.00007) * pow((tOut - 95),2) - .0085 * (tOut - 95) + 1);
		}
		compressorPower = capacityTotal / EER;
		capacitySensible = SHR * capacityTotal / 3.413 - fanHeat;
	} else {
		onTime = 0;					// reset on time counter
		SHR = 1;						// so latent capacity with fan on is calculated correctly
		compressorPower = 0;
		capacitySensible = 0;
	}
	
	// Track Coil Moisture
	if(SHR < 1) {
		capacityLatent = (1 - SHR) * capacityTotal / 3.413;			// latent capacity Watts
		coilMoisture += capacityLatent / 2501000 * dtau;				// condenstion in timestep kg/s * time
	} else {
		if(coilMoisture > 0) {													// if moisture on coil but no latcap, then we evaporate coil moisture until Mcoil = zero
			capacityLatent = -maxMoisture / 1800 * 2501000;				// evaporation capacity J/s - this is negative latent capacity
			coilMoisture -=  maxMoisture / 1800 * dtau;					// evaporation in timestep kg/s * time
		} else {
			capacityLatent = 0;
		}
	}

	if(coilMoisture < 0)
		coilMoisture = 0;
	if(coilMoisture > maxMoisture) {
		condensate = coilMoisture - maxMoisture;
		coilMoisture = maxMoisture;
	}
	
return compressorPower;
}

Dehumidifier::Dehumidifier(double capacity, double energyFactor, double rh) {
	capacityRated = capacity * 0.4732 / (24 * 60 * 60);	// convert pints/day to kg/s
	efficiency = 1000 / energyFactor;							// convert to Wh/kg
	setPoint = rh;
	onTime = 0;
}

bool Dehumidifier::run(double rhIn, double tIn) {
	const double deadBand = 0.05;	// 5% RH deadband
	const int initTime = 3;			// time to full capacity (minutes)
	// Curves from Winkler et. al., NREl Technical Report TP-5500-52791, December 2011
	// "Laboratory Test Report for Six ENERGY STAR® Dehumidifiers"
	// Capacity curve
	const double capA = -1.1625;
	const double capB = 0.022715;
	const double capC = -0.00011321;
	const double capD = 0.021111;
	const double capE = -6.9303E-05;
	const double capF = 0.00037884;
	// Energy Factor curve
	const double efA = -1.9022;
	const double efB = 0.063467;
	const double efC = -0.00062284;
	const double efD = 0.039540;
	const double efE = -0.00012564;
	const double efF = -0.00017672;
	
	
	// Control
	if(onTime > 0) {
		if(rhIn < setPoint - deadBand) {
			onTime = 0;
		} else {
			onTime++;
		}
	} else {
		if(rhIn > setPoint + deadBand) {
			onTime = 1;
		}
	}
	
	// Operate
	if(onTime > 0) {
		tIn -= C_TO_K;	// curves are in deg C and %RH
		rhIn *= 100;
		double capFTRH = capA + capB * tIn + capC * pow(tIn, 2) + capD * rhIn + capE * pow(rhIn, 2) + capF * tIn * rhIn;
		moisture = capacityRated * capFTRH; 
		double efFTRH = efA + efB * tIn + efC * pow(tIn, 2) + efD * rhIn + efE * pow(rhIn, 2) + efF * tIn * rhIn;
		power = moisture * 3600 * efficiency / efFTRH;		// kg/s * 3600 s/h * Wh/kg = Watts
		// adjust moisture for init time
		double capInit = (onTime < initTime) ? 1.0 - double(initTime - onTime) / double(initTime) : 1.0;
		moisture *= capInit;
		sensible = moisture * 2501000 + power;
//cout << " onTime=" << onTime << " capFTRH=" << capFTRH << " capInit=" << capInit << " efFTRH=" << efFTRH << endl;
	} else {
		// implement offTime here if we want to keep fan running
		moisture = 0;
		power = 0;
		sensible = 0;
	}
	return (onTime > 0);
}
	

