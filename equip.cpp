//#include "equip.h"
//#include "psychro.h"
//#include "constants.h"
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

class Compressor {
	private:
		double chargeCapDry;		// charge correction for capacity (dry coil)
		double chargeCapWet;		// charge correction for capacity (wet coil)
		double chargeEffDry;		// charge correction for EER (dry coil)
		double chargeEffWet;		// charge correction for EER (wet coil)
		double capacityRated;	// capacity at rated conditions (tons)
		double efficiencyRated;	// efficiency at rated conditions (kbtuh/kWh)
		double onTime;				// tracks the number of minutes the compressor is on for
		double coilMoisture;		// water on coil (kg)
		double maxMoisture;		// maximum amount of water that can be on coil (kg)
		const double rampTime = 3;	// time to reach full SHR (minutes)
		
	public:
		double SHR;
		double capacitySensible;	// sensible capacity at current conditions (watts)
		double capacityLatent;		// latent capacity at current conditions (watts)
		double condensate;			// condensed water removed (kg)

		Compressor(double EER, double tons, double charge=1.0) {
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
							
		double run(bool compOn, double hrReturn, double tReturn, double tOut, double fanFlow, double fanHeat, double mAH) {
			double hReturn;			// return air enthalpy (btu/lb)
			double capacityTotal;	// total capacity at current conditions (btuh)
			double EER;					// efficiency at current conditions (kbtuh/kWh)

			if(compOn) {	// compressor running
				onTime++;
				tOut = (tOut - 273.15) * (9.0 / 5.0) + 32;	// cap and eer equations use Tout in deg F
			
				hReturn = .24 * ((tReturn - 273.15) * 9 / 5 + 32) + hrReturn * (1061 + .444 * ((tReturn - 273.15) * 9 / 5 + 32));
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
					evapcap = latcap;
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
			
		return commpressorPower;
		}
};

class Fan {
	private:
		double power;
		double flow;		// was q
		double efficiency;
		
	public:
		double type;		// was oper
		double massFlow;	// was m
		double on; 			// Set on=1 in main program based on .oper
		
		
};