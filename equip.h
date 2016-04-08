#pragma once
#ifndef equip_h
#define equip_h

//using namespace std;

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
		
	public:
		double SHR;
		double capacitySensible;	// sensible capacity at current conditions (watts)
		double capacityLatent;		// latent capacity at current conditions (watts)
		double condensate;			// condensed water removed (kg)

		Compressor(double EER, double tons, double charge=1.0);
		double run(bool compOn, double hrReturn, double tReturn, double tOut, double fanFlow, double fanHeat, double mAH);
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

#endif