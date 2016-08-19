/* Unit tests for wood moisture function
	Leo Rainer (7/22/16)
*/
#include <sstream>
#include <iostream>
#include <fstream>
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include "wood.h"
//#include "gauss.h"
#include "constants.h"

using namespace std;

// to compile: g++ test_wood.cpp wood.cpp psychro.cpp gauss.cpp -o test_wood

int main() {
	double woodThick = 0.015;   // should match what is set in sub_heat for now
	double atticVolume = 300;
	double atticArea = 300;
	double roofPitch = 18;
	double tempInit = airTempRef;
	double mcInit = 0.15;

	// initialize attic moisture nodes
	WoodMoisture attic(woodThick, atticVolume, atticArea, roofPitch, tempInit, mcInit);
	
	double tempOut = airTempRef - 10;
	double RHOut = 50.0;
	double tempHouse = airTempRef;
	double RHHouse = 40.0;
	double airDensityOUT;
	double airDensityATTIC;
	double airDensityIN;
	int pressure = 101325;
	double H2, H4, H6;
   double mAttic = 1.0;
   double mCeiling = 0.5;
   double b[16];
   
   // initial temps
   b[0] = 293;	// attic
   b[1] = 293;	// inner north
   b[2] = 293;	// outer north
   b[3] = 293;	// inner south
   b[4] = 293;	// outer south
   b[5] = 293;	// bulk wood
   
   // mAttic = matticenvin - matticenvout;
      
   for(int m = 0; m < 2; m++) {
		// set node temperatures 
		attic.temperature[0] = b[3];
		attic.temperature[1] = b[1];
		attic.temperature[2] = b[5];
		attic.temperature[3] = b[0];
		attic.temperature[4] = (b[3] + b[4]) / 2;
		attic.temperature[5] = (b[1] + b[2]) / 2;
		attic.temperature[6] = b[5];
	
		// set heat transfer coefficients to natural value
		H4 = 3.2 * pow(abs(b[3] - b[0]), 1.0/3.0);	// inner south
		H2 = 3.2 * pow(abs(b[1] - b[0]), 1.0/3.0);	// inner north
		H6 = 3.2 * pow(abs(b[5] - b[0]), 1.0/3.0);	// bulk wood

		// set air densities
		airDensityOUT = airDensityRef * airTempRef / tempOut;
		airDensityATTIC = airDensityRef * airTempRef / b[0];
		airDensityIN = airDensityRef * airTempRef / tempHouse;
	
		attic.mass_cond_bal(tempOut, RHOut, tempHouse, RHHouse, airDensityOUT, airDensityATTIC, airDensityIN, pressure,
						  H4, H2, H6, mAttic,  mCeiling);
	
		for(int i = 0; i < 7; i++) {
			cout << m << ":" << i << "," << attic.temperature[i] << "," << attic.moistureContent[i] << "," << attic.PW[i] << "," << attic.mTotal[i] << endl;
			}
		}
}

