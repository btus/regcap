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
#include "gauss.h"
#include "constants.h"

using namespace std;

int main() {
double woodThick = 1;
double atticVolume = 100;
double atticArea = 10;
double roofPitch = 14;
double tempInit = airTempRef;
double mcInit = 0.15;

	// initialize attic moisture nodes
	WoodMoisture attic(woodThick, atticVolume, atticArea, roofPitch, tempInit, mcInit);
	
	double tempOut = airTempRef-10;
	double RHOut = 0.5;
	double tempHouse = airTempRef;
	double RHHouse = 0.4;
	double airDensityAttic = airDensityRef * airTempRef / airTempRef;
	double airDensityHouse = airDensityRef * airTempRef / tempHouse;
	double pressure;
	double hU0, hU1, hU2;
   double mAttic, mCeiling;
   
	mass_cond_bal(tempOut, RHOut, tempHouse, RHHouse, airDensityOut, airDensityAttic, airDensityHouse, pressure,
	              hU0, hU1, hU2, mAttic,  mCeiling);
}

