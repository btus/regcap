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

using namespace std;

void test_roofCp(double roofPitch, bool rowHouse, bool ridgeParallel);

int main() {
double woodThick;
double atticVolume;
double atticArea;
double roofPitch;
double tempInit;
double mcInit;

	WoodMoisture attic(woodThick, atticVolume, atticArea, roofPitch, tempInit, mcInit);
	
	
	mass_cond_bal(double tempOut, double RHOut, double tempHouse, double RHHouse,
                             double airDensityOut, double airDensityAttic, double airDensityHouse,
                             double pressure, double hU0, double hU1, double hU2,
                             double mAttic, double mCeiling);
}

void test_roofCp(double roofPitch, bool rowHouse, bool ridgeParallel) {
	int windAngle;
	double Cproof[4];
	double Cppitch[4];

	cout << "Roof Pitch =" << roofPitch << endl;
	for(windAngle = 0; windAngle < 360; windAngle += 5) {
		// from sub_atticLeak
		if(roofPitch < 10) {
			Cproof[0] = -.8;
			Cproof[1] = -.4;
		} else if(roofPitch > 30) {
			Cproof[0] = .3;
			Cproof[1] = -.5;
		} else {
			Cproof[0] = -.4;
			Cproof[1] = -.4;
		}
		if(rowHouse) {
			Cproof[2] = -.2;
			Cproof[3] = -.2;
		} else {
			// for isolated houses
			Cproof[2] = -.6;
			Cproof[3] = -.6;
		}
		if(!ridgeParallel) {
			swap(Cproof[2],Cproof[0]);
			swap(Cproof[3],Cproof[1]);
		}

		f_roofCpTheta(Cproof, windAngle, Cppitch, roofPitch);
		cout << windAngle << "," << Cppitch[0];
		//for(int i = 0; i < 4; i++) {
		//	cout << "," << Cppitch[i];
		//}
		cout << endl;
	}
}