/* Unit tests for functions in functions.cpp
	Leo Rainer (4/6/16)
*/
#include <sstream>
#include <iostream>
#include <fstream>
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include "functions.h"
#include "constants.h"

using namespace std;

void f_roofCpTheta(double* Cproof, int& windAngle, double* Cppitch, double& roofPitch);
void test_roofCp(double roofPitch, bool rowHouse, bool ridgeParallel);

int main() {
test_roofCp(5.0, false, true);
test_roofCp(20.0, false, true);
test_roofCp(35.0, false, true);
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