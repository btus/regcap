#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>		// RAD: so far used only for setprecission() in cmd output
#include <time.h>
#include <vector>
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include "psychro.h"
#include "constants.h"
#include "equip.h"

using namespace std;

int main() {

	Compressor centralAC(12.0, 3.0, 1.0);	// 12 SEER, 3 ton

	double hrReturn = 0.008;
	double tReturn = 273 + 20;
	double tOut = 273 + 30;
	double fanFlow = 400 * 3 * .0004719;
	double fanHeat = 500;
	double mAH = fanFlow * airDensityRef;

	for(int i = 0; i<20; i++) {
		float power = centralAC.run(1, hrReturn, tReturn, tOut, fanFlow, fanHeat, mAH);
		cout << "Power=" << power << " Sensible=" << centralAC.capacitySensible << " Latent=" << centralAC.capacityLatent << " condensate=" << centralAC.condensate << endl;
	}
	for(int i = 0; i<30; i++) {
		float power = centralAC.run(0, hrReturn, tReturn, tOut, fanFlow, fanHeat, mAH);
		cout << "Power=" << power << " Sensible=" << centralAC.capacitySensible << " Latent=" << centralAC.capacityLatent << " condensate=" << centralAC.condensate << endl;
	}
}