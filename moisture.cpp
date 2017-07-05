#include "moisture.h"
#include "psychro.h"
#include "constants.h"
#include "gauss.h"
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include <iostream>
#include <vector>

using namespace std;

/*
 * Moisture - Moisture class constructor
 * @param atticVolume - attic volume (m3)
 * @param supVolume - supply register volume (m3)
 * @param retVolume - return register volume (m3)
 * @param houseVolume - house volume (m3)
 * @param floorArea - conditioned floor area (m2)
 * @param atticArea - attic area (m2)
 * @param roofPitch - roof pitch (degrees)
 * @param ERV_TRE - ERV total recovery efficiency (fraction)
 * @param mcInit - initial wood moisture content (fraction - optional)
 *
 * The nodes are (0-5 wood, 6-10 air/mass):
 * 0. surface of south sheathing (facing the attic) - corresponding thermal node is (4)
 * 1. surface of north sheathing (facing the attic) - corresponding thermal node is (2)
 * 2. surface of the bulk wood in the attic - corresponding thermal node is (6)
 * 3. inside south sheathing - corresponding thermal node is ((4+5)/2)
 * 4. inside north sheathing - corresponding thermal node is ((2+3)/2)
 * 5. inside bulk wood - corresponding thermal node is (6)
 * 6. attic air - corresponding thermal node is (1)
 * 7. return air - corresponding thermal node is (12)
 * 8. supply air - corresponding thermal node is (15)
 * 9. house air - corresponding thermal node is (16)
 * 10. house materials - corresponding thermal node is (13)
 */
Moisture::Moisture(double atticVolume, double retVolume, double supVolume, double houseVolume, 
						double floorArea, double atticArea, double roofPitch, double mcInit) {
   int pressure = 101325;		// 1 atmosphere
   double sheathThick = 0.015;  // Roof sheathing thickness (m). match what is set in sub_heat for now
   double bulkThick = 0.038;    // Bulk wood thickness (m). this is low, should be .038 for 2x4 construction
   double surfaceThick = 0.005; // Surface wood layer thickness (m)
	double tempInit = airTempRef; // node initialization temperature (deg K)
	double RHInit = 50;	// initial air RH (%)

   A.resize(MOISTURE_NODES, vector<double>(MOISTURE_NODES+1, 0));		// set size of equation vectors to number of nodes (A contains both)
   PW.resize(MOISTURE_NODES, 0);
	deltaX[0] = sheathThick / 2;                // distance between the centers of the surface and inside wood layers. it is half the characteristic thickness of the wood member
	deltaX[1] = deltaX[0];                    // same as for the other sheathing surface
	deltaX[2] = bulkThick / 2;
	deltaX[3] = deltaX[0];
	deltaX[4] = deltaX[1];
	deltaX[5] = deltaX[2];

	area[0] = atticArea / 2 / cos(roofPitch * M_PI / 180);
	area[1] = area[0];
	area[2] = atticArea / 2;						// attic wood surface area. @TODO: sub_heat uses planArea * 1.5
	area[3] = area[0];
	area[4] = area[1];
	area[5] = area[2];

	volume[0] = area[0] * surfaceThick;
	volume[1] = area[1] * surfaceThick;
	volume[2] = area[2] * surfaceThick;  // was .001;
	volume[3] = area[0] * sheathThick - volume[0];
	volume[4] = area[1] * sheathThick - volume[1];
	volume[5] = area[2] * sheathThick - volume[2]; // was: 0.02 * area[2] - volume[2];	// this is close for 2x4 construction and may change a bit if 2x6 (I redid the calculations and for 2x4 this would be 0.017 but I dont think we can really justify that second significant digit!
	volume[6] = atticVolume;
	volume[7] = retVolume;
	volume[8] = supVolume;
	volume[9] = houseVolume;

	haHouse = 1 * 0.622 * .5 * floorArea / 186;	// moisture transport coefficient scales with floor area (kg/s) - 0.622 (dHR/dVP), 186 (area of std house in m2), 0.5 empirical coefficient (kg/s)
	massWHouse = 1 * 0.622 * 60 * floorArea;					// active mass containing moisture in the house (kg) - empirical

/*	volume[10] = 10; // @TODO@ need volume */

	// initialize wood nodes
	for(int i=0; i<6; i++) {
		tempOld[i] = tempInit;
		mTotal[i] = 0;
		moistureContent[i] = mcInit;
		PWOld[i] = calc_vapor_pressure(mcInit, tempInit, pressure);
		}
	// initialize air nodes
	for(int i=6; i<MOISTURE_NODES; i++) {
		tempOld[i] = tempInit;
		PWOld[i] = saturationVaporPressure(tempInit) * RHInit / 100;
		}

}							
							
/*
 * mass_cond_bal - Uses the results of the ventilation and heat transfer models to predict
 *                 moisture transport. Includes effects of surfaces at saturation pressure.
 * @param node_temps - array of node temperatures (deg K)
 * @param tempOut - Outdoor temperature (deg K)
 * @param RHOut - Outdoor relative humidity (%)
 * @param airDensityOut - Outdoor air density (kg/m3)
 * @param airDensityAttic - Attic air density (kg/m3)
 * @param airDensityHouse - Indoor air density (kg/m3)
 * @param airDensitySup - Supply air density (kg/m3)
 * @param airDensityRet - Return air density (kg/m3)
 * @param pressure - atmospheric pressure (Pa)
 * @param hU0 - surface heat transfer coefficient for node 0 (inner south) (J/sm2K) - was H4
 * @param hU1 - surface heat transfer coefficient for node 1 (inner north) (J/sm2K) - was H6
 * @param hU2 - surface heat transfer coefficient for node 2 (bulk wood) (J/sm2K) - was H10
 * @param mAtticIn - air mass flow into attic (kg/s)
 * @param mAtticOut - air mass flow out of attic (kg/s)
 * @param mCeiling - ceiling air mass flow (kg/s)
 * @param mHouseIn - air mass flow into house (kg/s)
 * @param mHouseOut - air mass flow out of house (kg/s)
 * @param mAH - air mass flow through the air handler (kg/s)
 * @param mRetOff - air mass flow through return with AH off (kg/s)
 * @param mRetLeak - air mass flow through return leaks (kg/s)
 * @param mRetReg - air mass flow through return register (kg/s)
 * @param mRetOut - air mass flow into return from outside due to HRV/ERV/Fan Cycler (kg/s)
 * @param mErvHouse - ERV mass flow rate from house (kg/s)
 * @param mSupOff - air mass flow through supply with AH off (kg/s)
 * @param mSupLeak - air mass flow through supply leaks (kg/s)
 * @param mSupReg - air mass flow through supply register (kg/s)
 * @param latcap - moisture removed by the AC coil (kg/s)
 * @param dhMoistRemv - moisture removed by the dehumidifier (kg/s)
 * @param latload - indoor latent load (kg/s)
 */
void Moisture::mass_cond_bal(double* node_temps, double tempOut, double RHOut,
                  double airDensityOut, double airDensityAttic, double airDensityHouse, double airDensitySup, double airDensityRet,
                  int pressure, double hU0, double hU1, double hU2,
                  double mAtticIn, double mAtticOut, double mCeiling, double mHouseIn, double mHouseOut,
                  double mAH, double mRetAHoff, double mRetLeak, double mRetReg, double mRetOut, double mErvHouse,
                  double mSupAHoff, double mSupLeak, double mSupReg, double latcap, double dhMoistRemv, double latload)
                                   {
	const double LEWIS23 = pow(0.919, 2.0/3.0);	// Lewis number for air and water vapor from ASHRAE 1989 p5-9 to the 2/3rds power
	const double DIFFCOEF = 3E-10; // diffusion coefficient for pine from Cunningham 1990 (m2/s)
	const int RWATER = 462;			// gas constant for water vapor (J/KgK)
	double PWOut;						// outdoor air vapor pressure (Pa)
	double hw;							// mass transfer coefficient for water vapor (m/s)
	
	// set node temperatures 
	temperature[0] = node_temps[3];	// inner south sheathing
	temperature[1] = node_temps[1];	// inner north sheathing
	temperature[2] = node_temps[5];	// bulk wood
	temperature[3] = (node_temps[3] + node_temps[4]) / 2;	// average of south sheathing
	temperature[4] = (node_temps[1] + node_temps[2]) / 2;	// average of north sheathing
	temperature[5] = node_temps[5];	// bulk wood
	temperature[6] = node_temps[0];	// attic air
	temperature[7] = node_temps[11];	// return air
	temperature[8] = node_temps[14];	// supply air
	temperature[9] = node_temps[15];	// house air
	temperature[10] = node_temps[12];	// house mass

	PWOut = saturationVaporPressure(tempOut) * RHOut / 100;

	//NODE  0 ON INSIDE OF SOUTH SHEATHING
	hw = hU0 / CpAir / LEWIS23 / airDensityAttic;
	kappa1[0] = calc_kappa_1(pressure, tempOld[0], moistureContent[0], volume[0]);
	kappa2[0] = calc_kappa_2(moistureContent[0], volume[0]);
	x1 = hw * area[0] / RWATER / temperature[0];
	x2 = DIFFCOEF * area[0] / RWATER / temperature[0] / deltaX[0];
	x3 = -hw * area[0] / RWATER / temperature[6];
	x4 = -DIFFCOEF * area[0] / RWATER / temperature[3] / deltaX[0];

	A[0][0] = kappa1[0] + x1 + x2;
	A[0][6] = x3;
	A[0][3] = x4;
	PWInit[0] = kappa1[0] * PWOld[0] - kappa2[0] * (temperature[0] - tempOld[0]);


	//NODE 1 INSIDE OF NORTH SHEATHING
	hw = hU1 / CpAir / LEWIS23 / airDensityAttic;
	kappa1[1] = calc_kappa_1(pressure, tempOld[1], moistureContent[1], volume[1]);
	kappa2[1] = calc_kappa_2(moistureContent[1], volume[1]);
	x5 = hw * area[1] / RWATER / temperature[1];
	x6 = DIFFCOEF * area[1] / RWATER / temperature[1] / deltaX[1];
	x7 = -hw * area[1] / RWATER / temperature[6];
	x8 = -DIFFCOEF * area[1] / RWATER / temperature[4] / deltaX[1];
	A[1][1] = kappa1[1] + x5 + x6;
	A[1][6] = x7;
	A[1][4] = x8;
	PWInit[1] = kappa1[1] * PWOld[1] - kappa2[1] * (temperature[1] - tempOld[1]);

	//NODE 2 IS OUTSIDE mass OF WOOD IN ATTIC JOISTS AND TRUSSES
   hw = hU2 / CpAir / LEWIS23 / airDensityAttic;
	kappa1[2] = calc_kappa_1(pressure, tempOld[2], moistureContent[2], volume[2]);
	kappa2[2] = calc_kappa_2(moistureContent[2], volume[2]);
	x9 = hw * area[2] / RWATER / temperature[2];
	x10 = DIFFCOEF * area[2] / RWATER / temperature[2] / deltaX[2];
	x11 = -hw * area[2] / RWATER / temperature[6];
	x12 = -DIFFCOEF * area[2] / RWATER / temperature[5] / deltaX[2];
	A[2][2] = kappa1[2] + x9 + x10;
	A[2][6] = x11;
	A[2][5] = x12;
	PWInit[2] = kappa1[2] * PWOld[2] - kappa2[2] * (temperature[2] - tempOld[2]);

	//NODE 3 IS  SOUTH SHEATHING
	kappa1[3] = calc_kappa_1(pressure, tempOld[3], moistureContent[3], volume[3]);
	kappa2[3] = calc_kappa_2(moistureContent[3], volume[3]);
	A[3][0] = -x2;
	A[3][3] = kappa1[3] - x4;
	PWInit[3] = kappa1[3] * PWOld[3] - kappa2[3] * (temperature[3] - tempOld[3]);

	//NODE 4 IS NORTH SHEATHING
	kappa1[4] = calc_kappa_1(pressure, tempOld[4], moistureContent[4], volume[4]);
	kappa2[4] = calc_kappa_2(moistureContent[4], volume[4]);
	A[4][1] = -x6;
	A[4][4] = kappa1[4] - x8;
	PWInit[4] = kappa1[4] * PWOld[4] - kappa2[4] * (temperature[4] - tempOld[4]);

	//NODE 5 IS INNER BULK WOOD
	kappa1[5] = calc_kappa_1(pressure, tempOld[5], moistureContent[5], volume[5]);
	kappa2[5] = calc_kappa_2(moistureContent[5], volume[5]);
	A[5][2] = -x10;
	A[5][5] = kappa1[5] - x12;
	PWInit[5] = kappa1[5] * PWOld[5] - kappa2[5] * (temperature[5] - tempOld[5]);

	//NODE 6 IS ATTIC AIR
	x13 = volume[6] / RWATER / temperature[6] / timeStep;
	x14 = -mAtticOut / airDensityAttic / RWATER / temperature[6];
	x15 = volume[6] * PWOld[6] / RWATER / temperature[6] / timeStep;
   x17 = mAtticIn * PWOut / airDensityOut / RWATER / tempOut;
	if(mCeiling < 0) {
		xn67 = (mRetAHoff + mRetLeak) / RWATER / temperature[7] / airDensityRet;
		xn68 = (mSupAHoff + mSupLeak) / RWATER / temperature[8] / airDensitySup;
		xn69 = mCeiling / RWATER / temperature[9] / airDensityHouse;
		x16 = 0;
		}
	else {
		xn67 = mRetLeak / RWATER / temperature[7] / airDensityRet;
		xn68 = mSupLeak / RWATER / temperature[8] / airDensitySup;
		xn69 = 0;
		x16 = (mCeiling + mSupAHoff + mRetAHoff) / RWATER / temperature[6] / airDensityAttic;
		}
	A[6][0] = -x1;
	A[6][1] = -x5;
	A[6][2] = -x9;
	A[6][6] = x13 - x3 - x7 - x11 + x14 + x16;
	A[6][7] = xn67;
	A[6][8] = xn68;
	A[6][9] = xn69;
	PWInit[6] = x15 + x17;

	//NODE 7 IS RETURN DUCT AIR
	xn7t = volume[7] / RWATER / temperature[7] / timeStep;
	xn7o = volume[7] * PWOld[7] / RWATER / temperature[7] / timeStep - mRetOut * PWOut / RWATER / tempOut / airDensityOut;
	if(mCeiling < 0) {
		xn7c = (mAH - mRetAHoff) / RWATER / temperature[7] / airDensityRet;
		xn76 = mRetLeak / RWATER / temperature[6] / airDensityAttic;
		xn79 = (mRetReg + mRetAHoff + mErvHouse) / RWATER / temperature[9] / airDensityHouse;
		}
	else {
		xn7c = (mAH + mRetAHoff) / RWATER / temperature[7] / airDensityRet;
		xn76 = (mRetLeak - mRetAHoff) / RWATER / temperature[6] / airDensityAttic;
		xn79 = (mRetReg + mErvHouse) / RWATER / temperature[9] / airDensityHouse;
		}
	A[7][6] = xn76;
	A[7][7] = xn7t + xn7c;
	A[7][9] = xn79;
	PWInit[7] = xn7o;

	//NODE 8 IS SUPPLY DUCT AIR
	xn8t = volume[8] / RWATER / temperature[8] / timeStep;
	xn8o = volume[8] * PWOld[8] / RWATER / temperature[8] / timeStep - latcap / 2501000;
	xn87 = -mAH / RWATER / temperature[7] / airDensityRet;
	if(mCeiling < 0) {
		xn8c = (mSupReg + mSupLeak - mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		xn86 = 0;
		xn89 = mSupAHoff / RWATER / temperature[9] / airDensityHouse;
		}
	else {
		xn8c = (mSupReg + mSupLeak + mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		xn86 = - mSupAHoff / RWATER / temperature[6] / airDensityAttic;
		xn89 = 0;
		}
	A[8][6] = xn86;
	A[8][7] = xn87;
	A[8][8] = xn8t + xn8c;
	A[8][9] = xn89;
	PWInit[8] = xn8o;

	//NODE 9 IS HOUSE AIR
	xn9t = volume[9] / RWATER / temperature[9] / timeStep + haHouse / pressure;
	xn9o = volume[9] * PWOld[9] / RWATER / temperature[9] / timeStep + mHouseIn * PWOut / RWATER / tempOut / airDensityOut - dhMoistRemv + latload;
	xn910 = -haHouse / pressure;
	if(mCeiling < 0) {
		xn9c = (-mHouseOut - mRetReg - mCeiling - mRetAHoff - mSupAHoff) / RWATER / temperature[9] / airDensityHouse;
		xn96 = 0;
		xn97 = 0;
		xn98 = -mSupReg / RWATER / temperature[8] / airDensitySup;
		}
	else {
		xn9c = (-mHouseOut - mRetReg) / RWATER / temperature[9] / airDensityHouse;
		xn96 = -mCeiling / RWATER / temperature[6] / airDensityAttic;
		xn97 = -mRetAHoff / RWATER / temperature[7] / airDensityRet;
		xn98 = (-mSupReg - mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		}
	A[9][6] = xn96;
	A[9][7] = xn97;
	A[9][8] = xn98;
	A[9][9] = xn9t + xn9c;
	A[9][10] = xn910;
	PWInit[9] = xn9o;

	//NODE 10 IS HOUSE MASS
	xn10t = massWHouse / pressure / timeStep + haHouse / pressure;
	xn10o = massWHouse * PWOld[10] / pressure / timeStep;
	xn109 = -haHouse / pressure;
	A[10][9] = xn109;
	A[10][10] = xn10t;
	PWInit[10] = xn10o;

   for (int i=0; i<MOISTURE_NODES; i++) {
       A[i][MOISTURE_NODES] = PWInit[i];
       }
   PW = gauss(A);
	
	// once the attic moisture nodes have been calculated assuming no condensation (as above)
	// then we call cond_bal to check for condensation and redo the calculations if necessasry
	// cond_bal performs an iterative scheme that tests for condensation
	cond_bal(pressure);

	for(int i=0; i<MOISTURE_NODES; i++) {
		tempOld[i] = temperature[i];
		PWOld[i] = abs(PW[i]);			// in case negative PW's are generated (should get set to 0??)
		}
}

void print_matrix(vector< vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            printf("%5e\t",A[i][j]);
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
    cout << endl;
}


/*
 * cond_bal - checks for condensation at the various nodes and corrects the mass balances
 *            to adapt for this condensation.
 * @param pressure - atomospheric pressure (Pa)
 */
void Moisture::cond_bal(int pressure) {
	double PWSaturation[MOISTURE_NODES];			// saturation vapor pressure at each node (Pa)
	double PWTest[MOISTURE_NODES];					// save original for test of convergence (Pa)
	double massCondensed[MOISTURE_NODES];			// the mass condensed or removed (if negative) at each node (kg)
	double fluxTotalOld = 0;
	bool hasCondensedMass[MOISTURE_NODES];			// flag indicates that node has condensed mass
	bool redoMassBalance;								// flag for condensed mass loop exit
	bool PWOutOfRange;									// flag for PW loop exit
	int outIter;											// number of outer loop iterations
	int inIter;												// number of inner loop iterations

	// Calculate saturation vapor pressure for each node (moved from mass_cond_bal)
	for(int i=0; i<MOISTURE_NODES; i++) {
		PWSaturation[i] = saturationVaporPressure(temperature[i]);
		}
			
//MASSTOOUT = 0
	outIter = 0;
	do {
		outIter++;
		inIter = 0;
		redoMassBalance = false;		// redo calc if a node has gone from having condensed mass @ PWS

		do {		// this loops until every PW is within  0.1 Pa of the previous iteration
			inIter++;
			if(outIter > 1 || inIter > 1)
				cout << "iterations: " << outIter << ":" << inIter << endl;
			for(int i=0; i<MOISTURE_NODES; i++)
				PWTest[i] = PW[i];

	   	// NODE 0:
			if(PW[0] < PWSaturation[0] && mTotal[0] <= 0) {
				// Checking for remaining condensed mass from previous iteration
				// There may still be mass condensed at a node even when PW<PWS 
				mTotal[0] = 0;
				massCondensed[0] = 0;
				double a00 = kappa1[0] + x1 + x2;
				PW[0] = -1 / a00 * (x3 * PW[6] + x4 * PW[3]) + PWInit[0] / a00;
				hasCondensedMass[0] = false;
				}
			else {
				// If there is still condensed mass, then PW is set equal to PWS
				// hasCondensedMass is set for this node to indicate that this node has condensed mass
				// to exchange moisture rather than changing the moisture contant and use the PW/MC relationship
				PW[0] = PWSaturation[0];
				hasCondensedMass[0] = true;
				}

   		// NODE 1:
			if(PW[1] < PWSaturation[1] && mTotal[1] <= 0) {
				mTotal[1] = 0;
				massCondensed[1] = 0;
				double a11 = kappa1[1] + x5 + x6;
				PW[1] = -1 / a11 * (x7 * PW[6] + x8 * PW[4]) + PWInit[1] / a11;
				hasCondensedMass[1] = false;
				}
			else {
            PW[1] = PWSaturation[1];
				hasCondensedMass[1] = true;
				}

   		// NODE 2:
			if(PW[2] < PWSaturation[2] && mTotal[2] <= 0) {
				mTotal[2] = 0;
				massCondensed[2] = 0;
				double a22 = kappa1[2] + x9 + x10;
				PW[2] = -1 / a22 * (x11 * PW[6] + x12 * PW[5]) + PWInit[2] / a22;
				hasCondensedMass[2] = false;
				}
			else {
            PW[2] = PWSaturation[2];
				hasCondensedMass[2] = true;
				}

   		// NODE 3:
			if(PW[3] < PWSaturation[3] && mTotal[3] <= 0) {
				mTotal[3] = 0;
				massCondensed[3] = 0;
				double a44 = kappa1[3] - x4;
				PW[3] = -1 / a44 * (-x2 * PW[0]) + PWInit[3] / a44;
				hasCondensedMass[3] = false;
				}
			else {
            PW[3] = PWSaturation[3];
				hasCondensedMass[3] = true;
				}

   		// NODE 4:
			if(PW[4] < PWSaturation[4] && mTotal[4] <= 0) {
				mTotal[4] = 0;
				massCondensed[4] = 0;
				double a55 = kappa1[4] - x8;
				PW[4] = -1 / a55 * (-x6 * PW[1]) + PWInit[4] / a55;
				hasCondensedMass[4] = false;
				}
			else {
            PW[4] = PWSaturation[4];
				hasCondensedMass[4] = true;
				}

   		// NODE 5:
			if(PW[5] < PWSaturation[5] && mTotal[5] <= 0) {
				mTotal[5] = 0;
				massCondensed[5] = 0;
				double a66 = kappa1[5] - x12;
				PW[5] = -1 / a66 * (-x10 * PW[2]) + PWInit[5] / a66;
				hasCondensedMass[5] = false;
				}
			else {
            PW[5] = PWSaturation[5];
				hasCondensedMass[5] = true;
				}

			// NODE 6: the attic node is treated differently as we do not have accumulating mass of moisture for the attic air - 
			//         i.e., there is no mTotal
			if(PW[6] < PWSaturation[6]) {
				massCondensed[6] = 0;
				double a33 = x13 - x3 - x7 - x11 + x14 + x16;
				PW[6] = -1 / a33 * (-x1 * PW[0] - x5 * PW[1] - x9 * PW[2] + xn67 * PW[7] + xn68 * PW[8] + xn69 * PW[9]) + PWInit[6] / a33;
				hasCondensedMass[6] = false;
				}
			else {
            PW[6] = PWSaturation[6];
				hasCondensedMass[6] = true;
				}
			
			// Other air nodes - do we need to recalc PW? just fix PW at saturation for now as there is no place to put the moisture
			for(int i=7; i<MOISTURE_NODES; i++) {
				//hasCondensedMass[i] = false;
				if(PW[i] > PWSaturation[i]) {
				   cout << "Air node " << i << " > saturation: " << PW[i] << "> " << PWSaturation[i] << endl;
            	PW[i] = PWSaturation[i];
            	}
            }

			PWOutOfRange = false;
			for(int i=0; i<MOISTURE_NODES; i++) {
				if(abs(PWTest[i] - PW[i]) > 0.1) {
					// it is possible that this 0.1Pa criterion is too “tight” and we may relax this 
					// depending on the solution time requirements and program stability.
					PWOutOfRange = true;
					}
				}
			} while (PWOutOfRange);
		
		// the following is the special case for condensation at the attic air node
		if(hasCondensedMass[6]) {
			double fluxTotal = 0;
			double fluxTo[4];
			massCondensed[6] = timeStep * (-A[6][6] * PW[6] + x1 * PW[0] + x5 * PW[1] + x9 * PW[2] + PWInit[6]);
			// This mass must be distributed over the other attic surfaces
			// this may push them over pws so the iterative process must be repeated
			// by calculating a new pw for each of the surface nodes based on their
			// moisture content including this transferred from the attic air
			// here the fluxes of water from the air to the other surface nodes in contact with attic air is determined
			fluxTo[0] = -(x1 * PW[0] + x3 * PW[6]);
			fluxTo[1] = -(x5 * PW[1] + x7 * PW[6]);
			fluxTo[2] = -(x9 * PW[2] + x11 * PW[6]);
			fluxTo[3] = x14 * PW[6];		// this is the moisture in attic air going to and from outside/house
			for(int i=0; i < 4; i++) {
				if(fluxTo[i] > 0)
					fluxTotal = fluxTotal + fluxTo[i];
				}

			if(abs(fluxTotalOld - fluxTotal) >= .001) { 
				// if not converged, then the node moisture fluxes are recalculated to include the extra moisture from the attic air
				for(int i=0; i < 3; i++) {
					if(fluxTo[i] > 0) {
						moistureContent[i] = moistureContent[i] + (fluxTo[i] / fluxTotal * massCondensed[6]) / rhoWood / volume[i];
						PW[i] = calc_vapor_pressure(moistureContent[i], temperature[i], pressure);
						if(PW[i] > PWSaturation[i]) {
							PW[i] = PWSaturation[i];
							double mcSat = mc_cubic(PW[i], pressure, temperature[i]);
							mTotal[i] = mTotal[i] + (moistureContent[i] - mcSat) * volume[i] * rhoWood;
							moistureContent[i] = mcSat;
							}
						}
					}
				//      IF FLUXTOOUT > 0 THEN
				//         MASSTOOUT = FLUXTOOUT / FLUXTOT * MCOND(4)
				//      END IF
				fluxTotalOld = fluxTotal;
				redoMassBalance = true;
				continue;	// GOTO REDOIT
				}
			}
		//ENDITER:

		// Calculate MC's for wood nodes based on these new PW's
		// this uses the inverse of cleary's relationship
		for(int i=0; i < 6; i++) {
			moistureContent[i] = max(mc_cubic(PW[i], pressure, temperature[i]), 0.032);
			if(hasCondensedMass[i]) {
				switch (i) {
				case 0:
					massCondensed[0] = timeStep * (-kappa1[0] * (PW[0] - PWOld[0]) - kappa2[0] * (temperature[0] - tempOld[0]) - x1 * PW[0] - x3 * PW[6] - x2 * PW[0] - x4 * PW[3]);
					break;
				case 1:
					massCondensed[1] = timeStep * (-kappa1[1] * (PW[1] - PWOld[1]) - kappa2[1] * (temperature[1] - tempOld[1]) - x5 * PW[1] - x7 * PW[6] - x6 * PW[1] - x8 * PW[4]);
					break;
				case 2:
					massCondensed[2] = timeStep * (-kappa1[2] * (PW[2] - PWOld[2]) - kappa2[2] * (temperature[2] - tempOld[2]) - x9 * PW[2] - x11 * PW[6] - x10 * PW[2] - x12 * PW[5]);
					break;
				case 3:
					massCondensed[3] = timeStep * (-kappa1[3] * (PW[3] - PWOld[3]) - kappa2[3] * (temperature[3] - tempOld[3]) + x2 * PW[0] + x4 * PW[3]);
					break;
				case 4:
					massCondensed[4] = timeStep * (-kappa1[4] * (PW[4] - PWOld[4]) - kappa2[4] * (temperature[4] - tempOld[4]) + x6 * PW[1] + x8 * PW[4]);
					break;
				case 5:
					massCondensed[5] = timeStep * (-kappa1[5] * (PW[5] - PWOld[5]) - kappa2[5] * (temperature[5] - tempOld[5]) + x10 * PW[2] + x12 * PW[5]);
					break;
				}
				mTotal[i] = mTotal[i] + massCondensed[i];
				if(mTotal[i] < 0) {
					// Need to estimate how much MC (and PW) change. A reduction in PW means mass condensed
					// even more -ve so an iterative scheme is not needed
					double moistureReduction = mTotal[i] / (volume[i] * rhoWood);
					moistureContent[i] = max(moistureContent[i] + moistureReduction, 0.032);
					PW[i] = calc_vapor_pressure(moistureContent[i], temperature[i], pressure);
					redoMassBalance = true;	// Redo mass balance with this PW because this node is no longer at saturation. 
													// this may become unnecessary at shorter time steps when the error due 
													// to lost condensed mass not accounted for is small
					}
				}
			}
		// Set air node RHs (%)
		for(int i=6; i<MOISTURE_NODES; i++) {
			moistureContent[i] = PW[i] / PWSaturation[i] * 100;
			}
		// house mass
		//moistureContent[MOISTURE_NODES-1] = max(mc_cubic(PW[MOISTURE_NODES-1], pressure, temperature[MOISTURE_NODES-1]), 0.032);
		} while(redoMassBalance);
}

// humidity ratio constants (from Cleary 1985)
static const double B3 =  15.8;		// was MCA
static const double B4 = -0.0015;	// was MCB
static const double B5 =  0.053;		// was MCC
static const double B6 = -0.184;		// was MCD
static const double B7 =  0.233;		// was MCE


/*
 * calc_kappa_1 - calculates kappa1: the partial of wood moisture content with respect to
 *						pressure (equation 4-15) times wood mass divided by tau
 * @param pressure - atomospheric pressure (Pa)
 * @param temp     - temperature (deg K)
 * @param mc       - moisture content
 * @param volume   - wood volume (m3)
 * @return kappa1  - (dWmc/dP)(Mw/t)
 */
double Moisture::calc_kappa_1(int pressure, double temp, double mc, double volume) {
	double k1;
	k1 = 1 / (pressure / 0.622 * exp((temp - C_TO_K) / B3) * (B5 + 2 * B6 * mc + 3 * B7 * pow(mc,2)));
	k1 = volume * rhoWood * k1 / timeStep;
//cout << "kappa1 in: p=" << pressure << " temp=" << temp << " mc=" << mc << " volume=" << volume << " out: " << k1 << endl;
	return (k1);
}

/*
 * calc_kappa_2 - calculates kappa2: the partial of wood moisture content with respect to
 *						temperature (equation 4-16) times wood mass divided by tau
 * @param mc		- moisture content
 * @param volume  - wood volume (m3)
 * @return kappa2 - (dWmc/dT)(Mw/t)
 */
double Moisture::calc_kappa_2(double mc, double volume) {
	double k2;
	k2 = (B4 + B5 * mc + B6 * pow(mc,2) + B7 * pow(mc,3)) / (-B3 * (B5 + 2 * B6 * mc + 3 * B7 * pow(mc,2)));
	k2 = volume * rhoWood * k2 / timeStep;
	return (k2);
}

/*
 * mc_cubic - cubic solver to get MC from PW
 * @param pw		 - vapor pressure
 * @param pressure - atmospheric pressure (Pa)
 * @param temp     - temperature (deg K)
 * @return mc		 - wood moisture content
 */
double Moisture::mc_cubic(double pw, int pressure, double temp) {
	double W = 0.622 * pw / (pressure - pw);
	const double a1 = B6 / B7;
	const double a2 = B5 / B7;
	double a3 = (B4 - W / exp((temp - C_TO_K) / B3)) / B7;
	double q = (3 * a2 - pow(a1, 2)) / 9;
	double r = (9 * a1 * a2 - 27 * a3 - 2 * pow(a1, 3)) / 54;
	double disc = pow(q, 3) + pow(r, 2);
	double s = pow(r + pow(disc, 0.5), 1.0/3.0);
	double t = -pow(abs(r - pow(disc,0.5)), 1.0/3.0);
//cout << "mc_cubic in: pw=" << pw << " pressure=" << pressure << " temp=" << temp  << " out: " << s + t - a1 / 3 << endl;
	return (s + t - a1 / 3);
	}

/*
 * calc_pw - calculates vapor pressure based on moisture content and temperature
 * @param mc		 - wood moisture content
 * @param temp     - temperature (deg K)
 * @param pressure - atmospheric pressure (Pa)
 * @return pw      - vapor pressure
 */
double Moisture::calc_vapor_pressure(double mc, double temp, int pressure) {
	double w = exp((temp - C_TO_K) / B3) * (B4 + B5 * mc + B6 * pow(mc, 2) + B7 * pow(mc, 3));
	return (w * pressure / (0.622 * (1 + w / 0.622)));
	}
