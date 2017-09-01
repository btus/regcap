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
 * @param sheathArea - roof sheathing area (m2)
 * @param bulkArea - bulk wood area (m2)
 * @param roofInsThick - thickness of interior roof insulation (m)
 * @param roofExtRval - R-value of exterior roof insulation (m2k/W)
 * @param mcInit - initial wood moisture content (fraction - optional)
 *
 * The nodes are (0-5 wood, 6-10 air/mass):
 * 0. surface of south sheathing (facing the attic) - corresponding thermal node is (3)
 * 1. surface of north sheathing (facing the attic) - corresponding thermal node is (1)
 * 2. surface of the bulk wood in the attic - corresponding thermal node is (5)
 * 3. inside south sheathing - corresponding thermal node is ((3+4)/2)
 * 4. inside north sheathing - corresponding thermal node is ((1+2)/2)
 * 5. inside bulk wood - corresponding thermal node is (5)
 * 6. attic air - corresponding thermal node is (0)
 * 7. return air - corresponding thermal node is (11)
 * 8. supply air - corresponding thermal node is (14)
 * 9. house air - corresponding thermal node is (15)
 * 10. house materials - corresponding thermal node is (12)
 * 11. interior south roof insulation ((3+17)/2)
 * 11. interior north roof insulation ((1+16)/2)
 */
Moisture::Moisture(double atticVolume, double retVolume, double supVolume, double houseVolume, 
						double floorArea, double sheathArea, double bulkArea, double roofInsThick, double roofExtRval, double mcInit) {
   int pressure = 101325;		// 1 atmosphere
   double sheathThick = 0.015;  // Roof sheathing thickness (m).
   double bulkThick = 0.013;    // Bulk wood thickness (m).
   double sheathSurfThick = 0.003; // Sheathing surface wood layer thickness (m)
   double bulkSurfThick = 0.001; // Bulk surface wood layer thickness (m)
	double tempInit = airTempRef; // node initialization temperature (deg K)
	double RHInit = 50;	// initial air RH (%)

	if(roofInsThick > 0) {
		moisture_nodes = 13;
		deltaX[11] = roofInsThick/2;
		deltaX[12] = deltaX[11];
		area[11] = sheathArea;
		area[12] = area[11];
		volume[11] = area[11] * roofInsThick;
		volume[12] = area[12] * roofInsThick;
		}
	else {
		moisture_nodes = 11;
		}
	if(roofExtRval > 0) {
	   roofInsulRatio = (1/roofExtRval)/(1/roofExtRval + 1/((thickSheathing / kWood) + 0.06));
	   }
	else {
	   roofInsulRatio = 1;
	   }
   A.resize(moisture_nodes, vector<double>(moisture_nodes+1, 0));		// set size of equation vectors to number of nodes (A contains both)
   PW.resize(moisture_nodes, 0);

	deltaX[0] = sheathThick / 2;                // distance between the centers of the surface and inside wood layers. it is half the characteristic thickness of the wood member
	deltaX[1] = deltaX[0];                    // same as for the other sheathing surface
	deltaX[2] = bulkThick / 2;
	deltaX[3] = deltaX[0];
	deltaX[4] = deltaX[1];
	deltaX[5] = deltaX[2];

	area[0] = sheathArea;
	area[1] = sheathArea;
	area[2] = bulkArea;
	area[3] = area[0];
	area[4] = area[1];
	area[5] = area[2];

	volume[0] = area[0] * sheathSurfThick;
	volume[1] = area[1] * sheathSurfThick;
	volume[2] = area[2] * bulkSurfThick;
	volume[3] = area[0] * sheathThick - volume[0];
	volume[4] = area[1] * sheathThick - volume[1];
	volume[5] = area[2] * bulkThick - volume[2];
	volume[6] = atticVolume;
	volume[7] = retVolume;
	volume[8] = supVolume;
	volume[9] = houseVolume;

	density[0] = densitySheathing;
	density[1] = densitySheathing;
	density[2] = densityWood;
	density[3] = densitySheathing;
	density[4] = densitySheathing;
	density[5] = densityWood;

	haHouse = 1 * 0.622 * .5 * floorArea / 186;	// moisture transport coefficient scales with floor area (kg/s) - 0.622 (dHR/dVP), 186 (area of std house in m2), 0.5 empirical coefficient (kg/s)
	massWHouse = 1 * 0.622 * 60 * floorArea;					// active mass containing moisture in the house (kg) - empirical

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
 * @param hU0 - surface heat transfer coefficient for node 0/11 (inner south) (J/sm2K) - was H4
 * @param hU1 - surface heat transfer coefficient for node 1/12 (inner north) (J/sm2K) - was H6
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
	const double LEWIS = 0.919;	// Lewis number for air and water vapor from ASHRAE 1989 p5-9
	const double diffCoefWood = 3E-10; // diffusion coefficient for pine from Cunningham 1990 (m2/s)
	const double diffCoefIns = 2.12E-5; // diffusion coefficient for fiberglass (106 perm-in) (m2/s)
	const int RWATER = 462;			// gas constant for water vapor (J/KgK)
	double PWOut;						// outdoor air vapor pressure (Pa)
	double hw0, hw1, hw2;			// mass transfer coefficients for water vapor (m/s)

	// set node temperatures 
	temperature[0] = node_temps[3];	// surface of south sheathing
	temperature[1] = node_temps[1];	// surface of north sheathing
	temperature[2] = node_temps[5];	// surface of bulk wood
	temperature[3] = calc_inter_temp(node_temps[3], node_temps[4], roofInsulRatio);	// interior of south sheathing
	temperature[4] = calc_inter_temp(node_temps[1], node_temps[2], roofInsulRatio);	// interior of north sheathing
	temperature[5] = node_temps[5];	// bulk wood
	temperature[6] = node_temps[0];	// attic air
	temperature[7] = node_temps[11];	// return air
	temperature[8] = node_temps[14];	// supply air
	temperature[9] = node_temps[15];	// house air
	temperature[10] = node_temps[12];	// house mass
	if(moisture_nodes > 11) {
		temperature[11] = calc_inter_temp(node_temps[3], node_temps[17], 1);	// middle of south interior insulation
		temperature[12] = calc_inter_temp(node_temps[1], node_temps[16], 1);	// middle of north interior insulation
		}

	PWOut = saturationVaporPressure(tempOut) * RHOut / 100;
	
	// water vapor mass transfer coefficients
	double hu_conv = 1 / CpAir / pow(LEWIS, 2.0/3.0) / airDensityAttic;
	hw0 = hU0 * hu_conv;
	hw1 = hU1 * hu_conv;
   hw2 = hU2 * hu_conv;

	//NODE  0 inside of south sheathing
	kappa1[0] = calc_kappa_1(pressure, tempOld[0], moistureContent[0], volume[0] * density[0]);
	kappa2[0] = calc_kappa_2(moistureContent[0], volume[0] * density[0]);
	x30 = diffCoefWood * area[0] / RWATER / temperature[0] / deltaX[0];
	x03 = -diffCoefWood * area[0] / RWATER / temperature[3] / deltaX[0];
	A[0][0] = kappa1[0] + x30;
	A[0][3] = x03;
	PWInit[0] = kappa1[0] * PWOld[0] - kappa2[0] * (temperature[0] - tempOld[0]);
	if(moisture_nodes > 11) {		// there is interior insulation
		x110 = diffCoefIns * area[0] / RWATER / temperature[0] / deltaX[11];
		x011 = -diffCoefIns * area[0] / RWATER / temperature[11] / deltaX[11];
		A[0][0] += x110;
		A[0][11] = x011;
		}
	else {		// exposed to attic air
		x60 = hw0 * area[0] / RWATER / temperature[0];
		x06 = -hw0 * area[0] / RWATER / temperature[6];
		A[0][0] += x60;
		A[0][6] = x06;
		}

	//NODE 1 inside of north sheathing
	kappa1[1] = calc_kappa_1(pressure, tempOld[1], moistureContent[1], volume[1] * density[1]);
	kappa2[1] = calc_kappa_2(moistureContent[1], volume[1] * density[1]);
	x41 = diffCoefWood * area[1] / RWATER / temperature[1] / deltaX[1];
	x14 = -diffCoefWood * area[1] / RWATER / temperature[4] / deltaX[1];
	A[1][1] = kappa1[1] + x41;
	A[1][4] = x14;
	PWInit[1] = kappa1[1] * PWOld[1] - kappa2[1] * (temperature[1] - tempOld[1]);
	if(moisture_nodes > 11) {		// there is interior insulation
		x121 = diffCoefIns * area[1] / RWATER / temperature[1] / deltaX[12];
		x112 = -diffCoefIns * area[1] / RWATER / temperature[12] / deltaX[12];
		A[1][1] += x121;
		A[1][12] = x112;
		}
	else {		// exposed to attic air
		x61 = hw1 * area[1] / RWATER / temperature[1];
		x16 = -hw1 * area[1] / RWATER / temperature[6];
		A[1][1] += x61;
		A[1][6] = x16;
		}

	//NODE 2 IS OUTSIDE mass OF WOOD IN ATTIC JOISTS AND TRUSSES
	kappa1[2] = calc_kappa_1(pressure, tempOld[2], moistureContent[2], volume[2] * density[2]);
	kappa2[2] = calc_kappa_2(moistureContent[2], volume[2] * density[2]);
	x62 = hw2 * area[2] / RWATER / temperature[2];
	x52 = diffCoefWood * area[2] / RWATER / temperature[2] / deltaX[2];
	x26 = -hw2 * area[2] / RWATER / temperature[6];
	x25 = -diffCoefWood * area[2] / RWATER / temperature[5] / deltaX[2];
	A[2][2] = kappa1[2] + x62 + x52;
	A[2][6] = x26;
	A[2][5] = x25;
	PWInit[2] = kappa1[2] * PWOld[2] - kappa2[2] * (temperature[2] - tempOld[2]);

	//NODE 3 IS  SOUTH SHEATHING
	kappa1[3] = calc_kappa_1(pressure, tempOld[3], moistureContent[3], volume[3] * density[3]);
	kappa2[3] = calc_kappa_2(moistureContent[3], volume[3] * density[3]);
	A[3][0] = -x30;
	A[3][3] = kappa1[3] - x03;
	PWInit[3] = kappa1[3] * PWOld[3] - kappa2[3] * (temperature[3] - tempOld[3]);

	//NODE 4 IS NORTH SHEATHING
	kappa1[4] = calc_kappa_1(pressure, tempOld[4], moistureContent[4], volume[4] * density[4]);
	kappa2[4] = calc_kappa_2(moistureContent[4], volume[4] * density[4]);
	A[4][1] = -x41;
	A[4][4] = kappa1[4] - x14;
	PWInit[4] = kappa1[4] * PWOld[4] - kappa2[4] * (temperature[4] - tempOld[4]);

	//NODE 5 IS INNER BULK WOOD
	kappa1[5] = calc_kappa_1(pressure, tempOld[5], moistureContent[5], volume[5] * density[5]);
	kappa2[5] = calc_kappa_2(moistureContent[5], volume[5] * density[5]);
	A[5][2] = -x52;
	A[5][5] = kappa1[5] - x25;
	PWInit[5] = kappa1[5] * PWOld[5] - kappa2[5] * (temperature[5] - tempOld[5]);

	//NODE 6 IS ATTIC AIR
	x66 = volume[6] / RWATER / temperature[6] / timeStep;
	x6out = -mAtticOut / airDensityAttic / RWATER / temperature[6];
	if(mCeiling < 0) {
		x67 = (mRetAHoff + mRetLeak) / RWATER / temperature[7] / airDensityRet;
		x68 = (mSupAHoff + mSupLeak) / RWATER / temperature[8] / airDensitySup;
		x69 = mCeiling / RWATER / temperature[9] / airDensityHouse;
		}
	else {
		x67 = mRetLeak / RWATER / temperature[7] / airDensityRet;
		x68 = mSupLeak / RWATER / temperature[8] / airDensitySup;
		x69 = 0;
		x66 += (mCeiling + mSupAHoff + mRetAHoff) / RWATER / temperature[6] / airDensityAttic;
		}
	A[6][6] = x66 - x26 + x6out;
	A[6][7] = x67;
	A[6][8] = x68;
	A[6][9] = x69;
	PWInit[6] = volume[6] * PWOld[6] / RWATER / temperature[6] / timeStep + mAtticIn * PWOut / airDensityOut / RWATER / tempOut;
	if(moisture_nodes > 11) {
		x611 = -hw0 * area[11] / RWATER / temperature[11];
		x116 = hw0 * area[11] / RWATER / temperature[6];
		x612 = -hw1 * area[12] / RWATER / temperature[12];
		x126 = hw1 * area[12] / RWATER / temperature[6];
		A[6][2] = -x62;
		A[6][6] += x116 + x126;
		A[6][11] = x611;
		A[6][12] = x612;
		}
	else {
		A[6][0] = -x60;
		A[6][1] = -x61;
		A[6][2] = -x62;
		A[6][6] += -x06 - x16;
		}

	//NODE 7 IS RETURN DUCT AIR
	A[7][7] = volume[7] / RWATER / temperature[7] / timeStep;
	if(mCeiling < 0) {
		A[7][7] += (mAH - mRetAHoff) / RWATER / temperature[7] / airDensityRet;
		A[7][6] = mRetLeak / RWATER / temperature[6] / airDensityAttic;
		A[7][9] = (mRetReg + mRetAHoff + mErvHouse) / RWATER / temperature[9] / airDensityHouse;
		}
	else {
		A[7][7] += (mAH + mRetAHoff) / RWATER / temperature[7] / airDensityRet;
		A[7][6] = (mRetLeak - mRetAHoff) / RWATER / temperature[6] / airDensityAttic;
		A[7][9] = (mRetReg + mErvHouse) / RWATER / temperature[9] / airDensityHouse;
		}
	PWInit[7] = volume[7] * PWOld[7] / RWATER / temperature[7] / timeStep - mRetOut * PWOut / RWATER / tempOut / airDensityOut;

	//NODE 8 IS SUPPLY DUCT AIR
	A[8][8] = volume[8] / RWATER / temperature[8] / timeStep;
	A[8][7] = -mAH / RWATER / temperature[7] / airDensityRet;
	if(mCeiling < 0) {
		A[8][8] += (mSupReg + mSupLeak - mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		A[8][6] = 0;
		A[8][9] = mSupAHoff / RWATER / temperature[9] / airDensityHouse;
		}
	else {
		A[8][8] += (mSupReg + mSupLeak + mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		A[8][6] = - mSupAHoff / RWATER / temperature[6] / airDensityAttic;
		A[8][9] = 0;
		}
	PWInit[8] = volume[8] * PWOld[8] / RWATER / temperature[8] / timeStep - latcap / 2501000;

	//NODE 9 IS HOUSE AIR
	A[9][9] = volume[9] / RWATER / temperature[9] / timeStep + haHouse / pressure;
	A[9][10] = -haHouse / pressure;
	if(mCeiling < 0) {
		A[9][9] += (-mHouseOut - mRetReg - mCeiling - mRetAHoff - mSupAHoff) / RWATER / temperature[9] / airDensityHouse;
		A[9][6] = 0;
		A[9][7] = 0;
		A[9][8] = -mSupReg / RWATER / temperature[8] / airDensitySup;
		}
	else {
		A[9][9] += (-mHouseOut - mRetReg) / RWATER / temperature[9] / airDensityHouse;
		A[9][6] = -mCeiling / RWATER / temperature[6] / airDensityAttic;
		A[9][7] = -mRetAHoff / RWATER / temperature[7] / airDensityRet;
		A[9][8] = (-mSupReg - mSupAHoff) / RWATER / temperature[8] / airDensitySup;
		}
	PWInit[9] = volume[9] * PWOld[9] / RWATER / temperature[9] / timeStep + mHouseIn * PWOut / RWATER / tempOut / airDensityOut - dhMoistRemv + latload;

	//NODE 10 IS HOUSE MASS
	A[10][9] = -haHouse / pressure;
	A[10][10] = massWHouse / pressure / timeStep + haHouse / pressure;
	PWInit[10] = massWHouse * PWOld[10] / pressure / timeStep;

	if(moisture_nodes > 11) {
		//NODE 11  - inside of south roof insulation
		A[11][11] = volume[11] / RWATER / temperature[11] / timeStep - x611 - x011;
		A[11][6] = -x116;
		A[11][0] = -x110;
		PWInit[11] = volume[11] * PWOld[11] / RWATER / temperature[11] / timeStep;

		//NODE 12  - inside of north roof insulation
		A[12][12] = volume[12] / RWATER / temperature[12] / timeStep - x612 - x112;
		A[12][6] = -x126;
		A[12][1] = -x121;
		PWInit[12] = volume[12] * PWOld[12] / RWATER / temperature[12] / timeStep;
		}


   for (int i=0; i<moisture_nodes; i++) {
       A[i][moisture_nodes] = PWInit[i];
       }
   PW = gauss(A);

	// once the attic moisture nodes have been calculated assuming no condensation (as above)
	// then we call cond_bal to check for condensation and redo the calculations if necessasry
	// cond_bal performs an iterative scheme that tests for condensation
	cond_bal(pressure);

	for(int i=0; i<moisture_nodes; i++) {
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
	for(int i=0; i<moisture_nodes; i++) {
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
			for(int i=0; i<moisture_nodes; i++)
				PWTest[i] = PW[i];

	   	// NODE 0:
			if(PW[0] < PWSaturation[0] && mTotal[0] <= 0) {
				// Checking for remaining condensed mass from previous iteration
				// There may still be mass condensed at a node even when PW<PWS 
				mTotal[0] = 0;
				massCondensed[0] = 0;
				hasCondensedMass[0] = false;
				if(moisture_nodes > 11) {
					double a00 = kappa1[0] + x110 + x30;
					PW[0] = -1 / a00 * (x011 * PW[11] + x03 * PW[3]) + PWInit[0] / a00;
					}
				else {
					double a00 = kappa1[0] + x60 + x30;
					PW[0] = -1 / a00 * (x06 * PW[6] + x03 * PW[3]) + PWInit[0] / a00;
					}
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
				hasCondensedMass[1] = false;
				if(moisture_nodes > 11) {
					double a11 = kappa1[1] + x121 + x41;
					PW[1] = -1 / a11 * (x112 * PW[12] + x14 * PW[4]) + PWInit[1] / a11;
					}
				else {
					double a11 = kappa1[1] + x61 + x41;
					PW[1] = -1 / a11 * (x16 * PW[6] + x14 * PW[4]) + PWInit[1] / a11;
					}
				}
			else {
            PW[1] = PWSaturation[1];
				hasCondensedMass[1] = true;
				}

   		// NODE 2:
			if(PW[2] < PWSaturation[2] && mTotal[2] <= 0) {
				mTotal[2] = 0;
				massCondensed[2] = 0;
				double a22 = kappa1[2] + x62 + x52;
				PW[2] = -1 / a22 * (x26 * PW[6] + x25 * PW[5]) + PWInit[2] / a22;
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
				double a44 = kappa1[3] - x03;
				PW[3] = -1 / a44 * (-x30 * PW[0]) + PWInit[3] / a44;
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
				double a55 = kappa1[4] - x14;
				PW[4] = -1 / a55 * (-x41 * PW[1]) + PWInit[4] / a55;
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
				double a66 = kappa1[5] - x25;
				PW[5] = -1 / a66 * (-x52 * PW[2]) + PWInit[5] / a66;
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
				hasCondensedMass[6] = false;
				if(moisture_nodes > 11) {
					double a33 = x66 + x116 + x126 - x26 + x6out;
					PW[6] = -1 / a33 * (x611 * PW[11] + x612 * PW[12] - x62 * PW[2] + x67 * PW[7] + x68 * PW[8] + x69 * PW[9]) + PWInit[6] / a33;
					}
				else {
					double a33 = x66 - x06 - x16 - x26 + x6out;
					PW[6] = -1 / a33 * (-x60 * PW[0] - x61 * PW[1] - x62 * PW[2] + x67 * PW[7] + x68 * PW[8] + x69 * PW[9]) + PWInit[6] / a33;
					}
				}
			else {
            PW[6] = PWSaturation[6];
				hasCondensedMass[6] = true;
				}
			
			// Other air nodes - do we need to recalc PW? just fix PW at saturation for now as there is no place to put the moisture
			for(int i=7; i<moisture_nodes; i++) {
				//hasCondensedMass[i] = false;
				if(PW[i] > PWSaturation[i]) {
				   cout << "Air node " << i << " > saturation: " << PW[i] << "> " << PWSaturation[i] << endl;
            	PW[i] = PWSaturation[i];
            	}
            }

			PWOutOfRange = false;
			for(int i=0; i<moisture_nodes; i++) {
				if(abs(PWTest[i] - PW[i]) > 0.1) {
					// it is possible that this 0.1Pa criterion is too “tight” and we may relax this 
					// depending on the solution time requirements and program stability.
					PWOutOfRange = true;
					}
				}
			} while (PWOutOfRange);
		
		// the following is the special case for condensation at the attic air node
		if(hasCondensedMass[6]) {
			// This mass must be distributed over the other attic surfaces
			// this may push them over pws so the iterative process must be repeated
			// by calculating a new pw for each of the surface nodes based on their
			// moisture content including this transferred from the attic air
			// here the fluxes of water from the air to the other surface nodes in contact with attic air is determined
			double fluxTotal = 0;
			double fluxTo[4];
			if(moisture_nodes > 11) {
				massCondensed[6] = timeStep * (-A[6][6] * PW[6] + x611 * PW[1] + x612 * PW[12] + x62 * PW[2] + PWInit[6]);
				}
			else {
				massCondensed[6] = timeStep * (-A[6][6] * PW[6] + x60 * PW[0] + x61 * PW[1] + x62 * PW[2] + PWInit[6]);
				}
			fluxTo[0] = -(x60 * PW[0] + x06 * PW[6]);
			fluxTo[1] = -(x61 * PW[1] + x16 * PW[6]);
			fluxTo[2] = -(x62 * PW[2] + x26 * PW[6]);
			fluxTo[3] = x6out * PW[6];		// this is the moisture in attic air going to and from outside/house
			for(int i=0; i < 4; i++) {
				if(fluxTo[i] > 0)
					fluxTotal = fluxTotal + fluxTo[i];
				}

			if(abs(fluxTotalOld - fluxTotal) >= .001) { 
				// if not converged, then the node moisture fluxes are recalculated to include the extra moisture from the attic air
				for(int i=0; i < 3; i++) {
					if(fluxTo[i] > 0) {
						moistureContent[i] = moistureContent[i] + (fluxTo[i] / fluxTotal * massCondensed[6]) / density[i] / volume[i];
						PW[i] = calc_vapor_pressure(moistureContent[i], temperature[i], pressure);
						if(PW[i] > PWSaturation[i]) {
							PW[i] = PWSaturation[i];
							double mcSat = mc_cubic(PW[i], pressure, temperature[i]);
							mTotal[i] = mTotal[i] + (moistureContent[i] - mcSat) * volume[i] * density[i];
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
					if(moisture_nodes > 11) {
						massCondensed[0] = timeStep * (-kappa1[0] * (PW[0] - PWOld[0]) - kappa2[0] * (temperature[0] - tempOld[0]) - x110 * PW[0] - x011 * PW[11] - x30 * PW[0] - x03 * PW[3]);
						}
					else {
						massCondensed[0] = timeStep * (-kappa1[0] * (PW[0] - PWOld[0]) - kappa2[0] * (temperature[0] - tempOld[0]) - x60 * PW[0] - x06 * PW[6] - x30 * PW[0] - x03 * PW[3]);
						}
					break;
				case 1:
					if(moisture_nodes > 11) {
						massCondensed[1] = timeStep * (-kappa1[1] * (PW[1] - PWOld[1]) - kappa2[1] * (temperature[1] - tempOld[1]) - x121 * PW[1] - x112 * PW[12] - x41 * PW[1] - x14 * PW[4]);
						}
					else {
						massCondensed[1] = timeStep * (-kappa1[1] * (PW[1] - PWOld[1]) - kappa2[1] * (temperature[1] - tempOld[1]) - x61 * PW[1] - x16 * PW[6] - x41 * PW[1] - x14 * PW[4]);
						}
					break;
				case 2:
					massCondensed[2] = timeStep * (-kappa1[2] * (PW[2] - PWOld[2]) - kappa2[2] * (temperature[2] - tempOld[2]) - x62 * PW[2] - x26 * PW[6] - x52 * PW[2] - x25 * PW[5]);
					break;
				case 3:
					massCondensed[3] = timeStep * (-kappa1[3] * (PW[3] - PWOld[3]) - kappa2[3] * (temperature[3] - tempOld[3]) + x30 * PW[0] + x03 * PW[3]);
					break;
				case 4:
					massCondensed[4] = timeStep * (-kappa1[4] * (PW[4] - PWOld[4]) - kappa2[4] * (temperature[4] - tempOld[4]) + x41 * PW[1] + x14 * PW[4]);
					break;
				case 5:
					massCondensed[5] = timeStep * (-kappa1[5] * (PW[5] - PWOld[5]) - kappa2[5] * (temperature[5] - tempOld[5]) + x52 * PW[2] + x25 * PW[5]);
					break;
				}
				mTotal[i] = mTotal[i] + massCondensed[i];
				if(mTotal[i] < 0) {
					// Need to estimate how much MC (and PW) change. A reduction in PW means mass condensed
					// even more -ve so an iterative scheme is not needed
					double moistureReduction = mTotal[i] / (volume[i] * density[i]);
					moistureContent[i] = max(moistureContent[i] + moistureReduction, 0.032);
					PW[i] = calc_vapor_pressure(moistureContent[i], temperature[i], pressure);
					redoMassBalance = true;	// Redo mass balance with this PW because this node is no longer at saturation. 
													// this may become unnecessary at shorter time steps when the error due 
													// to lost condensed mass not accounted for is small
					}
				}
			}
		// Set air node RHs (%)
		for(int i=6; i<moisture_nodes; i++) {
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
 * @param mass     - wood mass (kg)
 * @return kappa1  - (dWmc/dP)(Mw/t)
 */
double Moisture::calc_kappa_1(int pressure, double temp, double mc, double mass) {
	double k1;
	k1 = 1 / (pressure / 0.622 * exp((temp - C_TO_K) / B3) * (B5 + 2 * B6 * mc + 3 * B7 * pow(mc,2)));
	k1 = mass * k1 / timeStep;
	return (k1);
}

/*
 * calc_kappa_2 - calculates kappa2: the partial of wood moisture content with respect to
 *						temperature (equation 4-16) times wood mass divided by tau
 * @param mc		- moisture content
 * @param mass    - wood mass (kg)
 * @return kappa2 - (dWmc/dT)(Mw/t)
 */
double Moisture::calc_kappa_2(double mc, double mass) {
	double k2;
	k2 = (B4 + B5 * mc + B6 * pow(mc,2) + B7 * pow(mc,3)) / (-B3 * (B5 + 2 * B6 * mc + 3 * B7 * pow(mc,2)));
	k2 = mass * k2 / timeStep;
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

/*
 * calc_inter_temp - calculates the interior temperature of the roof sheathing
 * @param temp1		- roof sheathing interior surface temp
 * @param temp2    	- roof exterior surface temp
 * @param insRatio	- ratio of exterior insulation U-val to sheathing U-val
 */
double Moisture::calc_inter_temp(double temp1, double temp2, double insRatio) {
   double tInter = temp1 * (1 - insRatio) + temp2 * insRatio;
   return (temp1 + tInter) / 2;
   }