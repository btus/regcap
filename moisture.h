#pragma once
#ifndef moisture_h
#define moisture_h
#include <vector>

using namespace std;

const double rhoWood = 400.;		// Wood density (@TODO function.cpp uses 500 - need to make consistant)
const double timeStep = 60.;		// Timestep (1 minute)
const int MOISTURE_NODES = 11;   // Currently 11 nodes

class Moisture {
	private:
		
		double deltaX[MOISTURE_NODES];						// Node thickness
		double area[MOISTURE_NODES];							// Node area
		double volume[MOISTURE_NODES];						// Node volume
		double kappa1[MOISTURE_NODES], kappa2[MOISTURE_NODES];
		double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17;
		double xn67, xn68, xn69, xn7t, xn7o, xn7c, xn76, xn79, xn8t, xn8o, xn8c, xn86, xn87, xn89;
		double xn9t, xn9o, xn9c, xn96, xn97, xn98, xn910, xn109, xn10t, xn10o;
		vector< vector<double> > A;
		double PWOld[MOISTURE_NODES];							// previous time step vapor pressure (Pa)
      double PWInit[MOISTURE_NODES];						// initial vapor pressure (was B() in BASIC code) (Pa)
		double temperature[MOISTURE_NODES];					// Node temperature (deg K)
		double tempOld[MOISTURE_NODES];						// previous time step temperature (deg K)
		double haHouse;											// haHouse is the moisture transport coefficient of the house mass (kg/s)
		double massWHouse;										// active mass of moisture in the house (kg)

		void cond_bal(int pressure);
		double calc_kappa_1(int pressure, double temp, double mc, double volume);
		double calc_kappa_2(double mc, double volume);
      double mc_cubic(double pw, int pressure, double temp);
      double calc_vapor_pressure(double mc, double temp, int pressure);
		
	public:
		double moistureContent[MOISTURE_NODES];			// Node moisture content (%)
		vector <double> PW;										// Node vapor pressure (Pa). vector so it can be passed to gauss()
		double mTotal[MOISTURE_NODES];						// Node mass of condensed water (kg)

		Moisture(double atticVolume, double retVolume, double supVolume, double houseVolume, double floorArea,
					double atticArea, double roofPitch, double mcInit=0.15);
		void mass_cond_bal(double* node_temps, double tempOut, double RHOut,
               double airDensityOut, double airDensityAttic, double airDensityHouse, double airDensitySup, double airDensityRet,
               int pressure, double hU0, double hU1, double hU2,
               double mAtticIn, double mAtticOut, double mCeiling, double mHouseIn, double mHouseOut,
               double mAH, double mRetAHoff, double mRetLeak, double mRetReg, double mRetOut,
               double mSupAHoff, double mSupLeak, double mSupReg, double latcap, double dhMoistRemv, double latload);
};

void print_matrix(vector< vector<double> > A);

#endif