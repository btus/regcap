#pragma once
#ifndef wood_h
#define wood_h
#include <vector>

using namespace std;

class WoodMoisture {
	private:
		static const int rhoWood = 400;		// Wood density (@TODO function.cpp uses 500 - need to make consistant)
		static const int timeStep = 60;		// Timestep (1 minute)
		static const int MOISTURE_NODES = 7; // Currently 7 nodes
		
		double deltaX[MOISTURE_NODES];						// Node thickness
		double area[MOISTURE_NODES];							// Node area
		double volume[MOISTURE_NODES];						// Node volume
		double kappa1[MOISTURE_NODES], kappa2[MOISTURE_NODES];
		double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17;
		vector< vector<double> > A;
		double PWOld[MOISTURE_NODES];							// previous time step vapor pressure (Pa)
      double PWInit[MOISTURE_NODES];						// initial vapor pressure (was B() in BASIC code) (Pa)
		double temperature[MOISTURE_NODES];					// Node temperature (deg K)
		double tempOld[MOISTURE_NODES];						// previous time step temperature (deg K)

		void cond_bal(int pressure);
		double calc_kappa_1(int pressure, double temp, double mc, double volume);
		double calc_kappa_2(double mc, double volume);
      double mc_cubic(double pw, int pressure, double temp);
      double calc_vapor_pressure(double mc, double temp, int pressure);
		
	public:
		double moistureContent[MOISTURE_NODES];			// Node moisture content (%)
		vector <double> PW;										// Node vapor pressure (Pa). vector so it can be passed to gauss()
		double mTotal[MOISTURE_NODES];						// Node mass of condensed water (kg)

		WoodMoisture(double atticVolume, double atticArea, double roofPitch, double mcInit);
      void mass_cond_bal(double* node_temps, double tempOut, double RHOut, double RHHouse,
                             double airDensityOut, double airDensityAttic, double airDensityHouse,
                             int pressure, double hU0, double hU1, double hU2,
                             double mAttic, double mCeiling);
};


#endif