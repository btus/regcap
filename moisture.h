#pragma once
#ifndef moisture_h
#define moisture_h
#include <vector>

//using namespace std;

const double timeStep = 60.;		// Timestep (1 minute)
const int MOISTURE_NODES = 13;   // Max number of nodes

class Moisture {
	private:
		int moisture_nodes;										// Number of moisture nodes
		double deltaX[MOISTURE_NODES];						// Node thickness
		double area[MOISTURE_NODES];							// Node area
		double volume[MOISTURE_NODES];						// Node volume
		double density[6];										// Wood node density
		double kappa1[MOISTURE_NODES], kappa2[MOISTURE_NODES];
		double x60, x30, x06, x03, x61, x41, x16, x14, x62, x52, x26, x25, x66, x6out;
		double x67, x68, x69, x011, x110, x112, x121, x611, x116, x612, x126;
		std::vector< std::vector<double> > A;
		double PWOld[MOISTURE_NODES];							// previous time step vapor pressure (Pa)
      double PWInit[MOISTURE_NODES];						// initial vapor pressure (was B() in BASIC code) (Pa)
		double tempOld[MOISTURE_NODES];						// previous time step temperature (deg K)
		double haHouse;											// haHouse is the moisture transport coefficient of the house mass (kg/s)
		double massWHouse;										// active mass of moisture in the house (kg)
		double roofInsulRatio;									// ratio of exterior insulation U-val to sheathing U-val
		
		void cond_bal(int pressure);
		double calc_kappa_1(int pressure, double temp, double mc, double mass);
		double calc_kappa_2(double mc, double mass);
      double mc_cubic(double pw, int pressure, double temp);
      double calc_vapor_pressure(double mc, double temp, int pressure);
		double calc_inter_temp(double temp1, double temp2, double insRatio);
		
	public:
		double temperature[MOISTURE_NODES];					// Node temperature (deg K)
		double moistureContent[MOISTURE_NODES];			// Node moisture content (%)
		std::vector <double> PW;										// Node vapor pressure (Pa). vector so it can be passed to gauss()
		double mTotal[MOISTURE_NODES];						// Node mass of condensed water (kg)
		int saturated_minutes[MOISTURE_NODES];				// Number of minutes node vapor pressure is above saturation
		int total_in_iter, total_out_iter;				// Total number of inner and outer iterations

		Moisture(double atticVolume, double retDiameter, double retLength, double supDiameter, double supLength, double houseVolume,
					 double floorArea, double sheathArea, double bulkArea, double roofInsThick, double roofExtRval, double mcInit=0.15);
		void mass_cond_bal(double* node_temps, weatherData weather,
               double airDensityOut, double airDensityAttic, double airDensityHouse, double airDensitySup, double airDensityRet,
               double hU0, double hU1, double hU2,
               double mAtticIn, double mAtticOut, double mCeiling, double mHouseIn, double mHouseOut,
               double mAH, double mRetAHoff, double mRetLeak, double mRetReg, double mRetOut, double mErvHouse,
               double mSupAHoff, double mSupLeak, double mSupReg, double latcap, double dhMoistRemv, double latload);
};

void print_matrix(std::vector< std::vector<double> > A);

#endif