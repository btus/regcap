#pragma once
#ifndef wood_h
#define wood_h

//using namespace std;

class WoodMoisture {
	private:
		static const int rhoWood = 400;		// Wood density (@TODO function.cpp uses 500 - need to make consistant)
		static const int timeStep = 3600;	// Timestep (1 hour)
		
		double deltaX[7];							// Node thickness
		double area[7];							// Node area
		double volume[7];							// Node volume
		double kappa1[7], kappa2[7];
		double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17;
		double A[7][7] = {0};
		double PWOld[7];
		double tempOld[7];

		void cond_bal(double pressure);
		double calc_kappa_1(double pressure, double temp, double mc, double volume);
		double calc_kappa_2(double mc, double volume);
        double mc_cubic(double pw, double pressure, double temp);
        double calc_vapor_pressure(double mc, double temp, double pressure);
		
	public:
		double temperature[7];
		double moistureContent[7];
		double PW[7];
		double mTotal[7] = {0};

		WoodMoisture(double woodThick, double atticVolume, double atticArea, double roofPitch, double tempInit, double mcInit);
        void mass_cond_bal(double tempOut, double RHOut, double tempHouse, double RHHouse,
                             double airDensityOut, double airDensityAttic, double airDensityHouse,
                             double pressure, double hU0, double hU1, double hU2,
                             double mAttic, double mCeiling);
};


#endif