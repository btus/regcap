#pragma once
#ifndef wood_h
#define wood_h

//using namespace std;

class WoodMoisture {
	private:
		static const int rhoWood = 400;		// Wood density (@TODO function.cpp uses 500 - need to make consistant)
		static const int timeStep = 3600;	// Timestep (1 hour)
		// humidity ratio constants (from Cleary 1985)
		static const double B3 =  15.8;		// was MCA
		static const double B4 = -0.0015;	// was MCB
		static const double B5 =  0.053;		// was MCC
		static const double B6 = -0.184;		// was MCD
		static const double B7 =  0.233;		// was MCE
		
		double deltaX[7];							// Node thickness
		double area[7];							// Node area
		double volume[7];							// Node volume
		double kappa1[7], kappa2[7];
		double x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17;
		double A[7][7] = {0};
		double PWOld[7];
		double tempOld[7];

		double cond_bal();
		double calc_kappa_1(double pressure, double temp, double mc);
		double calc_kappa_2(double mc);
		
	public:
		double temperature[7];
		double moistureContent[7];
		double PW[7];
		double mTotal[7] = {0};

		WoodMoisture();
		double mass_cond_bal();
};


#endif