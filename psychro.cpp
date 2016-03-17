#include "constants.h"
#include "psychro.h"

#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

double saturationVaporPressure(double temp) {
	//Coefficients for saturation vapor pressure over ice -100 to 0C. ASHRAE HoF.
	const double C1 = -5.6745359E+03;
	const double C2 = 6.3925247E+00;
	const double C3 = -9.6778430E-03;
	const double C4 = 6.2215701E-07;
	const double C5 = 2.0747825E-09;
	const double C6 = -9.4840240E-13;
	const double C7 = 4.1635019E+00;

	//Coefficients for saturation vapor pressure over liquid water 0 to 200C. ASHRAE HoF.
	const double C8 = -5.8002206E+03;
	const double C9 = 1.3914993E+00;
	const double C10 = -4.8640239E-02;
	const double C11 = 4.1764768E-05;
	const double C12 = -1.4452093E-08;
	const double C13 = 6.5459673E+00;

	//Calculate Saturation Vapor Pressure, Equations 5 and 6 in 2009 ASHRAE HoF 1.2
	if(temp <= C_TO_K){
		return exp((C1/temp)+(C2)+(C3*temp)+(C4*pow(temp, 2))+(C5*pow(temp, 3))+(C6*pow(temp, 4))+(C7*log(temp)));
	} else{
		return exp((C8/temp)+(C9)+(C10*temp)+(C11*pow(temp, 2))+(C12*pow(temp, 3))+(C13*log(temp)));
	}
}
