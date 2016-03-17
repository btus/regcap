#include <iostream>
#include <fstream>
#include "weather.h"
#include "psychro.h"
#include "constants.h"

#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

weatherData readOneMinuteWeather(ifstream& file) {
	int day;
	weatherData result;
	double wd;
	file >> day
		>> result.directNormal
		>> result.globalHorizontal
		>> result.dryBulb
		>> result.humidityRatio
		>> result.windSpeed
		>> wd
		>> result.pressure
		>> result.skyCover;

	result.dryBulb += C_TO_K;					// Convert to deg K
	result.windDirection = int(wd);			// Wind direction is float in 1 minute files
	result.pressure *= 1000;					// Convert reference pressure to [Pa]
	result.skyCover /= 10;						// Converting cloud cover index to decimal fraction

	return result;
}
