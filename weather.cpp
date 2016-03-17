#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "weather.h"
#include "psychro.h"
#include "constants.h"

#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

std::istream& operator>>(std::istream& str,CSVRow& data)
{
    data.readNextRow(str);
    return str;
}   

double interpolate(double begin, double end, int step) {
	return begin + (end - begin) / 60 * step;
}

weatherData interpolateWind(weatherData begin, weatherData end, int step) {
	weatherData result;
	
	if(begin.windSpeed == 0) {
		result.windDirection = end.windDirection;
		result.windSpeed = interpolate(0, end.windSpeed, step);
	}
	else if(end.windSpeed == 0) {
		result.windDirection = begin.windDirection;
		result.windSpeed = interpolate(begin.windSpeed, 0, step);
	}
	else {
		double ns1 = cos(begin.windDirection * M_PI / 180) * begin.windSpeed;
		double ew1 = sin(begin.windDirection * M_PI / 180) * begin.windSpeed;
		double ns2 = cos(end.windDirection * M_PI / 180) * end.windSpeed;
		double ew2 = sin(end.windDirection * M_PI / 180) * end.windSpeed;
		double ns = interpolate(ns1, ns2, step);
		double ew = interpolate(ew1, ew2, step);
	
		result.windSpeed = pow((pow(ns, 2) + pow(ew, 2)), 0.5);
		result.windDirection = atan2(ew, ns) * 180 / M_PI;
		if(result.windDirection < 0)
			result.windDirection += 360;
	}
	
	return result;
}

weatherData interpWeather(weatherData begin, weatherData end, int minute) {
	weatherData result;
	
	result = interpolateWind(begin, end, minute);
	result.directNormal = interpolate(begin.directNormal, end.directNormal, minute);
	result.globalHorizontal = interpolate(begin.globalHorizontal, end.globalHorizontal, minute);
	result.dryBulb = interpolate(begin.dryBulb, end.dryBulb, minute);
	result.dewPoint = interpolate(begin.dewPoint, end.dewPoint, minute);
	result.relativeHumidity = interpolate(begin.relativeHumidity, end.relativeHumidity, minute);
	result.pressure = interpolate(begin.pressure, end.pressure, minute);
	result.skyCover = interpolate(begin.skyCover, end.skyCover, minute);
	//result.humidityRatio = calcHumidityRatio(result.dewPoint, result.pressure);
	result.humidityRatio = interpolate(begin.humidityRatio, end.humidityRatio, minute);
	
	return result;
}

weatherData readTMY3(ifstream& file) {
	CSVRow row;
	weatherData result;
	file >> row;
	if(row.size() == 71) {
		result.directNormal = row[7];
		result.globalHorizontal = row[4];
		result.dryBulb = row[31] + C_TO_K;		// convert from C to K
		result.dewPoint = row[34] + C_TO_K;		// convert from C to K
		result.relativeHumidity = row[37];
		result.windSpeed = row[46];
		result.windDirection = row[43];
		result.pressure = row[40] * 100;			// convert from mbar to Pa
		result.skyCover = row[25] / 10;			// convert to decimal fraction
		result.humidityRatio = calcHumidityRatio(result.dewPoint, result.pressure);
	}
	return result;
}	

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
