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
	result.humidityRatio = interpolate(begin.humidityRatio, end.humidityRatio, minute);
	
	return result;
}

weatherData readTMY3(ifstream& file) {
	CSVRow row;
	weatherData result;
	file >> row;
	if(row.size() >= 68) {							// don't set if at EOF
		result.directNormal = row[7];
		result.globalHorizontal = row[4];
		result.dryBulb = row[31] + C_TO_K;		// convert from C to K
		result.dewPoint = row[34] + C_TO_K;		// convert from C to K
		result.relativeHumidity = row[37];
		result.windSpeed = row[46];
		result.windDirection = row[43];
		result.pressure = row[40] * 100;			// convert from mbar to Pa
		result.skyCover = row[25] / 10;			// convert to decimal fraction
		result.humidityRatio = calcHumidityRatio(saturationVaporPressure(result.dewPoint), result.pressure);
	}
	return result;
}	

weatherData readEPW(ifstream& file) {
	CSVRow row;
	weatherData result;
	file >> row;
	if(row.size() >= 35) {							// hourly entry rows each have 35 columns. 10 is the header row. 
		result.directNormal = row[14];
		result.globalHorizontal = row[13];
		result.dryBulb = row[6] + C_TO_K;		// convert from C to K
		result.dewPoint = row[7] + C_TO_K;		// convert from C to K
		result.relativeHumidity = row[8];
		result.windSpeed = row[21];
		result.windDirection = row[20];
		result.pressure = row[9];			
		result.skyCover = row[22] / 10.;			// convert to decimal fraction
		result.humidityRatio = calcHumidityRatio(saturationVaporPressure(result.dewPoint), result.pressure);
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

/*
* surfaceInsolation()
*
* Calculates solar irradiance incident on a surface from direct beam, diffuse, and reflected insolation
* From ASHRAE HOF (2009 SI) chapter 14
* @param direct - direct beam insolation (W/m2)
* @param diffuse - total diffuse insolation (W/m2)
* @param beta - solar altitude angle (radians)
* @param sigma - tilt angle of surface (radians)
* @param gamma - surface-solar azimuth angle (phi - psi) (radians)
* @return surface solar insolation (W/m2)
*/

double surfaceInsolation(double direct, double diffuse, double beta, double sigma, double gamma) {
	double cosTheta;		// cosine of angle of incidence, eqn (22)
	double Y;				// ratio of vertical surface to clear-sky diffuse, eqn (28)
	double Eb = 0;			// direct beam component, eqn (26)
	double Ed;				// diffuse component, 2001 IP, ch29, table 14
	double Er;				// reflected component, eqn (31)
	
	if(sigma < M_PI / 2) {	// roofs
		cosTheta = cos(beta) * cos(gamma) * sin(sigma) + sin(beta) * cos(sigma);
		if(acos(cosTheta) < M_PI / 2) {
			Eb = cosTheta * direct;
		}
		// this is the diffuse component from the 2001 IP edition, pg 29.14.
		// Equation (29) in 2009 results in multipliers greater than 1 which does not make sense.
		Ed = diffuse * (1 + cos(sigma)) / 2;
		Er = (direct * sin(beta) + diffuse) * 0.2 * (1 - cos(sigma)) / 2;
		//cout << "Roof: Eb=" << Eb << " Ed=" << Ed << " Er=" << Er << endl;
	}
	else {						// walls
		cosTheta = cos(beta) * cos(gamma);
		Y = max(0.45,.55 + .437 * cosTheta + .313 * pow(cosTheta, 2));
		if(acos(cosTheta) < M_PI / 2) {
			Eb = cosTheta * direct;
		}
		Ed = diffuse * Y;
		Er = (direct * sin(beta) + diffuse) * 0.2 * 0.5;
		//cout << "Wall: Eb=" << Eb << " Ed=" << Ed << " Er=" << Er << endl;
	}
	return Eb + Ed + Er;
} 