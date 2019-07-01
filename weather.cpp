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

/*
 * Weather - Weather file class constructor
 * @param fileName - name of weather file to open
 */
Weather::Weather(int terrain, double eaveHeight) {
	// Terrain where the house is located (for wind shelter etc.) See ASHRAE Fundamentals 2009 F24.3
	double layerThickness;		// Atmospheric boundary layer thickness [m]

	switch (terrain) {
	case 1:		// Large city centres, at least 50% of buildings are higher than 25m
		windPressureExp = 0.33;
		layerThickness = 460;
		break;
	case 2:		// Urban and suburban areas, wooded areas
		windPressureExp = 0.22;
		layerThickness = 370;
		break;
	case 3:		// Open terrain with scattered obsructions e.g. meteorlogical stations
		windPressureExp = 0.14;
		layerThickness = 270;
		break;
	case 4:		// Completely flat e.g. open water
		windPressureExp = 0.10;
		layerThickness = 210;
		break;
	}

	// Wind speed correction to adjust met wind speed to that at building eaves height
	const double metExponent = 0.14;		// Power law exponent of the wind speed profile at the met station
	const double metThickness = 270;		// Atmospheric boundary layer thickness at met station [m]
	const double metHeight = 10;			// Height of met station wind measurements [m]
	windSpeedCorrection = pow((metThickness / metHeight), metExponent) * pow((eaveHeight / layerThickness), windPressureExp);
	//double windSpeedCorrection = pow((80/ 10), .15) * pow((h / 80), .3);		// The old wind speed correction
	}

void Weather::open(string fileName) {
	weatherFile.open(fileName);
	if(!weatherFile) { 
		throw fileName; 
	}

	// Determine weather file type and read in header
	CSVRow row;					// row of data to read from TMY3 csv

	weatherFile >> row;
	if(row.size() > 2 && row.size() < 10) {			//TMY3, was >2
		siteID = row[0];
		//string siteName = row[1];
		//string State = row[2];
		timeZone = row[3];
		latitude = row[4];
		longitude = row[5];
		elevation = row[6];
		weatherFile >> row;					// drop header
		begin = readTMY3();	// duplicate first hour for interpolation of 0->1
		end = begin;
		type = 1;
	}
	
	else if(row.size() == 10) {			//EnergyPlus weather file EPW
		siteID = row[5];
		timeZone = row[8];
		latitude = row[6];
		longitude = row[7];
		elevation = row[9];
		for(int i = 0; i < 7; ++i) {	// drop rest of header lines
			weatherFile >> row;
			}
		begin = readEPW();	// duplicate first hour for interpolation of 0->1
		end = begin;
		type = 2;
	}
	
	else {							// 1 minute data
		weatherFile.close();
		weatherFile.open(fileName);
		weatherFile >> latitude >> longitude >> timeZone >> elevation;
		siteID = 0;
		type = 0;
	}

}


double Weather::interpolate(double b, double e, int step) {
	return b + (e - b) / 60 * step;
}

weatherData Weather::interpolateWind(int step) {
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

weatherData Weather::interpWeather(int minute) {
	weatherData result;
	
	result = interpolateWind(minute);
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

weatherData Weather::readTMY3() {
	CSVRow row;
	weatherData result;
	weatherFile >> row;
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

weatherData Weather::readEPW() {
	CSVRow row;
	weatherData result;
	weatherFile >> row;
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

weatherData Weather::readOneMinuteWeather() {
	int day;
	weatherData result;
	double wd;
	weatherFile >> day
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

weatherData Weather::readMinute(int minute) {
weatherData current;

	// Read in or interpolate weather data
	if(type == 0) {  // minute
		current = readOneMinuteWeather();
		}
	else {
		current = interpWeather(minute);
		}
	if(current.windSpeed < 1)							// Minimum wind velocity allowed is 1 m/s to account for non-zero start up velocity of anenometers
		current.windSpeed = 1;							// Wind speed is never zero
	current.windSpeedLocal = current.windSpeed * windSpeedCorrection;		// Correct met wind speed to speed at building eaves height [m/s]
	return current;
	}

void Weather::nextHour() {
	begin = end;
	switch (type) {
	case 1:
		end = readTMY3();
		break;
	case 2:
		end = readEPW();
		break;
		}
	}

void Weather::close() {
	weatherFile.close();
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