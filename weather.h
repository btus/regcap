#pragma once
#ifndef weather_h
#define weather_h

using namespace std;

struct weatherData {
	int directNormal;				// Direct normal irradiance (Wh/m2)
	int globalHorizontal;		// Global horizontal irradiance (Wh/m2)
	double dryBulb;				// Dry-bulb temperature (deg K)
	double dewPoint;				// Dew-point temperature (deg K)
	double relativeHumidity;	// Relative humidity (%)
	double humidityRatio;		// Humidity Ratio (g/g)
	double windSpeed;				// Wind speed (m/s)
	int windDirection;			// wind direction (degrees from North)
	int pressure;					// Station pressure (Pa)
	double skyCover;				// Total sky cover (fraction)
};

weatherData readOneMinuteWeather(ifstream& file);

#endif

