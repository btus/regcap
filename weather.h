#pragma once
#ifndef weather_h
#define weather_h
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

//using namespace std;

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
	double inTdb;
	double inHR;
	double supplyTdb;
	double qAH;
	int julianMinute;
};

// Modified to use double instead of string from:
// http://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
class CSVRow
{
    public:
        double const& operator[](std::size_t index) const
        {
            return m_data[index];
        }
        std::size_t size() const
        {
            return m_data.size();
        }
        void readNextRow(std::istream& str)
        {
            std::string         line;
            std::getline(str,line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream,cell,','))
            {
                m_data.push_back(atof(cell.c_str()));
            }
        }
    private:
        std::vector<double>    m_data;
};

std::istream& operator>>(std::istream& str,CSVRow& data);
double interpolate(double begin, double end, int step);
weatherData interpolateWind(weatherData begin, weatherData end, int step);
weatherData interpWeather(weatherData begin, weatherData end, int minute);
weatherData readTMY3(std::ifstream& file);
weatherData readEPW(std::ifstream& file);
weatherData readOneMinuteWeather(std::ifstream& file);
weatherData readRealWeather(std::ifstream& file);
double surfaceInsolation(double direct, double diffuse, double beta, double sigma, double gamma);

#endif

