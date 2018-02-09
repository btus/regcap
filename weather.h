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

class Weather {
	private:
		weatherData begin;		// first hour weather data
		weatherData end;			// second hour weather data
		ifstream weatherFile;	// weather file

		double interpolate(double b, double e, int step);
		weatherData interpolateWind(int step);
		weatherData readTMY3();
		weatherData readEPW();
		weatherData readOneMinuteWeather();
		weatherData interpWeather(int minute);
		
	public:
		int type;					// type of weather file (0 - 1 minute, 1 - TMY3 , 2 - EPW)
		double siteID;
		double latitude;
		double longitude;
		int timeZone;
		double elevation;
		
		Weather();
		void open(string fileName);
		weatherData readMinute(int minute);
		void nextHour();
		void close();
};

double surfaceInsolation(double direct, double diffuse, double beta, double sigma, double gamma);
std::istream& operator>>(std::istream& str,CSVRow& data);

#endif

