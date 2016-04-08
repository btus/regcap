#pragma once
#ifndef psychro_h
#define psychro_h

//using namespace std;

double saturationVaporPressure(double temp);
double calcHumidityRatio(double dewpoint, double pressure);
double KtoF(double degK);

#endif