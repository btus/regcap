#pragma once
#ifndef psychro_h
#define psychro_h

//using namespace std;

double saturationVaporPressure(double temp);
double calcHumidityRatio(double pw, double pressure);
double KtoF(double degK);
double calcHfgAir(double temp);

#endif