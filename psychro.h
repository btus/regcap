#pragma once
#ifndef psychro_h
#define psychro_h

//using namespace std;

double saturation_vapor_pressure(double temp);
double calcHumidityRatio(double dewpoint, double pressure);
double KtoF(double degK);
double calcHfgAir(double temp);

#endif