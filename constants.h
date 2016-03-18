#pragma once
#ifndef constant_h
#define constant_h

// ============================== CONSTANTS ===============================================
const double C_TO_K = 273.15;
const double airDensityRef = 1.20411;	// Reference air density at 20 deg C at sea level [kg/m3]
const double g = 9.81;						// Acceleration due to gravity (m/s^2)
const int ArraySize = 16;
const double airTempRef = C_TO_K + 20;	// Reference room temp [K] = 20 deg C
const double SIGMA = 5.6704E-08;			// STEFAN-BOLTZMANN CONST (W/m^2/K^4)
const double CpAir = 1005.7;				// specific heat of air [j/kg K]

#endif