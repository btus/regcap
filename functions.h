#pragma once
#ifndef functions_h
#define functions_h

//#include <iostream>
//#include <string>
//#include <vector>

//using namespace std;

// additional data types
struct atticVent_struct {
	int wall;
	double h;
	double A;
	double n;
	double m;
	double dP;
};

struct soffit_struct {
	double h;
	double m;
	double dP;
};

struct winDoor_struct {
	int wall;
	double Bottom;
	double Top;
	double High;
	double Wide;
	double m;
	double mIN;
	double mOUT;
	double dPtop;
	double dPbottom;
};

struct fan_struct {
	double power;
	double q;
	double m;
	double on; // Set on=1 in main program based on .oper
	double oper;
};

struct pipe_struct {
	int wall;
	double h;
	double A;
	double n;
	double m;
	double dP;
	double Swf;
	double Swoff;
};

struct flue_struct {
	double flueC;
	double flueHeight;
	double flueTemp;
};

//ASHRAE 62.2-2016 Infiltration and Relative Dose Functions

void sub_infiltrationModel (

	double& envC, //Envelope leakage coefficient, m3/s/Pa^n
	double& envPressureExp, //Envelope pressure exponent
	double& G, //Wind speed multiplier
	double& s, //Shelter Factor
	double& Cs, //Stack coefficient
	double& Cw, //Wind coefficient
	double& windSpeed, //Wind speed, corrected for site conditions in main.cpp, m/s
	double& dryBulb, //Outside temp, K
	double& ventSum, //Sum of the larger of the mechanical inflows and outflows, ACH
	double& houseVolume, //House volume, m3
	double& windSpeedCorrection, //Correction factor used in main.cpp for windspeed
	double& Q_wind, //Wind driven airflow, L/s
	double& Q_stack, //Stack pressure driven airflow, L/s
	double& Q_infiltration, //Total infiltration airflow, combined wind and stack airflows, L/s
	double& Q_total, //Total airflow combined infiltration and mechanical, L/s
	double wInfil,
	int InfCalc
	
	);
	
double sub_relativeExposure (

	double& Aeq, //Qtot calculated according to 62.2-2016 without infiltration factor, ACH. 
	double& Q_total, //Total airflow combined infiltration and mechanical, L/s
	double& relExp_old, //Relative exposure from the prior time-step. 
	double dtau, //Simulation timestep in seconds (60). 
	double& houseVolume //House volume, m3
	//double& relExp //relative exposure, per 62.2-2016
	
	); 	
	
double sub_moldIndex(

	//Variables passed back and forth with main.cpp
	int SensitivityClass, //Material sensitivity class, determined in Table 6.1.1. 0 = VerySensitive, 1 = Sensitive.
	double MoldIndex_old, //MoldIndex from the prior hour.
	double SurfTemp_K, //Material surface temperature, K.
	double SurfPw,	//Material surface partial vapor pressure, Pa
	int& Time_decl //MoldIndex decline time, hr

	);	
	
double sub_Pollutant (

	double outdoorConc, 
	double indoorConc, 
	double indoorSource, 
	double houseVolume, 
	double qHouse, 
	double qDeposition,
	double penetrationFactor, 
	double qAH, 
	double AHflag,
	double filterEfficiency 
	);


// Additional functions

void sub_heat ( 
	weatherData weather, 
	double& mCeiling, 
	double& AL4, 
	double& ssolrad, 
	double& nsolrad, 
	double* tempOld, 
	double& atticVolume, 
	double& houseVolume, 
	double* b,
	double& floorArea, 
	double& roofPitch, 
	double& ductLocation, 
	double& mSupReg, 
	double& mRetReg, 
	double& mRetLeak, 
	double& mSupLeak, 
	double& mAH, 
	double& supRval, 
	double& retRval, 
	double& supDiameter, 
	double& retDiameter, 
	double& supThickness, 
	double& retThickness, 
	double& supVel, 
	double& retVel, 
	double& uaSolAir,
	double& uaTOut, 
	double& matticenvin, 
	double& matticenvout, 
	double& mHouseIN, 
	double& mHouseOUT, 
	double& planArea, 
	double& mSupAHoff, 
	double& mRetAHoff, 
	double& solgain, 
	double& tsolair, 
	double& mFanCycler, 
	double& roofPeakHeight, 
	double& retLength,
	double& supLength,
	int& roofType,
	double roofExtRval,
	double roofIntRval,
	double ceilRval,
	double gableEndRval,
	int& AHflag, 
	double& mERV_AH,
	double& ERV_SRE,
	double& mHRV,
	double& HRV_ASE,
	double& mHRV_AH,
	double& capacityc,
	double& capacityh,
	double& evapcap,
	double& internalGains,
	//int bsize,
	double& airDensityIN,
	double& airDensityOUT,
	double& airDensityATTIC,
	double& airDensitySUP,
	double& airDensityRET,
	int& numStories,
	double& storyHeight,
	double dhSensibleGain,
	double& H2,
	double& H4,
	double& H6,
	double bulkArea,
	double sheathArea
);

void sub_houseLeak ( 
	int& AHflag,
	int& leakIterations, 
	double& U, 
	int& windAngle, 
	double& tempHouse, 
	double& tempAttic, 
	double& tempOut, 
	double& envC, 
	double& n, 
	double& eaveHeight, 
	double leakFracCeil, 
	double leakFracFloor,
	double leakFracWall, 
	int& numFlues, 
	flue_struct* flue, 
	double* wallFraction, 
	double* floorFraction, 
	double* Sw, 
	double& flueShelterFactor, 
	int& numWinDoor, 
	winDoor_struct* winDoor, 
	int& numFans, 
	fan_struct* fan, 
	int& numPipes, 
	pipe_struct* Pipe, 
	double& mIN, 
	double& mOUT, 
	double& Pint, 
	double& mFlue, 
	double& mCeiling, 
	double* mFloor, 
	double& atticC, 
	double& dPflue, 
	int& Crawl, 
	double& Hfloor, 
	bool rowHouse, 
	double* soffitFraction, 
	double& Patticint, 
	double* wallCp, 
	//double& airDensityRef, 
	//double& airTempRef, 
	double& mSupReg, 
	double& mAH, 
	double& mRetLeak, 
	double& mSupLeak, 
	double& mRetReg, 
	double& mHouseIN, 
	double& mHouseOUT, 
	double& supC, 
	double& supn, 
	double& retC, 
	double& retn, 
	double& mSupAHoff, 
	double& mRetAHoff, 
	//double& Aeq,
	double& airDensityIN,
	double& airDensityOUT,
	double& airDensityATTIC,
	double& houseVolume,
	double& windPressureExp
);

void sub_atticLeak ( 
	int& leakIterations, 
	double& U, 
	int& windAngle, 
	double& tempHouse, 
	double& tempOut, 
	double& tempAttic, 
	double& atticC, 
	double& atticPressureExp, 
	double& eaveHeight, 
	double& roofPeakHeight, 
	double& flueShelterFactor, 
	double* Sw, 
	int& numAtticVents, 
	atticVent_struct* atticVent, 
	soffit_struct* soffit, 
	double& mAtticIN, 
	double& mAtticOUT, 
	double& Patticint, 
	double& mCeiling, 
	bool rowHouse, 
	double* soffitFraction, 
	double& roofPitch, 
	bool roofPeakPerpendicular, 
	int& numAtticFans, 
	fan_struct* atticFan, 
	//double& airDensityRef, 
	//double& airTempRef, 
	double& mSupReg, 
	double& mRetLeak, 
	double& mSupLeak, 
	double& matticenvin, 
	double& matticenvout, 
	double& mSupAHoff, 
	double& mRetAHoff,
	double& airDensityIN,
	double& airDensityOUT,
	double& airDensityATTIC
	);

void sub_filterLoading (
	int& MERV,
	int& loadingRate,
	int& BPMflag,
	double& A_qAH_heat,
	double& A_qAH_cool,
	double& A_wAH_heat, 
	double& A_wAH_cool, 
	double& A_DL, 
	double& k_qAH,
	double& k_wAH,
	//double& k_qAH_heat, 
	//double& k_qAH_cool, 
	//double& k_wAH_heat, 
	//double& k_wAH_cool, 
	double& k_DL,
	double& qAH_heat0, 
	double& qAH_cool0,
	double& qAH_low
	);

double saturationVaporPressure (double temp);

#endif
