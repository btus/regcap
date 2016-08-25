#pragma once
#ifndef functions_h
#define functions_h

#include <iostream>
#include <string>
#include <vector>

using namespace std;

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

// Additional functions

void sub_heat ( 
	double& tempOut, 
	//double& airDensityRef, 
	//double& airTempRef, 
	double& mCeiling, 
	double& AL4, 
	double& windSpeed, 
	double& ssolrad, 
	double& nsolrad, 
	double* tempOld, 
	double& atticVolume, 
	double& houseVolume, 
	double& sc, 
	double* b,
	int& ERRCODE, 
	double& TSKY, 
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
	double& supArea, 
	double& retArea, 
	double& supThickness, 
	double& retThickness, 
	double& supVolume, 
	double& retVolume, 
	double& supCp, 
	double& retCp, 
	double& supVel, 
	double& retVel, 
	double& suprho, 
	double& retrho, 
	int& pRef, 
	double& HROUT, 
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
	double& h, 
	//double& Whouse,
	double& retLength,
	double& supLength,
	int& roofType,
	double& M1,
	double& M12,
	double& M15,
	double& M16,
	double& roofRval,
	double& rceil,
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
	double dhSensibleGain
);

void sub_moisture ( 
	double* HR, 
	double& Mw1, 
	double& Mw2, 
	double& Mw3, 
	double& Mw4, 
	double& Mw5, 
	double& matticenvout, 
	double& mCeiling, 
	double& mSupAHoff, 
	double& mRetAHoff, 
	double& matticenvin, 
	double& HROUT, 
	double& mSupLeak, 
	double& mAH, 
	double& mRetReg, 
	double& mRetLeak, 
	double& mSupReg, 
	double& latcap, 
	double& mHouseIN, 
	double& mHouseOUT, 
	double& latentLoad, 
	double& mFanCycler, 
	double& mHRV_AH,
	double& mERV_AH, 
	double& ERV_TRE,
	double& MWha,
	double& airDensityIN,
	double& airDensityOUT,
	double dhMoistureRemoved
);

void sub_houseLeak ( 
	int& AHflag,
	int& flag, 
	double& U, 
	int& windAngle, 
	double& tempHouse, 
	double& tempAttic, 
	double& tempOut, 
	double& envC, 
	double& n, 
	double& h, 
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
	double& Aeq,
	double& airDensityIN,
	double& airDensityOUT,
	double& airDensityATTIC,
	double& houseVolume,
	double& windPressureExp
);

void sub_atticLeak ( 
	int& flag, 
	double& U, 
	int& windAngle, 
	double& tempHouse, 
	double& tempOut, 
	double& tempAttic, 
	double& atticC, 
	double& atticPressureExp, 
	double& h, 
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
