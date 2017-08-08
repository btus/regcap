#include "functions.h"
#include "psychro.h"
#include "constants.h"
#include "gauss.h"
#include <iomanip> // RAD: so far used only for setprecission() in cmd output
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif

using namespace std;

// ============================= FUNCTIONS ==============================================================
void f_CpTheta(double CP[4][4], int& windAngle, double* wallCp);

void f_flueFlow(double& tempHouse, double& flueShelterFactor, double& dPwind, double& dPtemp, double& eaveHeight, double& Pint, int& numFlues, flue_struct* flue, double& mFlue,
	double& airDensityOUT, double& airDensityIN, double& dPflue, double& tempOut, double& houseVolume, double& windPressureExp);

void f_floorFlow3(double& Cfloor, double& Cpfloor, double& dPwind, double& Pint,
	double& n, double& mFloor, double& airDensityOUT, double& airDensityIN, double& Hfloor, double& dPtemp);

void f_ceilingFlow(int& AHflag, double& Patticint, double& h, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& mCeiling, double& atticC, double& airDensityATTIC,
	double& airDensityIN, double& tempAttic, double& tempHouse, double& tempOut, double& airDensityOUT,
	double& mSupAHoff, double& mRetAHoff, double& supC, double& supn, double& retC, double& retn, double CCeiling);

void f_wallFlow3(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& wallCp,
	double& n, double& Cwall, double& eaveHeight, double& Pint, double& dPtemp, double& dPwind,
	double& mWallIn, double& mWallOut, double& Hfloor);

void f_fanFlow(fan_struct& fan, double& airDensityOUT, double& airDensityIN);

void f_pipeFlow(double& airDensityOUT, double& airDensityIN, double& CP, double& dPwind, double& dPtemp,
	double& Pint, pipe_struct& Pipe, double& tempHouse, double& tempOut);

void f_winDoorFlow(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& eaveHeight,
	double& wallCp, double& n, double& Pint, double& dPtemp, double& dPwind, winDoor_struct& winDoor);

void f_roofCpTheta(double* Cproof, int& windAngle, double* Cppitch, double& roofPitch);

double f_neutralLevel(double dPtemp, double dPwind, double Pint, double Cpr, double h);

void f_roofFlow(double& tempAttic, double& tempOut, double& airDensityATTIC, double& airDensityOUT, double& Cpr,
	double& atticPressureExp, double& Croof, double& roofPeakHeight, double& Patticint, double& dPtemp, double& dPwind,
	double& mRoofIn, double& mRoofOut, double& eaveHeight);

void f_atticVentFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	atticVent_struct& atticVent, double& tempAttic, double& tempOut);

void f_soffitFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	soffit_struct& soffit, double& soffitFraction, double& atticC, double& atticPressureExp, double& tempAttic, double& tempOut);

void f_atticFanFlow(fan_struct& atticFan, double& airDensityOUT, double& airDensityATTIC);

double heatTranCoef(double tempi, double tempa, double velocity);

double radTranCoef(double emissivity, double tempi, double tempj, double shapeFactor, double areaRatio);



//ASHRAE 62.2-2016 Infiltration and Relative Dose Functions

void sub_infiltrationModel(

	double& envC, //Envelope leakage coefficient, m3/s/Pa^n
	double& envPressureExp, //Envelope pressure exponent
	double& G, //Wind speed multiplier. This assumes a terrain factor of 2 (urban/suburban) and does not modify for other terrains.  
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
	
	) 
	
	{
	
	double TMYwindSpeed; //Un-corrected wind speed, so we can adjust wind speed based on 62.2-2016, m/s.
	double U; //Wind speed adjusted according to 62.2-2016, m/s
	double phi; //Additivity coefficient
	double ventSum_flow; //Mechanical airflow, L/s
	double envC_LitersPerSec; //Envelope leakage coefficient, L/s/Pa^n
	
	TMYwindSpeed = windSpeed / windSpeedCorrection;
	
	U = G * TMYwindSpeed;
	
	envC_LitersPerSec = envC * 1000; // [L/s/Pa^n] 
	
	ventSum_flow = 1000 * ventSum * houseVolume / 3600; // Calculate ventSum in airflow (not ACH), L/s

	Q_wind = envC_LitersPerSec * Cw * pow( s * U, 2 * envPressureExp); //Calculate wind driven airflow, L/s
	
	Q_stack = envC_LitersPerSec * Cs * pow(abs((C_TO_K + 20) - dryBulb), envPressureExp); //Calculation stack airflow, assumes 68F (20C indoor temp), L/s.
	
	if(InfCalc == 0){ //Use annual infiltration estimate from 62.2-2016.
		Q_infiltration = wInfil;
	} else { //Use real-time infiltration calculation.
		Q_infiltration = sqrt(pow(Q_wind, 2) + pow(Q_stack, 2)); //Calculate combined infiltration airflow using quadrature, L/s.
	}
	
	phi = Q_infiltration / (Q_infiltration + ventSum_flow); // Calculate the additivity coefficient
	
	Q_total = ventSum_flow + phi * Q_infiltration; //Calculate total airflow combined mechanical and natural, L/s.
	
	}
	
double sub_relativeExposure(

	double& Aeq, //Qtot calculated according to 62.2-2016 without infiltration factor, ACH. 
	double& Q_total, //Total airflow combined infiltration and mechanical, L/s
	double& relExp_old, //Relative exposure from the prior time-step.
	double dtau, //Simulation time step in seconds (60).
	double& houseVolume //House volume, m3
	//double& relExp
	
	) 
	
	{

	double Aeq_to_m3s; //Aeq (hr-1) converted to m3/s.
	double Q_total_to_m3s; //Q_total (L/s) converted to m3/s.
	double relExp; //relative exposure
	
	Q_total_to_m3s = (Q_total / 1000.);
	Aeq_to_m3s = Aeq * houseVolume / 3600.;

	if(Q_total_to_m3s == 0) {
		relExp = relExp_old + ((Aeq_to_m3s * dtau) / houseVolume);
	} else {
		relExp = (Aeq_to_m3s / Q_total_to_m3s) + (relExp_old - (Aeq_to_m3s / Q_total_to_m3s)) * exp(-1 * Q_total_to_m3s * dtau / houseVolume);
	}
	
	return(relExp);
}

double sub_moldIndex(

	//Mold Index calculation according to ASNI/ASHRAE Standard 160-2009 Addendum X (January 2016), 1st public review period.

	// Mold Index should be calculated for the following locations:
	
	// Node Identification (NOTE: The C++ array indices are one less than the node number e.g. node 1 = 0, node 2 = 1 etc.)
	// Node 2 is the Inner North Sheathing, b[1] (surface temps in b[] array).
	// Node 4 is the Inner South Sheathing, b[3]
	// Node 6 is all of the Wood (joists, trusses, etc.) lumped together, b[5]
	// Node 9 is the Inner Gable Wall (both lumped together), b[8]
	// Source for the surface RH values? (1) use the vapor pressure values translated to RH, or (2) calculate local surface RH based on attic air RH and surface temps. 

	// 	From the detailed attic wood moisture model, attic.temperature[n] (surface temperature node) and attic.PW[n] (surface vapor pressure node), where n is:
	//  0 - surface of south sheathing (facing the attic) - corresponding thermal node is (4)
	//  1 - surface of north sheathing (facing the attic) - corresponding thermal node is (2)
	//  2 - surface of the bulk wood in the attic  - corresponding thermal node is (6)
	//So there's no node in the wood attic model for the inner gable wall? (Node 9 thermal)

	//Variables passed back and forth with main.cpp
	int SensitivityClass, //Material sensitivity class, determined in Table 6.1.1. 0 = VerySensitive, 1 = Sensitive.
	double MoldIndex_old, //MoldIndex from the prior hour.
	double SurfTemp_K, //Material surface temperature, K.
	double SurfPw,	//Material surface partial vapor pressure, Pa
	int& Time_decl //MoldIndex decline time, hr

	)
	
	{
	//Variables used only internally in the function
	double A;
	double B;
	double C;
	double W;
	double k1; //Mold growth intensity factor.
	double k2; //Mold index attenuation factor.
	double k3; //Mold index decline coefficient, default to 0.25 (proposed Standard had 0.1)
	double MoldIndex_delta;
	double RH_critical; //critical surface RH for mold growth, % (0-100)
	double MoldIndex_Max;
	double SurfRH; //Material surface relative humidity, % (0-100)
	double SurfTemp; //Material surface temperature, C
	double MoldIndex_curr; //MoldIndex for the current hour.
	
	//Converting surface node values to RH and C.
	SurfRH  = 100 * (SurfPw / saturationVaporPressure(SurfTemp_K)); //Calculating surface node relative humidity (%, 0-100).
	SurfTemp = SurfTemp_K - C_TO_K; //Convert surface node temperature from K to C.

	//Calculated once each hour.
	
	//Selecting model coefficients
	
	k3 = 0.25; //
	//Coefficient Values from Table 6.1.2.
	if(SensitivityClass == 0) { //Very sensitive material
		A = 1.;
		B = 7.;
		C = 2.;
		W = 0.;
		if(MoldIndex_old < 1){
			k1 = 1.;
		} else {
			k1 = 2.;
		}	
	} else if (SensitivityClass == 1) { //Sensitive material.
		A = 0.3;
		B = 6.;
		C = 1.;
		W = 1.;
		if(MoldIndex_old < 1){
			k1 = 0.578;
		} else {
			k1 = 0.386;
		}
	}
	
	//Calculation of critical surface RH. Equation 6-2a.
	
	if(SurfTemp <= 20) {
		//RH_critical calculated for 'Sensitive' and 'Very Sensitive' material classes, which includes untreated wood, planed wood, paper products and wood-based boards.
		RH_critical = -0.00267 * pow(SurfTemp, 3) + 0.160 * pow(SurfTemp, 2) - 3.13 * SurfTemp + 100;
	} else {
		RH_critical = 80;
	}
	
	//Calculation of max mold index (for given surface temp / RH) and coefficient k2.
	
	MoldIndex_Max = A + B * ((RH_critical - SurfRH) / (RH_critical - 100)) - C * pow(((RH_critical - SurfRH) / (RH_critical - 100)), 2); //Equation 6-6
	
	k2 = (1 - exp(2.3 * (MoldIndex_old - MoldIndex_Max))); //Equation 6-5. 
	if(k2 <= 0) {
		k2 = 0;
	}
		
	//Calculation of MoldIndex_delta.	
	
	if(SurfTemp <= 0 || SurfRH <= RH_critical) { //Period of declining MoldIndex.
	
		Time_decl += 1; //Increment the Time_decl counter. Change of conditions from favorable to unfavorable.
	
		if(Time_decl <= 6) {
			MoldIndex_delta = -0.00133 * k3;
		} else if (Time_decl > 6 && Time_decl <= 24) {
			MoldIndex_delta = 0;	
		} else if(Time_decl > 24) {
			MoldIndex_delta = -0.000667 * k3;
		}	
		
	} else if(SurfTemp > 0 && SurfRH > RH_critical) { //Period of increasing MoldIndex.
		Time_decl = 0;
		MoldIndex_delta = (k1 * k2)/(168 * exp(-0.68 * log(SurfTemp) - 13.9 * log(SurfRH) + 0.14 * W + 66.02)); //Equation 6-4a
	}
	
	MoldIndex_curr = MoldIndex_old + MoldIndex_delta;
	
	if(MoldIndex_curr < 0) { //Set minimum value to 0. 
		MoldIndex_curr = 0;
	}
	
	return MoldIndex_curr;
	
}	
	


double sub_Pollutant (

	double outdoorConc, //Outdoor pollutant concentration, ug/m3
	double indoorConc, //Indoor pollutant concentration, ug/m3
	double indoorSource, //Indoor pollutant mass emission rate, ug/s
	double houseVolume, //m3
	double qHouse, //House airflow, equivalent to what a tracer gas would measure, m3/s
	double qDeposition, //Indoor deposition flow (loss) rate, m3/s. Default to 0 for no deposition. 
	double penetrationFactor, //Penetration factor is the fraction of infiltrating outside mass that gets inside. Default to 1 for no losses. 
	double qAH, //Central HVAC system airflow, m3/s
	double AHflag, //Air handler flag indicating status.
	double filterEfficiency //Filtration efficiency in the central HVAC system. Default to 0 for no filter. 
	)
	{
		//new concentration = total mass in house + net change in mass due to air exchange and indoor sources. Divided by the house volume.
		if(AHflag != 0){ //central AHU is ON.
			return 1 / houseVolume * (indoorConc * houseVolume + 60 * (qHouse * (-indoorConc + penetrationFactor * outdoorConc) + indoorSource - qDeposition * indoorConc - qAH * filterEfficiency * indoorConc));
		} else {
			return 1 / houseVolume * (indoorConc * houseVolume + 60 * (qHouse * (-indoorConc + penetrationFactor * outdoorConc) + indoorSource - qDeposition * indoorConc));
		}
}
					

//**********************************
// Functions definitions...

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
	double& skyCover, 
	double* x,
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
	double& innerNorthH,
	double& innerSouthH,
	double& bulkH,
	double bulkArea,
	double sheathArea,
	double& roofInsulRatio
) {
	vector<double> b;
	vector< vector<double> > A;
	double toldcur[ATTIC_NODES], area[ATTIC_NODES], mass[ATTIC_NODES];
	double heatCap[ATTIC_NODES], uVal[ATTIC_NODES], htCoef[ATTIC_NODES];
	double viewFactor[ATTIC_NODES][ATTIC_NODES] = {};
	double rtCoef[ATTIC_NODES][ATTIC_NODES];
	double woodThickness;
	double denShingles, cpShingles, Rshingles;
	double kWood;
	double kAir;
	double muAir;
	double characteristicVelocity;
	double HI;
	double woodEmissivity, roofAbsorptivity, roofEmissivity;	
	double gndCoef2, gndCoef4, skyCoef2, skyCoef4;
	double FRS, FG;
	double TSKY, PW;
	double TGROUND;
	int heatIterations;
	int roofInNorth, roofInSouth, attic_nodes;

	// Node 0 is the Attic Air
	// Node 1 is the Inner North Sheathing
	// Node 2 is the Outer North Sheathing
	// Node 3 is the Inner South Sheathing
	// Node 4 is the Outer South Sheathing
	// Node 5 is all of the Wood (joists, trusses, etc.) lumped together
	// Node 6 is the Ceiling of the House
	// Node 7 is the Floor of the Attic
	// Node 8 is the Inner Gable Wall (both lumped together)
	// Node 9 is the Outer Gable Wall (both lumped together)
	// Node 10 is the Return Duct Outer Surface
	// Node 11 is the Return Duct Air
	// Node 12 is The Mass of the House
	// Node 13 is the Supply Duct Outer Surface
	// Node 14 is the Supply Duct Air
	// Node 15 is the House Air (all one zone)
	// Node 16 is the Inner North Roof Insulation
	// Node 17 is the Inner South Roof Insulation

	if(roofIntRval > 0) {   // If there is interior insulation at the roof add nodes 16&17
      roofInNorth = 16;
      roofInSouth = 17;
      attic_nodes = 18;
      }
   else {                  // No insulation so interior nodes are the sheathing surface (1&3)
      roofInNorth = 1;
      roofInSouth = 3;
      attic_nodes = 16;
   }
   // set size of equation vectors to number of nodes (A contains both)
   A.resize(attic_nodes, vector<double>(attic_nodes+1, 0));
   b.resize(attic_nodes, 0);

	woodThickness = .015;		// thickness of sheathing material (m)
	woodEmissivity = .9;       // emissivity of building materials
	switch(roofType) {
		case 1:			// asphalt shingles
			roofAbsorptivity = .92;
			roofEmissivity = .91;
			Rshingles = .078;										// ASHRAE Fundamentals 2011 pg 26.7
			denShingles = 1100 * 2 * .005;					// asphalt shingles (factor of two because they overlap)
			cpShingles = 1260;									// CP asphalt shingles
			break;
		case 2: 			// red clay tile - edited for ConSol to be light brown concrete
			roofAbsorptivity = .58; 											// .67
			roofEmissivity = .9;
			Rshingles = .5;
			denShingles = 50;										// kg/m2
			cpShingles = 880;										// CP for tile roof
			break;
		case 3:			// low coating clay tile
			roofAbsorptivity = .5;
			roofEmissivity = .9;
			Rshingles = .5;
			denShingles = 50;
			cpShingles = 880;										// CP for tile roof
			break;
		case 4:			// asphalt shingles  & white coating
			roofAbsorptivity = .15;
			roofEmissivity = .91;
			Rshingles = .078;										// ASHRAE Fundamentals 2011 pg 26.7
			denShingles = 1100 * 2 * .005;					// asphalt shingles (factor of two because they overlap)
			cpShingles = 1260;									// CP asphalt shingles
			break;
	}

	PW = HROUT * pRef / (.621945 + HROUT);						// water vapor partial pressure pg 1.9 ASHRAE fundamentals 2009
	PW = PW / 1000 / 3.38;											// CONVERT TO INCHES OF HG
	TSKY = tempOut * pow((.55 + .33 * sqrt(PW)), .25);		// TSKY DEPENDS ON PW

	// Surface Area of Nodes
	area[1] = sheathArea;
	area[3] = area[1];
	
	// the following are commented out for ConSOl because cement tile is flat and does not have increased surface area
	if(roofType == 2 || roofType == 3) {
        area[2] = 1.5 * area[1];	   // tile roof has more surface area for convection heat transfer
	} else {
		area[2] = area[1];																
	}

   area[4] = area[2];
	area[5] = bulkArea;
	area[6] = planArea;									         // Ceiling
	area[7] = area[6];											   // Attic floor
	area[8] = planArea / 2 * tan(roofPitch * M_PI / 180);	// Total endwall area (assumes 2:1 aspect ratio)
	area[9] = area[8];
	area[10] = retArea;
	area[11] = M_PI * retLength * retDiameter;
	area[13] = supArea;
	area[14] = M_PI * supLength * supDiameter;
	
	// surface area of inside of house minus the end walls, roof and ceiling
	// Currently assuming two stories with heights of 2.5m and 3.0m
	//area[15] = 3 * pow(floorArea, .5) * 2 + 2.5 * pow(floorArea, .5) * 2 + 2 * floorArea;
	//area[15] = numStories * storyHeight * pow(planArea, .5) * 2 + (2 * (numStories -1)) * planArea;
	area[15] = 11 * pow(floorArea,0.5) + 2 * floorArea;	// Empirically derived relationship
	area[12] = 6 * area[15];									//  Surface area of everything in the house
	area[16] = area[1];
	area[17] = area[3];
	
	// masses
	mass[0] = atticVolume * airDensityATTIC;					     // mass of attic air
	mass[1] = .5 * area[1] * densitySheathing * woodThickness; // 1/2 OF TOTAL
	mass[2] = mass[1] + denShingles * area[1];					  // OTHER 1/2 OUTSIDE SHEATHING, WOOD TCOND W/MMC
	mass[3] = mass[1];
	mass[4] = mass[2];
	mass[5] = area[5] * 0.013 * densityWood;                   // 13mm equivalent thickness for 2x4 trusses
	mass[6] = .5 * 4 * area[6];									     // MASS OF JOISTS DRYWALL AND INSULATION
	mass[7] = mass[6];
	mass[8] = .5 * densityWood * woodThickness * area[8];
	mass[9] = mass[8];
	mass[10] = retLength * M_PI * (retDiameter + retThickness) * retThickness * retrho;	// retArea * retrho * retThickness
	mass[11] = retVolume * airDensityRET;

	//  maybe not - 02/2004 need to increase house mass with furnishings and their area: say 5000kg furnishings
	mass[12] = (storyHeight * pow(floorArea, .5) * 4 * 2000 * .01 + planArea * .05 * 2000);		// mass of walls (5 cm effctive thickness) + mass of slab (aso 5 cm thick)
	mass[13] = supLength * M_PI * (supDiameter + supThickness) * supThickness * suprho;			// supArea * suprho * supThickness
	mass[14] = supVolume * airDensitySUP;
	mass[15] = houseVolume * airDensityIN;
	if(roofIntRval > 0) {
	   mass[16] = area[16] * roofIntRval * 0.7;   // 0.7kg/m2/R-val = 20kg/m3 * 0.035W/mK for fiberglass
	   mass[17] = area[17] * roofIntRval * 0.7;
	   }

	// Specific heat capacities
	heatCap[0] = CpAir;
	heatCap[1] = 1210;													// CP plywood
	heatCap[3] = heatCap[1];
	heatCap[2] = cpShingles;
	heatCap[4] = heatCap[2];
	heatCap[5] = 1630;													// CP wood
	heatCap[6] = 1150;
	heatCap[7] = heatCap[6];
	heatCap[8] = heatCap[1];
	heatCap[9] = heatCap[8];
	heatCap[10] = retCp;												// input
	heatCap[11] = CpAir;
	heatCap[12] = 1300;												// combination of wood and drywall
	heatCap[13] = supCp;
	heatCap[14] = CpAir;
	heatCap[15] = CpAir;
	heatCap[16] = 1030;  // fiberglass: http://www.greenspec.co.uk/building-design/insulation-materials-thermal-properties/
	heatCap[17] = heatCap[16];

	// Thermal conductivities (k) [W/mK]and R-values [m2K/W] and the like
	kWood = 0.15;												// check with Iain about this
	//kAir = 0.02624;											// Thermal conductivity of air, now as function of air temperature
	kAir = 1.5207e-11 * pow(tempOld[15],3) - 4.8574e-8 * pow(tempOld[15],2) + 1.0184e-4 * tempOld[15] - 0.00039333;
	muAir = 0.000018462;										// Dynamic viscosity of air (mu) [kg/ms] Make temperature dependent  (this value at 300K)

	// U-values
	if(roofExtRval > 0)
	   roofInsulRatio = (1/roofExtRval)/(1/roofExtRval + 1/((woodThickness / kWood) + Rshingles));
	else
	   roofInsulRatio = 1;
	uVal[1] = 1 / ((woodThickness / kWood) + Rshingles + roofExtRval);
	uVal[2] = uVal[1];
	uVal[3] = uVal[1];
	uVal[4] = uVal[1];
	if(tempOld[15] > tempOld[0])
		uVal[6] = 1 / ceilRval + 0.085;   // From 2013 Res ACM table 2-2 (0.015 * 5.6783)
	else
		uVal[6] = 1 / ceilRval;
	uVal[7] = uVal[6];
	uVal[8] = 1 / gableEndRval;
	uVal[9] = uVal[8];
	// Inner Surface of Ducts
	// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
	// Note Use of HI notation
	// I think that the following may be an imperical relationship
	// Return Ducts
	HI = .023 * kAir / retDiameter * pow((retDiameter * airDensityRET * abs(retVel) / muAir), .8) * pow((CpAir * muAir / kAir), .4);
	uVal[10] = 1 / (retRval + 1/HI);
	// Supply Ducts
	HI = .023 * kAir / supDiameter * pow((supDiameter * airDensitySUP * supVel / muAir), .8) * pow((CpAir * muAir / kAir), .4);
	uVal[13] = 1 / (supRval + 1/HI);
	if(roofIntRval > 0) {
      uVal[16] = 1 / roofIntRval;
      uVal[17] = 1 / roofIntRval;
      }

	/* most of the surfaces in the attic undergo both natural and forced convection
	the overall convection is determined by the forced and natural convection coefficients
	to the THIRD power, adding them, and taking the cubed root.  This is Iain's idea
	and it seemed to work for him*/

	// Characteristic velocity
	characteristicVelocity = (matticenvin - matticenvout) / airDensityATTIC / AL4 / 2.0;
	if(characteristicVelocity == 0) {
		characteristicVelocity = abs(mCeiling) / airDensityATTIC / AL4 / 2.0;
		if(characteristicVelocity == 0)
			characteristicVelocity = .1;
	}

	// convection heat transfer coefficients
   htCoef[roofInNorth] = heatTranCoef(tempOld[roofInNorth], tempOld[0], characteristicVelocity);  // inner north sheathing
   htCoef[roofInSouth] = heatTranCoef(tempOld[roofInSouth], tempOld[0], characteristicVelocity);  // inner south sheathing
	htCoef[2] = heatTranCoef(tempOld[2], tempOut, windSpeed);                  // outer north sheathing
	htCoef[4] = heatTranCoef(tempOld[4], tempOut, windSpeed);                  // outer south sheathing
	htCoef[5] = heatTranCoef(tempOld[5], tempOld[0], characteristicVelocity);  // Wood (joists,truss,etc.)
   
	// Underside of Ceiling. Modified to use fixed numbers from ASHRAE Fundamentals ch.3 on 05/18/2000
	if(AHflag != 0)
		htCoef[6] = 9;
	else
		htCoef[6] = 6;
	// House Mass uses ceiling heat transfer coefficient as rest for house heat transfer coefficient	
	htCoef[12] = htCoef[6];

	htCoef[7] = heatTranCoef(tempOld[7], tempOld[0], characteristicVelocity);  // Attic Floor
	htCoef[8] = heatTranCoef(tempOld[8], tempOld[0], characteristicVelocity);  // Inner side of gable endwalls (lumped together)
	htCoef[9] = heatTranCoef(tempOld[9], tempOut, windSpeed);                 // Outer side of gable ends

	if(ductLocation == 1) { //  Ducts in the house
		// Outer Surface of Ducts
		if(AHflag != 0) {
			htCoef[10] = 9;
			htCoef[13] = htCoef[10];
		} else {
			htCoef[10] = 6;
			htCoef[13] = htCoef[10];
		}
	} else {
		htCoef[10] = heatTranCoef(tempOld[10], tempOld[0], characteristicVelocity);			// Outer Surface of Return Ducts
		htCoef[13] = heatTranCoef(tempOld[13], tempOld[0], characteristicVelocity);			// Outer Surface of Supply Ducts
	}

   // Pass back wood surface heat transfer coefficients for use by moisture routines
   innerNorthH = htCoef[1];
   innerSouthH = htCoef[3];
   bulkH = htCoef[5];


	// Radiation view factors
	/* 
	Only 5 nodes are involved in radiation transfer in the attic:
	North roof(1), South roof(3), Ceiling(7), Return duct(10), and Supply duct(13)
	The endwalls have a very small contribution to radiation exchange and are neglected.
	The wood may or may not contribute to radiation exchange, but their geometry is
	too complex to make any assumptions so it is excluded.
   */
   
   viewFactor[7][1] = 1 / 2.0;
   viewFactor[7][3] = viewFactor[7][1];
   viewFactor[1][7] = viewFactor[7][1] * area[7] / area[1];
   viewFactor[3][7] = viewFactor[1][7];
	if(ductLocation == 1) { //  Ducts in the house
		viewFactor[1][3] = (1 - viewFactor[1][7]);
		viewFactor[3][1] = (1 - viewFactor[3][7]);
	} else {			// ducts in the attic
		// 33.3% of each duct sees each sheathing surface (top third of duct)
		viewFactor[13][1] = 1 / 2.0;
		viewFactor[10][1] = viewFactor[13][1];
		viewFactor[13][3] = viewFactor[13][1];
		viewFactor[10][3] = viewFactor[13][1];

		viewFactor[1][13] = viewFactor[13][1] * (area[13] / 3) / area[1];
		viewFactor[1][10] = viewFactor[10][1] * (area[10] / 3) / area[1];
		viewFactor[3][13] = viewFactor[13][3] * (area[13] / 3) / area[3];
		viewFactor[3][10] = viewFactor[10][3] * (area[10] / 3) / area[3];
		viewFactor[1][3] = (1 - viewFactor[1][7] - viewFactor[1][10] - viewFactor[1][13]);
		viewFactor[3][1] = (1 - viewFactor[3][7] - viewFactor[3][10] - viewFactor[3][13]);

		// North Sheathing
		rtCoef[roofInNorth][10] = radTranCoef(woodEmissivity, tempOld[roofInNorth], tempOld[10], viewFactor[1][10], area[1]/(area[10]/3));
		rtCoef[roofInNorth][13] = radTranCoef(woodEmissivity, tempOld[roofInNorth], tempOld[13], viewFactor[1][13], area[1]/(area[13]/3));

		// South Sheathing
		rtCoef[roofInSouth][10] = radTranCoef(woodEmissivity, tempOld[roofInSouth], tempOld[10], viewFactor[3][10], area[3]/(area[10]/3));
		rtCoef[roofInSouth][13] = radTranCoef(woodEmissivity, tempOld[roofInSouth], tempOld[13], viewFactor[3][13], area[3]/(area[13]/3));

		// Return Ducts (note, No radiative exchange w/ supply ducts)
		rtCoef[10][roofInSouth] = radTranCoef(woodEmissivity, tempOld[10], tempOld[roofInSouth], viewFactor[10][3], area[10]/area[3]);
		rtCoef[10][roofInNorth] = radTranCoef(woodEmissivity, tempOld[10], tempOld[roofInNorth], viewFactor[10][1], area[10]/area[1]);

		// Supply Ducts (note, No radiative exchange w/ return ducts)
		rtCoef[13][roofInSouth] = radTranCoef(woodEmissivity, tempOld[13], tempOld[roofInSouth], viewFactor[13][3], area[13]/area[3]);
		rtCoef[13][roofInNorth] = radTranCoef(woodEmissivity, tempOld[13], tempOld[roofInNorth], viewFactor[13][1], area[13]/area[1]);
	}

	// North Sheathing
	rtCoef[roofInNorth][roofInSouth] = radTranCoef(woodEmissivity, tempOld[roofInNorth], tempOld[roofInSouth], viewFactor[1][3], area[1]/area[3]);
	rtCoef[roofInNorth][7] = radTranCoef(woodEmissivity, tempOld[roofInNorth], tempOld[7], viewFactor[1][7], area[1]/area[7]);

	// South Sheathing
	rtCoef[roofInSouth][roofInNorth] = radTranCoef(woodEmissivity, tempOld[roofInSouth], tempOld[roofInNorth], viewFactor[3][1], area[3]/area[1]);
	rtCoef[roofInSouth][7] = radTranCoef(woodEmissivity, tempOld[roofInSouth], tempOld[7], viewFactor[3][7], area[3]/area[7]);

	// Attic Floor
	rtCoef[7][roofInSouth] = radTranCoef(woodEmissivity, tempOld[7], tempOld[roofInSouth], viewFactor[7][3], area[7]/area[3]);
	rtCoef[7][roofInNorth] = radTranCoef(woodEmissivity, tempOld[7], tempOld[roofInNorth], viewFactor[7][1], area[7]/area[1]);

	// underside of ceiling
	rtCoef[6][12] = radTranCoef(woodEmissivity, tempOld[6], tempOld[12], 1, area[6]/area[12]);

	// Sky and ground radiation
	FRS = (1 - skyCover) * (180 - roofPitch) / 180;      	// ROOF-SKY SHAPE FACTOR
	FG = 1 - FRS;                            					// ROOF-GROUND SHAPE FACTOR
	TGROUND = tempOut;                           			// ASSUMING GROUND AT AIR TEMP
	if(skyCover < 1) {
		skyCoef4 = radTranCoef(roofEmissivity, tempOld[4], TSKY, FRS, 0);
		skyCoef2 = radTranCoef(roofEmissivity, tempOld[2], TSKY, FRS, 0);
	} else {
		skyCoef4 = 0;
		skyCoef2 = 0;
	}
	gndCoef4 = radTranCoef(roofEmissivity, tempOld[4], TGROUND, FG, 0);
	gndCoef2 = radTranCoef(roofEmissivity, tempOld[2], TGROUND, FG, 0);


	// ITERATION OF TEMPERATURES WITHIN HEAT SUBROUTINE
	// THIS ITERATES BETWEEN ALL TEMPERATURES BEFORE RETURNING TO MAIN PROGRAM
	heatIterations = 0;
	for(int i=0; i < attic_nodes; i++) {
		toldcur[i] = tempOld[i];
	}
	while(1) {

		heatIterations++;
		
		// reset array A to 0
		if(heatIterations > 1) {
			A.clear();
			A.resize(attic_nodes, vector<double>(attic_nodes+1));
		}

		// NODE 0 IS ATTIC AIR
		if(mCeiling >= 0) {
			// flow from attic to house
			A[0][0] = mass[0] * heatCap[0] / dtau + htCoef[7] * area[7] + htCoef[5] * area[5] + mCeiling * heatCap[0]
			        + mSupAHoff * heatCap[14] + mRetAHoff * heatCap[11] + htCoef[3] * area[3] + htCoef[1] * area[1]
			        + area[8] * htCoef[8] - matticenvout * heatCap[0] - mRetLeak * heatCap[0];
			b[0] = mass[0] * heatCap[0] * tempOld[0] / dtau + matticenvin * heatCap[0] * tempOut + mSupLeak * heatCap[0] * toldcur[14];
		} else {
			// flow from house to attic
			A[0][0] = mass[0] * heatCap[0] / dtau + htCoef[7] * area[7] + htCoef[5] * area[5] + htCoef[3] * area[3]
			        + htCoef[1] * area[1] + area[8] * htCoef[8] - matticenvout * heatCap[0] - mRetLeak * heatCap[0];
			b[0] = mass[0] * heatCap[0] * tempOld[0] / dtau - mCeiling * heatCap[0] * toldcur[15] - mSupAHoff * heatCap[14] * toldcur[14]
			     - mRetAHoff * heatCap[11] * toldcur[11] + matticenvin * heatCap[0] * tempOut + mSupLeak * heatCap[14] * toldcur[14];
		}
		A[0][1] = -htCoef[1] * area[1];
		A[0][3] = -htCoef[3] * area[3];
		A[0][5] = -htCoef[5] * area[5];
		A[0][7] = -htCoef[7] * area[7];
		A[0][8] = -htCoef[8] * area[8];
		if(ductLocation == 0) {		// duct surface conduction loss to attic
			A[0][0] +=  htCoef[13] * area[13] / 2 + htCoef[10] * area[10] / 2;
			A[0][10] = -htCoef[10] * area[10] / 2;
			A[0][13] = -htCoef[13] * area[13] / 2;
		}

		// NODE 1 IS INSIDE NORTH SHEATHING
   	if(roofIntRval > 0) {
         A[1][1] = mass[1] * heatCap[1] / dtau + uVal[16] * area[1] + area[1] * uVal[1];
         b[1] = mass[1] * heatCap[1] * tempOld[1] / dtau;
         A[1][2] = -area[1] * uVal[1];
         A[1][16] = -area[1] * uVal[16];
      }
      else {
         A[1][0] = -htCoef[1] * area[1];
         A[1][1] = mass[1] * heatCap[1] / dtau + htCoef[1] * area[1] + area[1] * uVal[1] + rtCoef[1][3] * area[1] + rtCoef[1][7] * area[1];
         b[1] = mass[1] * heatCap[1] * tempOld[1] / dtau;
         A[1][2] = -area[1] * uVal[1];
         A[1][3] = -rtCoef[1][3] * area[1];
         A[1][7] = -rtCoef[1][7] * area[1];

         if(ductLocation == 0) {			// duct surface radiation to sheathing
            A[1][1] += rtCoef[1][10] * area[1] + rtCoef[1][13] * area[1];
            A[1][10] = -rtCoef[1][10] * area[1];
            A[1][13] = -rtCoef[1][13] * area[1];
         }
      }

		// NODE 2 IS OUTSIDE NORTH SHEATHING
		A[2][1] = -area[1] * uVal[2];
		A[2][2] = mass[2] * heatCap[2] / dtau + htCoef[2] * area[2] + area[1] * uVal[2] + skyCoef2 * area[1] + gndCoef2 * area[1];
		b[2] = mass[2] * heatCap[2] * tempOld[2] / dtau + htCoef[2] * area[2] * tempOut + area[1] * nsolrad * roofAbsorptivity
		     + skyCoef2 * area[1] * TSKY + gndCoef2 * area[1] * TGROUND;

		// NODE 3 IS INSIDE SOUTH SHEATHING
   	if(roofIntRval > 0) {
         A[3][3] = mass[3] * heatCap[3] / dtau + uVal[17] * area[3] + area[3] * uVal[3];
         b[3] = mass[3] * heatCap[3] * tempOld[3] / dtau;
         A[3][4] = -area[3] * uVal[3];
         A[3][17] = -area[3] * uVal[17];
      }
      else {
         A[3][0] = -htCoef[3] * area[3];
         A[3][1] = -rtCoef[3][1] * area[3];
         A[3][3] = mass[3] * heatCap[3] / dtau + htCoef[3] * area[3] + area[3] * uVal[3] + rtCoef[3][1] * area[3] + rtCoef[3][7] * area[3];
         b[3] = mass[3] * heatCap[3] * tempOld[3] / dtau;
         A[3][4] = -area[3] * uVal[3];
         A[3][7] = -rtCoef[3][7] * area[3];

         if(ductLocation == 0) {			// duct surface radiation to sheathing
            A[3][3] += rtCoef[3][10] * area[3] + rtCoef[3][13] * area[3];
            A[3][10] = -rtCoef[3][10] * area[3];
            A[3][13] = -rtCoef[3][13] * area[3];
         }
      }

		// NODE 4 IS OUTSIDE SOUTH SHEATHING
		A[4][3] = -area[3] * uVal[4];
		A[4][4] = mass[4] * heatCap[4] / dtau + htCoef[4] * area[4] + area[3] * uVal[4] + skyCoef4 * area[3] + gndCoef4 * area[3];
		b[4] = mass[4] * heatCap[4] * tempOld[4] / dtau + htCoef[4] * area[4] * tempOut + area[3] * ssolrad * roofAbsorptivity
		     + skyCoef4 * area[3] * TSKY + gndCoef4 * area[3] * TGROUND;

		// NODE 5 IS MASS OF WOOD IN ATTIC I.E. JOISTS AND TRUSSES
		A[5][0] = -htCoef[5] * area[5];
		A[5][5] = mass[5] * heatCap[5] / dtau + htCoef[5] * area[5];
		b[5] = mass[5] * heatCap[5] * tempOld[5] / dtau;

		// NODE 6 ON INSIDE OF CEILING
		A[6][6] = mass[6] * heatCap[6] / dtau + htCoef[6] * area[6] + rtCoef[6][12] * area[6] + area[6] * uVal[6];
		b[6] = mass[6] * heatCap[6] / dtau * tempOld[6];
		A[6][7] = -area[6] * uVal[6];
		A[6][15] = -htCoef[6] * area[6];
		A[6][12] = -rtCoef[6][12] * area[6];

		// NODE 7 ON ATTIC FLOOR
		A[7][0] = -htCoef[7] * area[7];
		A[7][roofInNorth] = -rtCoef[7][roofInNorth] * area[7];
		A[7][roofInSouth] = -rtCoef[7][roofInSouth] * area[7];
		A[7][6] = -area[7] * uVal[7];
		A[7][7] = mass[7] * heatCap[7] / dtau + htCoef[7] * area[7] + rtCoef[7][roofInNorth] * area[7] + rtCoef[7][roofInSouth] * area[7] + area[7] * uVal[7];				// + HR8t11 * area[7] + HR8t14 * area[7]
		b[7] = mass[7] * heatCap[7] / dtau * tempOld[7];

		// NODE 8 IS INSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[8][0] = -htCoef[8] * area[8];
		A[8][8] = mass[8] * heatCap[8] / dtau + htCoef[8] * area[8] + area[8] * uVal[8];
		A[8][9] = -area[8] * uVal[8];
		b[8] = mass[8] * heatCap[8] * tempOld[8] / dtau;

		// NODE 9 IS OUTSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[9][8] = -area[9] * uVal[9];
		A[9][9] = mass[9] * heatCap[9] / dtau + htCoef[9] * area[9] + area[9] * uVal[9];
		b[9] = mass[9] * heatCap[9] * tempOld[9] / dtau + htCoef[9] * area[9] * tempOut;

		// NODE 10 Exterior Return Duct Surface
		// Remember that the fluid properties are evaluated at a constant temperature
		// therefore, the convection on the inside of the ducts is
		A[10][11] = -area[11] * uVal[10];
		b[10] = mass[10] * heatCap[10] * tempOld[10] / dtau;
		if(ductLocation == 1) {			// ducts in house
			A[10][10] = mass[10] * heatCap[10] / dtau + htCoef[10] * area[10] + area[11] * uVal[10];
			A[10][15] = -area[10] * htCoef[10];
		} else {
			A[10][0] = -area[10] * htCoef[10] / 2;
			A[10][roofInNorth] = -area[10] * rtCoef[10][roofInNorth] / 3;
			A[10][roofInSouth] = -area[10] * rtCoef[10][roofInSouth] / 3;
			A[10][10] += mass[10] * heatCap[10] / dtau + htCoef[10] * area[10] / 2 + area[11] * uVal[10]
			          + area[10] * rtCoef[10][roofInNorth] / 3 + area[10] * rtCoef[10][roofInSouth] / 3;
		}
		
		// NODE 11 Air in return duct
		A[11][10] = -area[11] * uVal[10];
		if(mCeiling >= 0) {
			// flow from attic to house
			A[11][11] = mass[11] * heatCap[11] / dtau + area[11] * uVal[10] + mAH * heatCap[11] + mRetAHoff * heatCap[11];
			b[11] = mass[11] * heatCap[11] * tempOld[11] / dtau + mRetAHoff * heatCap[0] * toldcur[0]
			      - mRetLeak * heatCap[0] * toldcur[0] - mRetReg * heatCap[0] * toldcur[15]
			      - mFanCycler * heatCap[0] * tempOut - mHRV_AH * heatCap[15] * ((1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15])
			      - mERV_AH * heatCap[15] * ((1-ERV_SRE) * tempOut + ERV_SRE * tempOld[15]);
			
		} else {
			// flow from house to attic
			A[11][11] = mass[11] * heatCap[11] / dtau + area[11] * uVal[10] + mAH * heatCap[11] - mRetAHoff * heatCap[11];
			b[11] = mass[11] * heatCap[11] * tempOld[11] / dtau - mRetAHoff * heatCap[15] * toldcur[15]
			      - mRetLeak * heatCap[0] * toldcur[0] - mRetReg * heatCap[0] * toldcur[15]
			      - mFanCycler * heatCap[0] * tempOut - mHRV_AH * heatCap[15] * ((1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15])
			      - mERV_AH * heatCap[15] * ((1-ERV_SRE) * tempOut + ERV_SRE * tempOld[15]);
		}

		// Node 12 is the mass of the structure of the house that interacts
		// with the house air to increase its effective thermal mass
		// 95% of solar gain goes to house mass, 5% to house air
		A[12][12] = mass[12] * heatCap[12] / dtau + htCoef[12] * area[12] + rtCoef[6][12] * area[6];
		A[12][15] = -htCoef[12] * area[12];
		A[12][6] = -rtCoef[6][12] * area[6];
		b[12] = mass[12] * heatCap[12] * tempOld[12] / dtau + .95 * solgain;

		// NODE 13 Exterior Supply Duct Surface
		b[13] = mass[13] * heatCap[13] * tempOld[13] / dtau;
		A[13][14] = -area[14] * uVal[13];
		if(ductLocation == 1) {			// ducts in house
			A[13][13] = mass[13] * heatCap[13] / dtau + htCoef[13] * area[13] + area[14] * uVal[13];
			A[13][15] = -area[13] * htCoef[13];
		} else {
			A[13][0] = -area[13] * htCoef[13] / 2;
			A[13][roofInNorth] = -area[13] * rtCoef[13][roofInNorth] / 3;
			A[13][roofInSouth] = -area[13] * rtCoef[13][roofInSouth] / 3;
			A[13][13] = mass[13] * heatCap[13] / dtau + htCoef[13] * area[13] / 2 + area[14] * uVal[13]
			          + area[13] * rtCoef[13][roofInNorth] / 3 + area[13] * rtCoef[13][roofInSouth] / 3;
		}

		// NODE 14 Air in SUPPLY duct
		// capacity is AC unit capacity in Watts
		// this is a sensible heat balance, the moisture is balanced in a separate routine.
		A[14][13] = -area[14] * uVal[13];
		if(mCeiling >= 0) {
			// flow from attic to house
			A[14][14] = mass[14] * heatCap[14] / dtau + area[14] * uVal[13] + mSupReg * heatCap[14]
			          + mSupLeak * heatCap[14] + mSupAHoff * heatCap[14];
			b[14] = mass[14] * heatCap[14] * tempOld[14] / dtau - capacityc + capacityh + evapcap
			      + mAH * heatCap[11] * toldcur[11] + mSupAHoff * heatCap[0] * toldcur[0];
		} else {
			// flow from house to attic
			A[14][14] = mass[14] * heatCap[14] / dtau + area[14] * uVal[13] + mSupReg * heatCap[14]
			          + mSupLeak * heatCap[14] - mSupAHoff * heatCap[14];
			b[14] = mass[14] * heatCap[14] * tempOld[14] / dtau - capacityc + capacityh + evapcap
			      + mAH * heatCap[11] * toldcur[11] - mSupAHoff * heatCap[15] * toldcur[15];
		}

		// NODE 15 AIR IN HOUSE
		// use solair tmeperature for house UA
		if(mCeiling >= 0) {
			// flow from attic to house
			A[15][15] = mass[15] * heatCap[15] / dtau + htCoef[6] * area[6] - mRetReg * heatCap[15]
			          - mHouseOUT * heatCap[15] + htCoef[12] * area[12] + uaSolAir + uaTOut;
			b[15] = mass[15] * heatCap[15] * tempOld[15] / dtau + (mHouseIN - mHRV) * heatCap[15] * tempOut
			      + mHRV * heatCap[15] * (( 1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15])
			      + uaSolAir * tsolair + uaTOut * tempOut + .05 * solgain + mSupReg * heatCap[0] * toldcur[14] 
				   + mCeiling * heatCap[0] * toldcur[0] + mSupAHoff * heatCap[14] * toldcur[14]
				   + mRetAHoff * heatCap[11] * toldcur[11] + internalGains + dhSensibleGain;
		} else {
			// flow from house to attic
			A[15][15] = mass[15] * heatCap[15] / dtau + htCoef[6] * area[6] - mCeiling * heatCap[15]
			          - mSupAHoff * heatCap[15] - mRetAHoff * heatCap[15] - mRetReg * heatCap[15]
			          - mHouseOUT * heatCap[15] + htCoef[12] * area[12] + uaSolAir + uaTOut;
			b[15] = mass[15] * heatCap[15] * tempOld[15] / dtau + (mHouseIN - mHRV) * heatCap[15] * tempOut
			      + mHRV * heatCap[15] * ((1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15]) + uaSolAir * tsolair
			      + uaTOut * tempOut + .05 * solgain + mSupReg * heatCap[0] * toldcur[14] + internalGains + dhSensibleGain;
		}
		A[15][6] = -htCoef[6] * area[6];
		A[15][12] = -htCoef[12] * area[12];
		if(ductLocation == 1) {
			// ducts in house
			A[15][15] += area[10] * htCoef[10] + area[13] * htCoef[13];
			A[15][10] = -area[10] * htCoef[10];
			A[15][13] = -area[13] * htCoef[13];
		}

   	if(roofIntRval > 0) {
         // NODE 16 IS INSIDE NORTH Insulation
         A[16][0] = -htCoef[16] * area[16];
         A[16][1] = -area[16] * uVal[16];
         A[16][16] = mass[16] * heatCap[16] / dtau + htCoef[16] * area[16] + area[16] * uVal[16] + rtCoef[16][17] * area[16] + rtCoef[16][7] * area[16];
         b[16] = mass[16] * heatCap[16] * tempOld[16] / dtau;
         A[16][17] = -rtCoef[16][17] * area[16];
         A[16][7] = -rtCoef[16][7] * area[16];

         if(ductLocation == 0) {			// duct surface radiation to insulation
            A[16][16] += rtCoef[16][10] * area[16] + rtCoef[16][13] * area[16];
            A[16][10] = -rtCoef[16][10] * area[16];
            A[16][13] = -rtCoef[16][13] * area[16];
         }

         // NODE 17 IS INSIDE SOUTH Insulation
         A[17][0] = -htCoef[17] * area[17];
         A[17][3] = -area[17] * uVal[17];
         A[17][17] = mass[17] * heatCap[17] / dtau + htCoef[17] * area[17] + area[17] * uVal[17] + rtCoef[17][16] * area[17] + rtCoef[17][7] * area[17];
         b[17] = mass[17] * heatCap[17] * tempOld[17] / dtau;
         A[17][16] = -rtCoef[17][16] * area[17];
         A[17][7] = -rtCoef[17][7] * area[17];

         if(ductLocation == 0) {			// duct surface radiation to insulation
            A[17][17] += rtCoef[17][10] * area[17] + rtCoef[17][13] * area[17];
            A[17][10] = -rtCoef[17][10] * area[17];
            A[17][13] = -rtCoef[17][13] * area[17];
         }
      }

		for (int i=0; i<attic_nodes; i++) {
			A[i][attic_nodes] = b[i];
		}
		b = gauss(A);

		if(abs(b[0] - toldcur[0]) < .1) {
			break;
		} else {
			for(int i=0; i < attic_nodes; i++) {
				toldcur[i] = b[i];
			}
		}
	} // END of DO LOOP
	for (int i=0; i<attic_nodes; i++) {
		x[i] = b[i];
		}
}


void sub_houseLeak (
	int& AHflag,
	int& leakIterations, 
	double& windSpeed, 
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
	) {
		double Cpwallvar = 0;		// This variable replaces the non-array wallCp var
		double CPvar = 0;			// This variable replaces the non-array CP var
		double CP[4][4];
		double dPint;
		double mWallIn, mWallOut;
		//double rhoi;
		//double rhoo;
		//double rhoa;
		//double fluePressureExp;
		double dPwind;
		double dPtemp;
		//double Cproof;
		double Cpwalls;
		//double Cpattic;		
		double Cpfloor;
		double Cwall;
		double Cfloor;
		
		mFlue = 0;
		mCeiling = 0;

		for(int i=0; i < 4; i++) {
			mFloor[i] = 0;			
			wallCp[i] = 0;
			for(int j=0; j < 4; j++) {
				CP[i][j] = 0;
			}
		}

		// the following are pressure differences common to all the flow equations
		dPwind = airDensityOUT / 2 * pow(windSpeed,2);
		dPtemp = airDensityOUT * g * (tempHouse - tempOut) / tempHouse;
		
		// the following are some typical pressure coefficients for rectangular houses
		
		//Cproof = -.4;

		for(int i=0; i < 4; i++) {
			CP[i][0] = .6;
			CP[i][1] = -.3;
		}
		// here the variation of each wall Cp with wind angle is accounted for:
		// for row houses:

		if(rowHouse) {
			CP[0][2] = -0.2;
			CP[0][3] = -0.2;
			CP[1][2] = -0.2;
			CP[1][3] = -0.2;
			CP[2][2] = -0.65;
			CP[2][3] = -0.65;
			CP[3][2] = -0.65;
			CP[3][3] = -0.65;
		} else {
			// for isolated houses
			for(int i=0; i < 4; i++) {
				CP[i][2] = -0.65;
				CP[i][3] = -0.65;
			}
		}

		f_CpTheta(CP, windAngle, wallCp);

		Cpwalls = 0;
		//Cpattic = 0;

		for(int i=0; i < 4; i++) {
			//Cpattic = Cpattic + Sw[i] * wallCp[i] * soffitFraction[i];
			Cpwalls = Cpwalls + Sw[i] * wallCp[i] * wallFraction[i];			// Shielding weighted Cp
		}

		//Cpattic = Cpattic + pow(flueShelterFactor, 2) * Cproof * soffitFraction[4];

		//if(leakIterations < 1) {        // Yihuan: delete the if condition for the flag
			Pint = 0;			// a reasonable first guess   
			dPint = 200;		// increased from 25 to account for economizer operation
			//dPint = 25;		// increased from 25 to account for economizer operation
		//} else {
			//dPint = .25;
	//}

		do {
			mIN = 0;
			mOUT = 0;
			
			if(numFlues) {					//FF: This IF behaves as if(numFlues != 0)
				f_flueFlow(tempHouse, flueShelterFactor, dPwind, dPtemp, eaveHeight, Pint, numFlues, flue, mFlue, airDensityOUT, airDensityIN, dPflue, tempOut, houseVolume, windPressureExp);

				if(mFlue >= 0) {
					mIN = mIN + mFlue;		// Add mass flow through flue
				} else {
					mOUT = mOUT + mFlue;
				}
			}
			
			

			if(leakFracFloor > 0) {
				if(Crawl == 1) {
					// for a crawlspace the flow is put into array position 1
					Cpfloor = Cpwalls;
					Cfloor = envC * leakFracFloor;

					f_floorFlow3(Cfloor, Cpfloor, dPwind, Pint, n, mFloor[0], airDensityOUT, airDensityIN, Hfloor, dPtemp);
					
					if(mFloor[0] >= 0) {
						mIN = mIN + mFloor[0];
					} else {
						mOUT = mOUT + mFloor[0];
					}

				} else {
					for(int i=0; i < 4; i++) {
						Cpfloor = Sw[i] * wallCp[i];
						Cfloor = envC * leakFracFloor * floorFraction[i];

						f_floorFlow3(Cfloor, Cpfloor, dPwind, Pint, n, mFloor[i], airDensityOUT, airDensityIN, Hfloor, dPtemp);
						
						if(mFloor[i] >= 0) {							
							mIN = mIN + mFloor[i];
						} else {
							mOUT = mOUT + mFloor[i];
						}
					}					
				}
			}
			
			if(leakFracCeil > 0) {

				double CCeiling = envC * leakFracCeil;
				f_ceilingFlow(AHflag, Patticint, eaveHeight, dPtemp, dPwind, Pint, envC, n, mCeiling, atticC, airDensityATTIC, airDensityIN, tempAttic, tempHouse, tempOut, airDensityOUT, mSupAHoff, mRetAHoff, supC, supn, retC, retn, CCeiling);

				if(mCeiling >= 0) {
					mIN = mIN + mCeiling + mSupAHoff + mRetAHoff;
				} else {
					mOUT = mOUT + mCeiling + mSupAHoff + mRetAHoff;
				}
			}

			
			if(leakFracWall > 0) {
				for(int i=0; i < 4; i++) {
					Cpwallvar = Sw[i] * wallCp[i];
					Cwall = envC * leakFracWall * wallFraction[i];
					
					f_wallFlow3(tempHouse, tempOut, airDensityIN, airDensityOUT, Cpwallvar, n, Cwall, eaveHeight, Pint, dPtemp, dPwind, mWallIn, mWallOut, Hfloor);
					
					mIN = mIN + mWallIn;
					mOUT = mOUT + mWallOut;
				}
			}

			for(int i=0; i < numFans; i++) {
				if(fan[i].on == 1) {						// for cycling fans they will somtimes be off and we don;t want ot include them

					f_fanFlow(fan[i], airDensityOUT, airDensityIN);					
					
					if(fan[i].m >= 0) {
						mIN = mIN + fan[i].m;
					} else {
						mOUT = mOUT + fan[i].m;
					}
				}
			}
			
			for(int i=0; i < numPipes; i++) {

				// FF: This if is to counter the 0 index out of bound present in BASIC version
				// example: if pipe[i].wall = 0 then it would look for Sw[0] in BASIC version, which defaults to 0 while in the C++ version
				// it would try to find Sw[-1] and cause error.

				if(Pipe[i].wall - 1 >= 0)
					CPvar = Sw[Pipe[i].wall - 1] * wallCp[Pipe[i].wall - 1];
				else
					CPvar = 0;

				f_pipeFlow(airDensityOUT, airDensityIN, CPvar, dPwind, dPtemp, Pint, Pipe[i], tempHouse, tempOut);
				
				if(Pipe[i].m >= 0) {
					mIN = mIN + Pipe[i].m;
				} else {
					mOUT = mOUT + Pipe[i].m;
				}
			}
			
			for(int i=0; i < numWinDoor; i++) {
				if(winDoor[i].wall-1 >= 0) {
					Cpwallvar = Sw[winDoor[i].wall-1] * wallCp[winDoor[i].wall-1];
					f_winDoorFlow(tempHouse, tempOut, airDensityIN, airDensityOUT, eaveHeight, Cpwallvar, n, Pint, dPtemp, dPwind, winDoor[i]);
					mIN = mIN + winDoor[i].mIN;
					mOUT = mOUT + winDoor[i].mOUT;
					}
			}

			// DUCT MASS FLOWS
			// Msup is flow out of supply registers plus leakage to inside
			// Mret is flow into return registers plus leakage to inside

			mIN = mIN + mSupReg;
			mOUT = mOUT + mRetReg; // Note Mret should be negative

			Pint = Pint - copysign(dPint, mIN + mOUT);
			dPint = dPint / 2;
		} while (dPint > .0001);
		//} while (dPint > .01);

		if(mCeiling >= 0) { // flow from attic to house
			mHouseIN = mIN - mCeiling - mSupReg - mSupAHoff - mRetAHoff;
			mHouseOUT = mOUT - mRetReg;
		} else {
			mHouseIN = mIN - mSupReg;
			mHouseOUT = mOUT - mCeiling - mRetReg - mSupAHoff - mRetAHoff;
		}
		// nopressure:
}

void sub_atticLeak ( 
	int& leakIterations, 
	double& windSpeed, 
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
) {

	double mRoofIn, mRoofOut;
	double wallCp[4];
	double Batto[4];
	double Cproof[4];
	double Cppitch[4];
	double CP[4][4];

	//double rhoi;
	//double rhoo;
	//double rhoa;
	double mAtticFloor;
	double dPwind;
	double dPtemp;
	double dPatticint;
	double Croof;
	double Cpr;
	double CPvar;

	for(int i=0; i < 4; i++) {
		wallCp[i] = 0;
		Batto[i] = 0;
		Cproof[i] = 0;
		Cppitch[i] = 0;
		for(int j=0; j < 4; j++) {
			CP[i][j] = 0;
		}
	}

	mAtticFloor = -mCeiling;

	// the following are pressure differences common to all the flow equations
	dPwind = airDensityOUT / 2 * pow(windSpeed, 2);
	dPtemp = airDensityOUT * g * (tempAttic - tempOut) / tempAttic;

	/*if(dPwind == 0 && dPtemp == 0) {
		Qnet = 0;
		return;
	}*/

	// the following are some typical pressure coefficients for pitched roofs
	if(roofPitch < 10) {
		Cproof[0] = -.8;
		Cproof[1] = -.4;
	} else if(roofPitch > 30) {
		Cproof[0] = .3;
		Cproof[1] = -.5;
	} else {
		Cproof[0] = -.4;
		Cproof[1] = -.4;
	}
	if(rowHouse) {
		Cproof[2] = -.2;
		Cproof[3] = -.2;
	} else {
		// for isolated houses
		Cproof[2] = -.6;
		Cproof[3] = -.6;
	}
	if(roofPeakPerpendicular) {
		Cproof[2] = Cproof[0];
		Cproof[3] = Cproof[1];
		if(rowHouse) {
			Cproof[0] = -.2;
			Cproof[1] = -.2;
		} else {
			// for isolated houses
			Cproof[0] = -.6;
			Cproof[1] = -.6;
		}
	}

	// ***************************
	// the following are some typical pressure coefficients for rectangular houses

	for(int i=0; i < 4; i++) {
		CP[i][0] = .6;
		CP[i][1] = -.3;
	}

	// here the variation of each wall Cp with wind angle is accounted for:
	// for row houses:

	if(rowHouse) {
		CP[0][2] = -.2;
		CP[0][3] = -.2;
		CP[1][2] = -.2;
		CP[1][3] = -.2;
		CP[2][2] = -.65;
		CP[2][3] = -.65;
		CP[3][2] = -.65;
		CP[3][3] = -.65;
	} else {
		// for isolated houses
		for(int i=0; i < 4; i++) {
			CP[i][2] = -.65;
			CP[i][3] = -.65;
		}
	}

	f_CpTheta(CP, windAngle, wallCp);

	// for pitched roof leaks
	f_roofCpTheta(Cproof, windAngle, Cppitch, roofPitch);
	
	if(leakIterations < 2) {
		Patticint = 0;            // a reasonable first guess
		dPatticint = 25;
	} else {
		dPatticint = .25;
	}

	do {
		mAtticIN = 0;
		mAtticOUT = 0;
		
		// the following section is for the pitched part of the roof where the two pitched faces are assumed to have the same leakage
		Croof = atticC * soffitFraction[4] / 2;
		
		// for first pitched part either front, above wall 1, or side above wall 3
		if(roofPeakPerpendicular) {
			Cpr = Cppitch[2] * Sw[2];
		} else {
			Cpr = Cppitch[0] * Sw[0];
		}

		// developed from wallflow3:
		f_roofFlow(tempAttic, tempOut, airDensityATTIC, airDensityOUT, Cpr, atticPressureExp, Croof, roofPeakHeight, Patticint, dPtemp, dPwind, mRoofIn, mRoofOut, eaveHeight);

		mAtticIN = mAtticIN + mRoofIn;
		mAtticOUT = mAtticOUT + mRoofOut;

		// for second pitched part either back, above wall 2, or side above wall 4
		if(roofPeakPerpendicular) {
			Cpr = Cppitch[3] * Sw[3];
		} else {
			Cpr = Cppitch[1] * Sw[1];
		}

		// developed from wallflow3:
		f_roofFlow(tempAttic, tempOut, airDensityATTIC, airDensityOUT, Cpr, atticPressureExp, Croof, roofPeakHeight, Patticint, dPtemp, dPwind, mRoofIn, mRoofOut, eaveHeight);

		mAtticIN = mAtticIN + mRoofIn;
		mAtticOUT = mAtticOUT + mRoofOut;

		for(int i=0; i < numAtticVents; i++) {

			// FF: This if is to counter the 0 index out of bound present in BASIC version
			// example: if atticVent[i].wall = 0 then it would look for Sw[0] in BASIC version, which defaults to 0 while in the C++ version
			// it would try to find Sw[-1] and cause error.
			if(atticVent[i].wall - 1 >= 0)
				CPvar = Sw[atticVent[i].wall - 1] * Cppitch[atticVent[i].wall - 1];
			else
				CPvar = 0;

			f_atticVentFlow(airDensityOUT, airDensityATTIC, CPvar, dPwind, dPtemp, Patticint, atticVent[i], tempAttic, tempOut);
			
			if(atticVent[i].m >= 0) {
				mAtticIN = mAtticIN + atticVent[i].m;
			} else {
				mAtticOUT = mAtticOUT + atticVent[i].m;
			}
		}
		// note that gable vents are the same as soffits
		for(int i=0; i < 4; i++) {
			CPvar = Sw[i] * wallCp[i];

			f_soffitFlow(airDensityOUT, airDensityATTIC, CPvar, dPwind, dPtemp, Patticint, soffit[i], soffitFraction[i], atticC, atticPressureExp, tempAttic, tempOut);

			if(soffit[i].m >= 0)
				mAtticIN = mAtticIN + soffit[i].m;
			else
				mAtticOUT = mAtticOUT + soffit[i].m;
		}
		
		for(int i=0; i < numAtticFans; i++) {
			if(atticFan[i].on == 1) {        // for cycling fans hey will somtimes be off and we don;t want ot include them

				f_atticFanFlow(atticFan[i], airDensityOUT, airDensityATTIC);

				if(atticFan[i].m >= 0)
					mAtticIN = mAtticIN + atticFan[i].m;
				else
					mAtticOUT = mAtticOUT + atticFan[i].m;
			}
		}

		// Attic floor flow is determined first by HOUSELEAK
		// this fixed flow rate will fix the Patticint that is then
		// passed back to HOUSELEAK as the interior pressure of the attic

		// note that mattic floor has a sign change so that inflow to the attic is positive for a negative mCeiling
		if(mAtticFloor >= 0)
			mAtticIN = mAtticIN + mAtticFloor - mSupAHoff - mRetAHoff;
		else
			mAtticOUT = mAtticOUT + mAtticFloor - mSupAHoff - mRetAHoff;

		mAtticIN = mAtticIN + mSupLeak;
		mAtticOUT = mAtticOUT + mRetLeak;

		Patticint = Patticint - copysign(dPatticint, mAtticIN + mAtticOUT);
		dPatticint = dPatticint / 2;
	} while(dPatticint > .0001);

	if(mAtticFloor >= 0) {
		matticenvin = mAtticIN - mAtticFloor - mSupLeak + mSupAHoff + mRetAHoff;
		matticenvout = mAtticOUT - mRetLeak;
	} else {
		matticenvin = mAtticIN - mSupLeak;
		matticenvout = mAtticOUT - mAtticFloor - mRetLeak + mSupAHoff + mRetAHoff;
	}
}

void sub_filterLoading ( 
	int& MERV,
	int& loadingRate,
	int& AHMotorType,
	double& A_qAH_heat,
	double& A_qAH_cool,
	double& A_wAH_heat, 
	double& A_wAH_cool, 
	double& A_DL, 
	double& k_qAH,
	double& k_wAH,
	double& k_DL,
	double& qAH_heat0, 
	double& qAH_cool0,
	double& qAH_low
	) {
		int MERV_n = 0;		// Index for filter loading coefficient array (do not change)
		
		switch (MERV) {		// Choice of MERV filter based on user input from main function
		case 0:		// Equivalent to MERV 5 (used as error trapping)
			MERV_n = 0;
			break;

		case 5:		// MERV 5 Rated Filter
			MERV_n = 0;
			break;
		
		case 8:		// MERV 8 Rated Filter
			MERV_n = 1;
			break;

		case 11:	// MERV 11 Rated Filter
			MERV_n = 2;
			break;

		case 16:	// MERV 16 Rated Filter
			MERV_n = 3;
			break;
		}

		double A_qAH_high_array[4];
		double A_qAH_low_array[4];
		double A_wAH_high_array[4];
		double A_wAH_low_array[4];
		double A_DL_array[4];
		
		double k_qAH_array[3][4];
		double k_wAH_array[3][4];		
		double k_DL_array[3][4];			

		if(AHMotorType == 0) {
			// =================================== PSC Coefficients ========================================
			// Initial changes to AH airflow, power draw and duct static pressure due to filter
			// 1x4 array of filter airflow resistance coefficients
			A_qAH_high_array[0] = 1.00;		// MERV 5
			A_qAH_high_array[1] = 0.93;		// MERV 8
			A_qAH_high_array[2] = 0.85;		// MERV 11
			A_qAH_high_array[3] = 0.73;		// MERV 16

			A_qAH_low_array[0] = A_qAH_high_array[0];
			A_qAH_low_array[1] = A_qAH_high_array[1];
			A_qAH_low_array[2] = A_qAH_high_array[2];
			A_qAH_low_array[3] = A_qAH_high_array[3];
				
			// 1x4 array of filter air handler power draw coefficients
			A_wAH_high_array[0] = 1.00;		// MERV 5
			A_wAH_high_array[1] = 0.96;		// MERV 8
			A_wAH_high_array[2] = 0.91;		// MERV 11
			A_wAH_high_array[3] = 0.84;		// MERV 16

			// Coefficients for PSC motor are the same
			A_wAH_low_array[0] = A_wAH_high_array[0];
			A_wAH_low_array[1] = A_wAH_high_array[1];
			A_wAH_low_array[2] = A_wAH_high_array[2];
			A_wAH_low_array[3] = A_wAH_high_array[3];

			// 1x4 array of filter Duct Leakage change coefficients
			A_DL_array[0] = 1.00;		// MERV 5
			A_DL_array[1] = 1.22;		// MERV 8
			A_DL_array[2] = 1.44;		// MERV 11
			A_DL_array[3] = 1.80;		// MERV 16
		
			// Changes to AH airflow, power draw and duct static pressure due to filter loading
			// 3x4 array of filter loading coefficients k_qAH[row][column]
			// k_qAH is the change in airflow rate through the air handler
			k_qAH_array[0][0] = -1.0e-6;	// [f(low)][MERV_5]
			k_qAH_array[1][0] = -3.0e-6;	// [f(med)][MERV_5]
			k_qAH_array[2][0] = -10.0e-6;	// [f(high)][MERV_5]
			
			k_qAH_array[0][1] = -1.0e-6;	// [f(low)][MERV_8]
			k_qAH_array[1][1] = -3.0e-6;	// [f(med)][MERV_8]
			k_qAH_array[2][1] = -10.0e-6;	// [f(high)][MERV_8]
	
			k_qAH_array[0][2] = -1.0e-6;	// [f(low)][MERV_11]
			k_qAH_array[1][2] = -5.0e-6;	// [f(med)][MERV_11]
			k_qAH_array[2][2] = -20.0e-6;	// [f(high)][MERV_11]

			k_qAH_array[0][3] = -1.0e-6;	// [f(low)][MERV_16]
			k_qAH_array[1][3] = -16.0e-6;	// [f(med)][MERV_16]
			k_qAH_array[2][3] = -30.0e-6;	// [f(high)][MERV_16]

			// 3x4 array of filter loading coefficients k_wAH[row][column]
			// k_wAH is the change in power draw of the air handler
			k_wAH_array[0][0] = -1.0e-6;	// [f(low)][MERV_5]
			k_wAH_array[1][0] = -4.0e-6;	// [f(med)][MERV_5]
			k_wAH_array[2][0] = -9.0e-6;	// [f(high)][MERV_5]

			k_wAH_array[0][1] = -1.0e-6;	// [f(low)][MERV_8]
			k_wAH_array[1][1] = -4.0e-6;	// [f(med)][MERV_8]
			k_wAH_array[2][1] = -9.0e-6;	// [f(high)][MERV_8]

			k_wAH_array[0][2] = -1.0e-6;	// [f(low)][MERV_11]
			k_wAH_array[1][2] = -7.0e-6;	// [f(med)][MERV_11]
			k_wAH_array[2][2] = -35.0e-6;	// [f(high)][MERV_11]

			k_wAH_array[0][3] = -1.0e-6;	// [f(low)][MERV_16]
			k_wAH_array[1][3] = -9.0e-6;	// [f(med)][MERV_16]
			k_wAH_array[2][3] = -60.0e-6;	// [f(high)][MERV_16]

			// 3x4 array of filter loading coefficients k_DL[row][column]
			// k_DL is the change in Return Duct Leakage (retLF)
			k_DL_array[0][0] = 5.0e-6;	// [f(low)][MERV_5]
			k_DL_array[1][0] = 10.0e-6;	// [f(med)][MERV_5]
			k_DL_array[2][0] = 20.0e-6;	// [f(high)][MERV_5]

			k_DL_array[0][1] = 5.0e-6;	// [f(low)][MERV_8]
			k_DL_array[1][1] = 10.0e-6;	// [f(med)][MERV_8]
			k_DL_array[2][1] = 20.0e-6;	// [f(high)][MERV_8]

			k_DL_array[0][2] = 5.0e-6;	// [f(low)][MERV_11]
			k_DL_array[1][2] = 10.0e-6;	// [f(med)][MERV_11]
			k_DL_array[2][2] = 20.0e-6;	// [f(high)][MERV_11]

			k_DL_array[0][3] = 5.0e-6;	// [f(low)][MERV_16]
			k_DL_array[1][3] = 20.0e-6;	// [f(med)][MERV_16]
			k_DL_array[2][3] = 50.0e-6;	// [f(high)][MERV_16]
			
		} else {
			// =================================== BPM Coefficients ========================================
			// Initial changes to AH airflow, power draw and duct static pressure due to filter
			// 1x4 array of filter airflow resistance coefficients
			A_qAH_high_array[0] = 1.00; // 0.81;		// MERV 5
			A_qAH_high_array[1] = 1.00; // 0.71;		// MERV 8
			A_qAH_high_array[2] = 1.00; // 0.62;		// MERV 11
			A_qAH_high_array[3] = 1.00; // 0.51;		// MERV 16

			A_qAH_low_array[0] = 1.00;		// MERV 5
			A_qAH_low_array[1] = 1.00;		// MERV 8
			A_qAH_low_array[2] = 1.00;		// MERV 11
			A_qAH_low_array[3] = 1.00;		// MERV 16
				
			// 1x4 array of filter coefficients for air handler power draw in HIGH speed mode
			A_wAH_high_array[0] = 1.00; // 1.03;		// MERV 5
			A_wAH_high_array[1] = 1.10; // 1.06;		// MERV 8
			A_wAH_high_array[2] = 1.15; // 1.09;		// MERV 11
			A_wAH_high_array[3] = 1.20; // 1.13;		// MERV 16

			// 1x4 array of filter coefficients for air handler power draw in LOW speed mode
			A_wAH_low_array[0] = 1.00;		// MERV 5
			A_wAH_low_array[1] = 1.10; // 1.37;		// MERV 8
			A_wAH_low_array[2] = 1.15; // 1.74;		// MERV 11
			A_wAH_low_array[3] = 1.20; // 2.07;		// MERV 16

			// 1x4 array of filter Duct Leakage change coefficients
			A_DL_array[0] = 1.00;			// MERV 5
			A_DL_array[1] = 1.22;			// MERV 8
			A_DL_array[2] = 1.44;			// MERV 11
			A_DL_array[3] = 1.80;			// MERV 16
		
			// Changes to AH airflow, power draw and duct static pressure due to filter loading
			// 3x4 array of filter loading coefficients k_qAH[row][column]
			// k_qAH is the change in airflow rate through the air handler
			k_qAH_array[0][0] = 0.00;		// [f(low)][MERV_5]
			k_qAH_array[1][0] = 0.00;		// [f(med)][MERV_5]
			k_qAH_array[2][0] = 0.00;		// [f(high)][MERV_5]
			
			k_qAH_array[0][1] = 0.00;		// [f(low)][MERV_8]
			k_qAH_array[1][1] = 0.00;		// [f(med)][MERV_8]
			k_qAH_array[2][1] = 0.00;		// [f(high)][MERV_8]
	
			k_qAH_array[0][2] = 0.00;		// [f(low)][MERV_11]
			k_qAH_array[1][2] = 0.00;		// [f(med)][MERV_11]
			k_qAH_array[2][2] = 0.00;		// [f(high)][MERV_11]

			k_qAH_array[0][3] = 0.00;		// [f(low)][MERV_16]
			k_qAH_array[1][3] = 0.00;		// [f(med)][MERV_16]
			k_qAH_array[2][3] = 0.00;		// [f(high)][MERV_16]

			// 3x4 array of filter loading coefficients k_wAH[row][column]
			// k_wAH is the change in power draw of the air handler
			k_wAH_array[0][0] = 1.0e-6;		// [f(low)][MERV_5]
			k_wAH_array[1][0] = 2.5e-6;		// [f(med)][MERV_5]
			k_wAH_array[2][0] = 5.0e-6;		// [f(high)][MERV_5]

			k_wAH_array[0][1] = 1.0e-6;		// [f(low)][MERV_8]
			k_wAH_array[1][1] = 2.5e-6;		// [f(med)][MERV_8]
			k_wAH_array[2][1] = 5.0e-6;		// [f(high)][MERV_8]

			k_wAH_array[0][2] = 1.0e-6;		// [f(low)][MERV_11]
			k_wAH_array[1][2] = 2.5e-6;		// [f(med)][MERV_11]
			k_wAH_array[2][2] = 5.0e-6;		// [f(high)][MERV_11]

			k_wAH_array[0][3] = 1.0e-6;		// [f(low)][MERV_16]
			k_wAH_array[1][3] = 2.5e-6;		// [f(med)][MERV_16]
			k_wAH_array[2][3] = 5.0e-6;		// [f(high)][MERV_16]

			// 3x4 array of filter loading coefficients k_DL[row][column]
			// k_DL is the change in Return Duct Leakage (retLF)
			k_DL_array[0][0] = 5.0e-6;		// [f(low)][MERV_5]
			k_DL_array[1][0] = 10.0e-6;		// [f(med)][MERV_5]
			k_DL_array[2][0] = 20.0e-6;		// [f(high)][MERV_5]

			k_DL_array[0][1] = 5.0e-6;		// [f(low)][MERV_8]
			k_DL_array[1][1] = 10.0e-6;		// [f(med)][MERV_8]
			k_DL_array[2][1] = 20.0e-6;		// [f(high)][MERV_8]

			k_DL_array[0][2] = 5.0e-6;		// [f(low)][MERV_11]
			k_DL_array[1][2] = 10.0e-6;		// [f(med)][MERV_11]
			k_DL_array[2][2] = 20.0e-6;		// [f(high)][MERV_11]

			k_DL_array[0][3] = 5.0e-6;		// [f(low)][MERV_16]
			k_DL_array[1][3] = 20.0e-6;		// [f(med)][MERV_16]
			k_DL_array[2][3] = 50.0e-6;		// [f(high)][MERV_16]
		}

		// Set filter loading coefficients depending on MERV and loading choices
		if(qAH_heat0 > qAH_cool0) {	// If heating AH airflow is larger than cooling AH airflow
			qAH_low = qAH_cool0;									// Set cooling AH speed to lowest AH speed
			
			A_qAH_heat = A_qAH_high_array[MERV_n];					// AH airflow coefficient
			A_qAH_cool = A_qAH_low_array[MERV_n];					// AH airflow coefficient
			A_wAH_heat = A_wAH_high_array[MERV_n];					// AH power coefficient
			A_wAH_cool = A_wAH_low_array[MERV_n];					// AH power coefficient
			A_DL  = A_DL_array[MERV_n];								// Return duct leakage coefficient
			
			k_qAH = k_qAH_array[loadingRate][MERV_n];
			k_wAH = k_wAH_array[loadingRate][MERV_n];
			k_DL = k_DL_array[loadingRate][MERV_n];					// Return duct leakage gradient

		} else {					// If cooling AH airflow is larger than heating AH airflow (normal)
			qAH_low = qAH_heat0;									// Set heating AH speed to lowest AH speed
			
			A_qAH_cool = A_qAH_high_array[MERV_n];					// AH airflow coefficient
			A_qAH_heat = A_qAH_low_array[MERV_n];					// AH airflow coefficient
			A_wAH_cool = A_wAH_high_array[MERV_n];					// AH power coefficient
			A_wAH_heat = A_wAH_low_array[MERV_n];					// AH power coefficient
			A_DL  = A_DL_array[MERV_n];								// Return duct leakage coefficient
			
			k_qAH = k_qAH_array[loadingRate][MERV_n];
			k_wAH = k_wAH_array[loadingRate][MERV_n];
			k_DL = k_DL_array[loadingRate][MERV_n];					// Return duct leakage gradient
		}
	}

void f_CpTheta(double CP[4][4], int& windAngle, double* wallCp) {
	// this function takes Cps from a single wind angle perpendicular to the
	// upwind wall and finds Cps for all the walls for any wind angle
	
	double CPvar = 0;
	double Theta = 0;
	
	for(int i=0; i < 4; i++) {
		Theta = windAngle * M_PI / 180;
		
		if(i == 1) {
			Theta = Theta - M_PI;
		} else if(i == 2) {
			Theta = Theta - .5 * M_PI;
		} else if(i == 3) {
			Theta = Theta - 1.5 * M_PI;
		}
		
		CPvar = (CP[i][0] + CP[i][1]) * pow(pow(cos(Theta),2),.25);
		CPvar = CPvar + (CP[i][0] - CP[i][1]) * cos(Theta) * pow(abs(cos(Theta)),-.25);
		CPvar = CPvar + (CP[i][2] + CP[i][3]) * pow(pow(sin(Theta),2),2);
		CPvar = CPvar + (CP[i][2] - CP[i][3]) * sin(Theta);
		
		wallCp[i] = CPvar / 2;
	}
}

void f_flueFlow(double& tempHouse, double& flueShelterFactor, double& dPwind, double& dPtemp, double& eaveHeight, double& Pint, int& numFlues, flue_struct* flue, double& mFlue,
	double& airDensityOUT, double& airDensityIN, double& dPflue, double& tempOut, double& houseVolume, double& windPressureExp) {

		// calculates flow through the flue

		// dkm: internal variable to calculate external variable mFlue (because there may be more than one flue)
		double massflue = 0;
		double CpFlue;
		double fluePressureExp = 0.5;
		//double P = 0.14;	// Wind pressure coefficient of flue (0.14 for 

		for(int i=0; i < numFlues; i++) {
			//CpFlue = -.5 * pow(flueShelterFactor,2) * pow((flue[i].flueHeight / h),(2 * P));
			CpFlue = -.5 * pow((flue[i].flueHeight / eaveHeight),(2 * windPressureExp));

			if(flue[i].flueTemp == -99) {
					
				dPflue = Pint - dPtemp * flue[i].flueHeight + dPwind * pow(flueShelterFactor,2) * CpFlue;
				if(dPflue >= 0) { 	// flow in through flue
					massflue = massflue + (airDensityOUT * flue[i].flueC * pow((airTempRef / tempHouse),(3 * fluePressureExp - 2)) * pow(dPflue,fluePressureExp));
				} else {			// flow out through flue
					massflue = massflue + (-airDensityIN * flue[i].flueC * pow((airTempRef / tempOut),(3 * fluePressureExp - 2)) * pow(-dPflue,fluePressureExp));
				}
			} else {
				// for a heated flue:  driving pressure correction:
				dPflue = Pint - dPtemp * flue[i].flueHeight + dPwind * pow(flueShelterFactor,2) * CpFlue - g * airDensityIN * flue[i].flueHeight * (1 - tempHouse / flue[i].flueTemp);
				// density-viscosity correction:
				if(dPflue >= 0) {
					massflue = massflue + (pow((airTempRef / flue[i].flueTemp),(3 * fluePressureExp - 2)) * airDensityOUT * flue[i].flueC * pow(dPflue,fluePressureExp));
				} else {
					massflue = massflue + (pow(-(airTempRef / flue[i].flueTemp),(3 * fluePressureExp - 2)) * airDensityIN * flue[i].flueC * pow(-dPflue,fluePressureExp));
				}
			}
		}

		mFlue = massflue;
}


void f_floorFlow3(double& Cfloor, double& Cpfloor, double& dPwind, double& Pint,
	double& n, double& mFloor, double& airDensityOUT, double& airDensityIN, double& Hfloor, double& dPtemp) {

		// calculates flow through floor level leaks
		double dPfloor = Pint + Cpfloor * dPwind - Hfloor * dPtemp;

		if(dPfloor >= 0)
			mFloor = airDensityOUT * Cfloor * pow(dPfloor,n);
		else
			mFloor = -airDensityIN * Cfloor * pow(-dPfloor,n);
}

void f_ceilingFlow(int& AHflag, double& Patticint, double& eaveHeight, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& mCeiling, double& atticC, double& airDensityATTIC,
	double& airDensityIN, double& tempAttic, double& tempHouse, double& tempOut, double& airDensityOUT,
	double& mSupAHoff, double& mRetAHoff, double& supC, double& supn, double& retC, double& retn, double CCeiling) {

		double dPceil = Pint - Patticint - airDensityOUT * g * ((tempHouse - tempOut) / tempHouse - (tempAttic - tempOut) / tempAttic) * eaveHeight;

		if(dPceil >= 0) {
			mCeiling = airDensityATTIC * CCeiling * pow(dPceil,n);
			if(AHflag == 0) {
				mSupAHoff = airDensityATTIC * supC * pow(dPceil,supn);
				mRetAHoff = airDensityATTIC * retC * pow(dPceil,retn);
				if(atticC == 0) {
					mSupAHoff = 0;
					mRetAHoff = 0;
					mCeiling = 0;
				}
			} else {
				mSupAHoff = 0;
				mRetAHoff = 0;
			}
		} else {
			mCeiling = -airDensityIN * CCeiling * pow(-dPceil,n);
			if(AHflag == 0) {
				mSupAHoff = -airDensityIN * supC * pow(-dPceil,supn);
				mRetAHoff = -airDensityIN * retC * pow(-dPceil,retn);
				if(atticC == 0) {
					mCeiling = 0;
					mSupAHoff = 0;
					mRetAHoff = 0;
				}
			} else {
				mSupAHoff = 0;
				mRetAHoff = 0;
			}
		}
}

// calculates the flow through a wall
void f_wallFlow3(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& wallCp,
	double& n, double& Cwall, double& eaveHeight, double& Pint, double& dPtemp, double& dPwind,
	double& mWallIn, double& mWallOut, double& Hfloor) {
		
		double Hwall = eaveHeight - Hfloor;
		double dPwalltop = Pint + dPwind * wallCp - dPtemp * eaveHeight;
		double dPwallbottom = Pint + dPwind * wallCp - dPtemp * Hfloor;
		double Bo = f_neutralLevel(dPtemp, dPwind, Pint, wallCp, eaveHeight);

		if(tempHouse == tempOut) {
			if(dPwallbottom > 0) {
				mWallIn = airDensityOUT * Cwall * pow(dPwallbottom, n);
				mWallOut = 0;
			} else {
				mWallIn = 0;
				mWallOut = -airDensityIN * Cwall * pow(-dPwallbottom, n);
			}
		} else {
			if(tempHouse > tempOut) {
				if(Bo <= 0) {
					mWallIn = 0;
					mWallOut = airDensityIN * Cwall / Hwall / dPtemp / (n + 1) * (dPwalltop * pow(abs(dPwalltop), n) - dPwallbottom * pow(abs(dPwallbottom), n));
				} else if(Bo >= 1) {
					mWallIn = -airDensityOUT * Cwall / Hwall / dPtemp / (n + 1) * (dPwalltop * pow(abs(dPwalltop), n) - dPwallbottom * pow(abs(dPwallbottom), n));
					mWallOut = 0;
				} else {
					mWallIn = airDensityOUT * Cwall / Hwall / dPtemp / (n + 1) * dPwallbottom * pow(abs(dPwallbottom), n);
					mWallOut = airDensityIN * Cwall / Hwall / dPtemp / (n + 1) * dPwalltop * pow(abs(dPwalltop), n);
				}
			} else {
				if(Bo <= 0) {
					mWallIn = -airDensityOUT * Cwall / Hwall / dPtemp / (n + 1) * (dPwalltop * pow(abs(dPwalltop), n) - dPwallbottom * pow(abs(dPwallbottom), n));
					mWallOut = 0;
				} else if(Bo >= 1) {
					mWallIn = 0;
					mWallOut = airDensityIN * Cwall / Hwall / dPtemp / (n + 1) * (dPwalltop * pow(abs(dPwalltop), n) - dPwallbottom * pow(abs(dPwallbottom), n));
				} else {
					mWallIn = -airDensityOUT * Cwall / Hwall / dPtemp / (n + 1) * dPwalltop * pow(abs(dPwalltop), n);
					mWallOut = -airDensityIN * Cwall / Hwall / dPtemp / (n + 1) * dPwallbottom * pow(abs(dPwallbottom), n);
				}
			}
		}
}

void f_fanFlow(fan_struct& fan, double& airDensityOUT, double& airDensityIN) {
	
	// calculates flow through ventilation fans
	if(fan.q < 0)
		fan.m = airDensityIN * fan.q;
	else
		fan.m = airDensityOUT * fan.q;
}

void f_pipeFlow(double& airDensityOUT, double& airDensityIN, double& CP, double& dPwind, double& dPtemp,
	double& Pint, pipe_struct& Pipe, double& tempHouse, double& tempOut) {
		
		// calculates flow through pipes
		// changed on NOV 7 th 1990 so Pipe.A is Cpipe
		// changed on Nov 9th 1990 for density & viscosity variation of C
		// this is representeed by the temperature ratios

		/* WJNT: f_pipeFlow seems to do the same as f_flueFlow with the exception of using tempOut instead of tempHouse while Pipe.dP >= 0
		will consider removing- check with Iain */

		Pipe.dP = Pint - dPtemp * Pipe.h + dPwind * CP;

		if(Pipe.dP >= 0)
			Pipe.m = airDensityOUT * Pipe.A * pow((airTempRef / tempOut), (3 * Pipe.n - 2)) * pow(Pipe.dP, Pipe.n);
		else
			Pipe.m = -airDensityIN * Pipe.A * pow((airTempRef / tempHouse), (3 * Pipe.n - 2)) * pow(-Pipe.dP, Pipe.n);
}

void f_winDoorFlow(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& eaveHeight,
	double& wallCp, double& n, double& Pint, double& dPtemp, double& dPwind, winDoor_struct& winDoor) {

		// calculates flow through open doors or windows

		double dT;
		double Awindoor;
		double dummy;
		double dummy1;
		double dummy2;
		double Viscosity;
		double Kwindow;
		double Kwindownew;
		double Topcrit;
		double Bottomcrit;
		double Pwindow;
		double ReH;

		dT = tempHouse - tempOut;
		winDoor.dPbottom = 2 * (dPwind * wallCp + Pint - winDoor.Bottom * dPtemp) / airDensityOUT;
		winDoor.dPtop = 2 * (dPwind * wallCp + Pint - winDoor.Top * dPtemp) / airDensityOUT;
		double Bo = f_neutralLevel(dPtemp, dPwind, Pint, wallCp, eaveHeight);

		// winDoor.High and winDoor.Wide are not yet read in as inputs
		if(dT == 0) {
			Awindoor = winDoor.High * winDoor.Wide;
			dummy = Pint + dPwind * wallCp;
			if(dummy >= 0) {
				winDoor.mIN = .6 * Awindoor * sqrt(dummy * airDensityOUT * 2);
				winDoor.mOUT = 0;
			} else {
				winDoor.mIN = 0;
				winDoor.mOUT = -.6 * Awindoor * sqrt(-dummy * airDensityIN * 2);
			}
		} else {
			dummy1 = winDoor.dPbottom;
			dummy2 = winDoor.dPtop;
			if(Bo * eaveHeight <= winDoor.Bottom) {
				Kwindow = .6;
				dummy = dummy2 * sqrt(abs(dummy2)) - dummy1 * sqrt(abs(dummy1));
				if(dT > 0) {
					winDoor.mIN = 0;
					winDoor.mOUT = sqrt(airDensityIN * airDensityOUT) * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
				} else {
					winDoor.mIN = -airDensityOUT * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
					winDoor.mOUT = 0;
				}
			} else if(Bo * eaveHeight > winDoor.Top) {
				Kwindow = .6;
				dummy = dummy1 * sqrt(abs(dummy1)) - dummy2 * sqrt(abs(dummy2));
				if(dT > 0) {
					winDoor.mIN = airDensityOUT * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
					winDoor.mOUT = 0;
				} else {
					winDoor.mIN = 0;
					winDoor.mOUT = -sqrt(airDensityOUT * airDensityIN) * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
				}
			} else {
				Viscosity = .0000133 + .0000009 * ((tempOut + tempHouse) / 2 - C_TO_K);
				Kwindow = .4 + .0045 * dT;
				Topcrit = winDoor.Top - .1 * winDoor.High;
				Bottomcrit = winDoor.Bottom + .1 * winDoor.High;
				if(Bo * eaveHeight > Topcrit) {
					Pwindow = Topcrit * dPtemp - dPwind * wallCp;
				} else if(Bo * eaveHeight < Bottomcrit) {
					Pwindow = Bottomcrit * dPtemp - dPwind * wallCp;
				} else {
					Pwindow = Pint;
				}
				dummy1 = 2 * (dPwind * wallCp + Pwindow - winDoor.Bottom * dPtemp) / airDensityOUT;
				dummy2 = 2 * (dPwind * wallCp + Pwindow - winDoor.Top * dPtemp) / airDensityOUT;
				do {
					Kwindownew = Kwindow;
					if(dT > 0) {
						winDoor.mIN = airDensityOUT * Kwindownew * winDoor.Wide * tempHouse / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
						winDoor.mOUT = sqrt(airDensityIN * airDensityOUT) * Kwindownew * winDoor.Wide * tempHouse / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
					} else {
						winDoor.mIN = -airDensityOUT * Kwindownew * winDoor.Wide * tempHouse / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
						winDoor.mOUT = -sqrt(airDensityIN * airDensityOUT) * Kwindownew * winDoor.Wide * tempHouse / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
					}
					ReH = 2 * abs(winDoor.mIN / airDensityOUT + winDoor.mOUT / airDensityIN) / winDoor.Wide / Viscosity;
					Kwindow = (.3 + .6 * .00003 * ReH) / (1 + .00003 * ReH);
				} while(!(abs(Kwindownew - Kwindow) < .001));

				if(Bo * eaveHeight > Topcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * winDoor.High) * (eaveHeight * Bo - Topcrit) + Kwindow;
				} else if(Bo * eaveHeight < Bottomcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * winDoor.High) * (Bottomcrit - eaveHeight * Bo) + Kwindow;
				}
				if(dT > 0) {
					winDoor.mIN = airDensityOUT * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
					winDoor.mOUT = sqrt(airDensityIN * airDensityOUT) * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
				} else {
					winDoor.mIN = -airDensityOUT * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy2 * sqrt(abs(dummy2));
					winDoor.mOUT = -sqrt(airDensityIN * airDensityOUT) * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy1 * sqrt(abs(dummy1));
				}
			}
		}

		winDoor.m = winDoor.mIN + winDoor.mOUT;
}

void f_roofCpTheta(double* Cproof, int& windAngle, double* Cppitch, double& roofPitch) {
	
	double Theta;
	double F, C, S;

	for(int i=0; i < 4; i++) {
		Theta = windAngle * M_PI / 180;
		if(i == 1) {
			Theta = Theta - M_PI;
		} else if(i == 2) {
			Theta = Theta - .5 * M_PI;
		} else if(i == 3) {
			Theta = Theta - 1.5 * M_PI;
		}

		S = sin(Theta);
		C = cos(Theta);
		if(roofPitch < 28)
			F = C;
		else
			F = pow(C,5);

		Cppitch[i] = (Cproof[0] + Cproof[1]) * pow(C,2);
		Cppitch[i] = Cppitch[i] + (Cproof[0] - Cproof[1]) * F;
		Cppitch[i] = Cppitch[i] + (Cproof[2] + Cproof[3]) * pow(S,2);
		Cppitch[i] = Cppitch[i] + (Cproof[2] - Cproof[3]) * S;
		Cppitch[i] = Cppitch[i] / 2;
	}
}

// calculates the neutral level for a surface
double f_neutralLevel(double dPtemp, double dPwind, double Pint, double Cpr, double h) {
	if(dPtemp != 0)
		return((Pint + dPwind * Cpr) / dPtemp / h);
	else
		return(0);
}


// calculates the flow through the pitched section of the roof
void f_roofFlow(double& tempAttic, double& tempOut, double& airDensityATTIC, double& airDensityOUT, double& Cpr,
	double& atticPressureExp, double& Croof, double& roofPeakHeight, double& Patticint, double& dPtemp, double& dPwind,
	double& mRoofIn, double& mRoofOut, double& eaveHeight) {
		
		double Hroof = roofPeakHeight - eaveHeight;
		double dProoftop = Patticint + dPwind * Cpr - dPtemp * roofPeakHeight;
		double dProofbottom = Patticint + dPwind * Cpr - dPtemp * eaveHeight;
		double bRoof = f_neutralLevel(dPtemp, dPwind, Patticint, Cpr, roofPeakHeight);

		if(tempAttic == tempOut) {
			if(dProofbottom > 0) {
				mRoofIn = airDensityOUT * Croof * pow(dProofbottom, atticPressureExp);
				mRoofOut = 0;
			} else {
				mRoofIn = 0;
				mRoofOut = -airDensityATTIC * Croof * pow(-dProofbottom, atticPressureExp);
			}
		} else {
			if(tempAttic > tempOut) {
				if(bRoof <= eaveHeight / roofPeakHeight) {
					mRoofIn = 0;
					mRoofOut = airDensityATTIC * Croof / Hroof / dPtemp / (atticPressureExp + 1) * (dProoftop * pow(abs(dProoftop), atticPressureExp) - dProofbottom * pow(abs(dProofbottom), atticPressureExp));
				} else if(bRoof >= 1) {
					mRoofIn = -airDensityOUT * Croof / Hroof / dPtemp / (atticPressureExp + 1) * (dProoftop * pow(abs(dProoftop), atticPressureExp) - dProofbottom * pow(abs(dProofbottom), atticPressureExp));
					mRoofOut = 0;
				} else {
					mRoofIn = airDensityOUT * Croof / Hroof / dPtemp / (atticPressureExp + 1) * dProofbottom * pow(abs(dProofbottom), atticPressureExp);
					mRoofOut = airDensityATTIC * Croof / Hroof / dPtemp / (atticPressureExp + 1) * dProoftop * pow(abs(dProoftop), atticPressureExp);
				}
			} else {
				if(bRoof <= eaveHeight / roofPeakHeight) {
					mRoofIn = -airDensityOUT * Croof / Hroof / dPtemp / (atticPressureExp + 1) * (dProoftop * pow(abs(dProoftop), atticPressureExp) - dProofbottom * pow(abs(dProofbottom), atticPressureExp));
					mRoofOut = 0;
				} else if(bRoof >= 1) {
					mRoofIn = 0;
					mRoofOut = airDensityATTIC * Croof / Hroof / dPtemp / (atticPressureExp + 1) * (dProoftop * pow(abs(dProoftop), atticPressureExp) - dProofbottom * pow(abs(dProofbottom), atticPressureExp));
				} else {
					mRoofIn = -airDensityOUT * Croof / Hroof / dPtemp / (atticPressureExp + 1) * dProoftop * pow(abs(dProoftop), atticPressureExp);
					mRoofOut = -airDensityATTIC * Croof / Hroof / dPtemp / (atticPressureExp + 1) * dProofbottom * pow(abs(dProofbottom), atticPressureExp);
				}
			}
		}
}

void f_atticVentFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	atticVent_struct& atticVent, double& tempAttic, double& tempOut) {

		// calculates flow through attic roof vents
		atticVent.dP = Patticint - dPtemp * atticVent.h + dPwind * CP;

		if(atticVent.dP >= 0)
			atticVent.m = airDensityOUT * atticVent.A * pow((airTempRef / tempOut), (3 * atticVent.n - 2)) * pow(atticVent.dP, atticVent.n);
		else
			atticVent.m = -airDensityATTIC * atticVent.A * pow((airTempRef / tempAttic), (3 * atticVent.n - 2)) * pow(-atticVent.dP, atticVent.n);
}

void f_soffitFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	soffit_struct& soffit, double& soffitFraction, double& atticC, double& atticPressureExp, double& tempAttic, double& tempOut) {
		
		// calculates flow through attic soffit vents and Gable end vents	
		soffit.dP = Patticint - dPtemp * soffit.h + dPwind * CP;

		if(soffit.dP >= 0)
			soffit.m = airDensityOUT * soffitFraction * atticC * pow((airTempRef / tempOut), (3 * atticPressureExp - 2)) * pow(soffit.dP, atticPressureExp);
		else
			soffit.m = -airDensityATTIC * soffitFraction * atticC * pow((airTempRef / tempAttic), (3 * atticPressureExp - 2)) * pow(-soffit.dP, atticPressureExp);
}

void f_atticFanFlow(fan_struct& atticFan, double& airDensityOUT, double& airDensityATTIC) {
	
	// calculates flow through ventilation fans
	if(atticFan.q < 0)
		atticFan.m = airDensityATTIC * atticFan.q;
	else
		atticFan.m = airDensityOUT * atticFan.q;
}


/*
* heatTranCoef()
*
* Calculates heat transfer coefficient (hT) using combination of natural and forced convection from Ford
* From Walker (1993) eqn 3-41, 3-55
* @param tempi - surface temperature (deg K)
* @param tempa - air temperature (deg K)
* @param velocity - air velocity (m/s)
* @return heat transfer coefficient (W/m2K)
*/

double heatTranCoef(double tempi, double tempa, double velocity) {
	double hNatural = 3.2 * pow(abs(tempi - tempa), 1.0 / 3.0); 		// Natural convection from Ford
	double tFilm = (tempi + tempa) / 2;                					// film temperature
	double hForced = (18.192 - .0378 * tFilm) * pow(velocity, 0.8);   // forced convection from Ford
	return pow((pow(hNatural, 3) + pow(hForced, 3)), 1.0 / 3.0);      // combine using cube
}

/*
* radTranCoef()
*
* Calculates radiation heat transfer coefficient (hR) using linearized solution from Holman.
* From Walker (1993) eqn 3-23
* @param emissivity - surface emissivity
* @param tempi - surface temperature of first surface (deg K)
* @param tempj - surface temperature of second surface (deg K)
* @param viewFactor - view factor between surfaces Fi-j
* @param areaRatio - ratio of surface areas (Ai/Aj)
* @return heat transfer coefficient (W/m2K)
*/

double radTranCoef(double emissivity, double tempi, double tempj, double viewFactor, double areaRatio) {
	double rT = (1 - emissivity) / emissivity + 1 / viewFactor + (1 - emissivity) / emissivity * areaRatio;
	return SIGMA * (tempi + tempj) * (pow(tempi, 2) + pow(tempj, 2)) / rT;
}
