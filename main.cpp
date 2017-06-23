#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>		// RAD: so far used only for setprecission() in cmd output
#include <time.h>
#include <vector>
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include "functions.h"
#include "weather.h"
#include "psychro.h"
#include "equip.h"
#include "moisture.h"
#include "constants.h"
#include "config/config.h"

using namespace std;

/* Additions/changes to REGCAP (WJNT)
Output files = .rco (minute-by-minute) and .rc2 (annual summary)
Input files = .csv
Occupancy scheduling including calculation of occupied dose and exposure
Dynamic exhuast fan scheduling - reads in new input files
Economizer overides thermostat control to prevent heating after pre-cooling (unless you use the 7-day running average
swtich between heating and cooling operation - then the economizer only runs on cooling days i.e. hcFlag = 2)
Economizer operation opens pressure relief
Aeq calculation options for no infiltration, default 62.2 infiltration credit, or Addendum N weather factors
Terrain wind and shelter coefficients are selectable
Running average thermostat
Filter loading effects (new subroutine) - initial performance and loading performance impacts
Passive stacks/flues can be flow limited to the mechanical 62.2 rate

RIVEC:
	Exposure limit equation added
	New RIVEC control agorithm (fan.oper = 50)
	RIVEC algorithm can include infiltration credits
	Controller can close passive stacks (Hybrid applications)

New inputs:
	Bathroom scheduling flag (now redundant in C++)
	Furnace AFUE efficiency (AFUE)
	Number of bedrooms (numBedrooms)
	Number of stories (numStories)
	Addendum N weather factors (weatherFactor)
	
	Node 1 is the Attic Air
	Node 2 is the Inner North Sheathing
	Node 3 is the Outer North Sheathing
	Node 4 is the Inner South Sheathing
	Node 5 is the Outer South Sheathing
	Node 6 is all of the Wood (joists, trusses, etc.) lumped together
	Node 7 is the Ceiling of the House
	Node 8 is the Floor of the Attic
	Node 9 is the Inner Gable Wall (both lumped together)
	Node 10 is the Outer Gable Wall (both lumped together)
	Node 11 is the Return Duct Outer Surface
	Node 12 is the Return Duct Air
	Node 13 is The Mass of the House
	Node 14 is the Supply Duct Outer Surface
	Node 15 is the Supply Duct Air
	Node 16 is the House Air (all one zone)

*/

// ============================= FUNCTIONS ==============================================================

// Main function
int main(int argc, char *argv[], char* envp[])
{ 	
	// Read in batch file name from command line
	string batchFileName = "";
	string configFileName = "";
	if ( (argc <= 1) || (argv[argc-1] == NULL) || (argv[argc-1][0] == '-') ) {
		cerr << "usage: " << argv[0] << " batch_file" << endl;
      return(1);
   }
   else {
      batchFileName = argv[argc-1];
   }

	// Open batch File ============================================================================================
	ifstream batchFile(batchFileName); 
	if(!batchFile) { 
		cout << "Cannot open batch file: " << batchFileName << endl;
		return 1; 
	}

	// read in config file
	batchFile >> configFileName;
	Config config(configFileName, envp);
	
	// File paths
	string inPath = config.pString("inPath");
	string outPath = config.pString("outPath");
	string weatherPath = config.pString("weatherPath");
	string schedulePath = config.pString("schedulePath");
	
	// output file control
	bool printMoistureFile = config.pBool("printMoistureFile");
	bool printFilterFile = config.pBool("printFilterFile");
	bool printOutputFile = config.pBool("printOutputFile");
	
	// configuration vars
	double atticMCInit = config.pDouble("atticMCInit");			// initial moisture content of attic wood (fraction)
	double dhDeadBand = config.pDouble("dhDeadBand");				// Dehumidifier dead band (+/- %RH)
	double cCapAdjustTime = config.pDouble("cCapAdjustTime");	// First minute adjustment of cooling capacity (fraction)
	double moistureModel = config.pDouble("moistureModel");

	// Simulation Batch Timing
	time_t startTime, endTime;
	time(&startTime);
	string runStartTime = ctime(&startTime);
	
	int simNum = 0;
	string simName = "";

	cout << "REGCAP++ Building Simulation Tool LBNL" << endl;

	// Main loop on each input file =======================================================
	while(batchFile >> simName)  {		
		simNum++;
		string inputFileName = inPath + simName + ".in";
		string outputFileName = outPath + simName + ".rco";
		string moistureFileName = outPath + simName + ".hum";
		string filterFileName = outPath + simName + ".fil";
		string summaryFileName = outPath + simName + ".rc2";

		// Declare structures
		winDoor_struct winDoor[10] = {0};
		fan_struct fan[10] = {0};
		fan_struct atticFan[10] = {0};
		pipe_struct Pipe[10] = {0};		
		atticVent_struct atticVent[10] = {0};
		soffit_struct soffit[4] = {0};
		flue_struct flue[6] = {0};

		//Declare arrays
		double Sw[4];
		double floorFraction[4]; 		// Fraction of leak in floor below wall 1, 2, 3 and 4
		double wallFraction[4]; 		// Fraction of leak in wall 1, 2, 3 and 4
		double Swinit[4][361];		
		double mFloor[4] = {0,0,0,0};
		double soffitFraction[5];
		double wallCp[4] = {0,0,0,0};
		double mechVentPower;
		double b[ATTIC_NODES];
		double tempOld[ATTIC_NODES];
		double HR[5];		
		double heatThermostat[24];
		double coolThermostat[24];
		int occupied[2][24];	    // Used for setting which hours of the weekday/weekend the house is occupied (1) or vacant (0)

		// Input file variables
		string weatherFileName;
		string fanScheduleFileName;
		string tstatFileName;
		string occupancyFileName;
		string shelterFileName;
		double envC;					// Envelope leakage coefficient
		double envPressureExp;		// Envelope Pressure Exponent
		double windSpeedMultiplier; //Wind speed multiplier (G), for 62.2-2016 infiltration calcs
		double shelterFactor; 			//Shelter Factor, for 62.2-2016 infiltration calcs
		double stackCoef; 				//Stack coefficient, for 62.2-2016 infiltration calcs
		double windCoef; 				//Wind coefficient, for 62.2-2016 infiltration calcs
		double eaveHeight;			// Eave Height [m]
		double R;						// Ceiling Floor Leakage Sum
		double X;						// Ceiling Floor Leakage Difference
		int numFlues;					// Number of flues/chimneys/passive stacks
		double flueShelterFactor;	// Shelter factor at the top of the flue (1 if the flue is higher than surrounding obstacles
		int numPipes;					// Number of passive vents but appears to do much the same as flues
		double Hfloor;
		string rowOrIsolated;		// House in a row (R) or isolated (any string other than R)
		bool rowHouse;					// Flag that house is in a row
		double houseVolume;			// Conditioned volume of house (m3)
		double floorArea;				// Conditioned floor area (m2)
		double planArea;				// Footprint of house (m2)
		double storyHeight;			// Story height (m)
		//double houseLength;			// Long side of house (m) NOT USED IN CODE
		//double houseWidth;			// Short side of house (m) NOT USED IN CODE
		double uaWall;					// UA of opaque wall elements (walls and doors) (W/K)
		double uaFloor;				// UA of floor or slab (no solar gain, not used for cooling load) (W/K)
		double uaWindow;				// UA of windows for conductive gain (W/K)
		int numWinDoor;
		int numFans;
		double windowWE;
		double windowN;
		double windowS;
		double winShadingCoef;
		double ceilRval_heat;
		double ceilRval_cool;
		double latentLoad;
		double internalGains1;
		double atticVolume;
		double atticC;
		double atticPressureExp;
		int numAtticVents;
		double roofPitch;
		string roofPeakOrient;		// Roof peak orientation, D = perpendicular to front of house (Wall 1), P = parrallel to front of house
		bool roofPeakPerpendicular;
		double roofPeakHeight;
		int numAtticFans;
		double roofRval;
		int roofType;
		double ductLocation;			
		double supThickness;
		double retThickness;
		double supRval;
		double retRval;
		double supLF;					//Supply duct leakage fraction (e.g., 0.01 = 1% leakage).					
		double retLF;					//Return duct leakage fraction (e.g., 0.01 = 1% leakage)
		double supLength;
		double retLength;
		double supDiameter;
		double retDiameter;
		double qAH_cool0;				// Cooling Air Handler air flow (m^3/s)
		double qAH_heat0;				// Heating Air Handler air flow (m^3/s)
		double supn;
		double retn;
		double supC;					// Supply leak flow coefficient
		double retC;					// Return leak flow coefficient
		//double buried;
		double capacityraw;
		double capacityari;
		double EERari;
		double hcapacity;				// Heating capacity [kBtu/h]
		double fanPower_heating0;	// Heating fan power [W]
		double fanPower_cooling0;	// Cooling fan power [W]
		double charge;
		double AFUE;				// Annual Fuel Utilization Efficiency for the furnace
		//int bathroomSchedule;	// Bathroom schedule file to use (1, 2 or 3)
		int numBedrooms;			// Number of bedrooms (for 62.2 target ventilation calculation)
		int numStories;			// Number of stories in the building (for Nomalized Leakage calculation)
		double weatherFactor;	// Weather Factor (w) (for infiltration calculation from ASHRAE 136)
		int terrain;				// 1 = large city centres, 2 = urban and suburban, 3 = open terrain, 4 = open sea
		int Crawl;					// Is there a crawlspace - 1=crawlspace (use first floorFraction), 0=not a crawlspace (use all floorFractions)
		double HRV_ASE;			// Apparent Sensible Effectiveness of HRV unit
		double ERV_SRE;			// Sensible Recovery Efficiency of ERV unit. SRE and TRE based upon averages from ERV units in HVI directory, as of 5/2015.
		double ERV_TRE;			// Total Recovery Efficiency of ERV unit, includes humidity transfer for moisture subroutine		
		double Aeq;				// 62.2-2016 Total Ventilation Rate, hr-1. 
		int InfCalc; 		//Index value, 0 = Annual infiltration rate for exposure calcs, 1 = Real-time infiltration estimates for exposure calcs.
		// Inputs to set filter type and loading rate
		int filterLoadingFlag;	// Filter loading flag = 0 (OFF) or 1 (ON)
		int MERV;					// MERV rating of filter (may currently be set to 5, 8, 11 or 16)
		int loadingRate;			// loading rate of filter, f (0,1,2) = (low,med,high)
		int AHMotorType;				// BPM (1) or PSC (0) air handler motor
		// Inputs to set humidity control
		int rivecFlagInd;		// Indicator variable that instructs a fan code 13 or 17 to be run by RIVEC controls. 1= yes, 0=no. Brennan.
		int OccContType;		// Type of occupancy control to be used in the RIVEC calculations. 1,2,...n Brennan.
		int AuxFanIndex; 	// //Index value that determines if auxiliary fans (dryer, kitchen and bath fans) are counted towards RIVEC relative exposure calculations
		double wCutoff;			// Humidity Ratio cut-off calculated as some percentile value for the climate zone. Brennan.
		double wDiffMaxNeg;		// Maximum average indoor-outdoor humidity differene, when wIn < wOut. Climate zone average.
		double wDiffMaxPos;		// Maximum average indoor-outdoor humidity differene, when wIn > wOut. Climate zone average.
		double W25[12]; 			// 25th percentiles for each month of the year, per TMY3
		double W75[12]; 			// 75th percentiles for each month of the year, per TMY3
		int FirstCut;				// Monthly Indexes assigned based on climate zone
		int SecondCut; 			// Monthly Indexes assigned based on climate zone
		double doseTarget;		// Targeted dose value
		double HiDose;				// Variable high dose value for real-time humidity control, based on worst-case large, low-occupancy home. 
		int HiMonths[3];
		int LowMonths[3];
		double HiMonthDose;
		double LowMonthDose;
		double dhCapacity;		// Dehumidifier capacity (pints/day)
		double dhEnergyFactor;	// Dehumidifier energy factor (L/kWh)
		double dhSetPoint;		// Dehumidifier set point (%RH)
		string endOfFile;
		
		// Zeroing the variables to create the sums for the .ou2 file
		long int minuteYear = 1;
		int endrunon = 0;
		double Mcoil = 0;
		double SHR = 0;		
		double meanOutsideTemp = 0;
		double meanAtticTemp = 0;
		double meanHouseTemp = 0;
		double meanHouseACH = 0;
		double meanFlueACH = 0;
		double gasTherm = 0;
		double AH_kWh = 0;
		double compressor_kWh = 0;
		double mechVent_kWh = 0;
		double furnace_kWh = 0;
		double dehumidifier_kWh = 0;
		double RHind60 = 0; //index value (0 or 1) if RHhouse > 60 
		double RHtot60 = 0; //cumulative sum of index value (0 or 1) if RHhouse > 60 
		double RHexcAnnual60 = 0; //annual fraction of the year where RHhouse > 60
		double RHind70 = 0; //index value (0 or 1) if RHhouse > 70 
		double RHtot70 = 0; //cumulative sum of index value (0 or 1) if RHhouse > 70
		double RHexcAnnual70 = 0; //annual fraction of the year where RHhouse > 70
		double RHHouse = 50, RHAttic = 50;  
		double latitude;
		double longitude;
		int timeZone;
		double altitude;

		// Open moisture output file
		ofstream moistureFile;
		if(printMoistureFile) {
			moistureFile.open(moistureFileName);
			if(!moistureFile) { 
				cout << "Cannot open moisture file: " << moistureFileName << endl;
				return 1; 
			}

			if(moistureModel == 1)
				for(int i=0; i<6; i++) {
					moistureFile << "MC" << i << "\tmTotal" << i << "\t";
					}
			moistureFile << "HROut\tHRHouse\tHRAttic\tHRSupply\tHRReturn\t";
			moistureFile << "RHHouse\tRHAttic\tTempHouse\ttempAttic" << endl;
		}

		// [START] Read in Building Inputs =========================================================================================================================
		ifstream buildingFile(inputFileName); 
		if(!buildingFile) { 
			cout << "Cannot open input file: " << inputFileName << endl;
			return 1; 
		}

		buildingFile >> weatherFileName;
		weatherFileName = weatherPath + weatherFileName;
		buildingFile >> fanScheduleFileName;
		fanScheduleFileName = schedulePath + fanScheduleFileName;
		buildingFile >> tstatFileName;
		tstatFileName = schedulePath + tstatFileName;
		buildingFile >> occupancyFileName;
		occupancyFileName = schedulePath + occupancyFileName;
		buildingFile >> shelterFileName;
		shelterFileName = schedulePath + shelterFileName;

		buildingFile >> envC;
		buildingFile >> envPressureExp;
		buildingFile >> windSpeedMultiplier;
		buildingFile >> shelterFactor;
		buildingFile >> stackCoef;
		buildingFile >> windCoef;
		buildingFile >> eaveHeight;
		buildingFile >> R;
		buildingFile >> X;
		buildingFile >> numFlues;
		for(int i=0; i < numFlues; i++) {
			buildingFile >> flue[i].flueC;
			buildingFile >> flue[i].flueHeight;
			buildingFile >> flue[i].flueTemp;
		}

		// =========================== Leakage Inputs ==============================
		for(int i=0; i < 4; i++) {
			buildingFile >> wallFraction[i];
		}
		for(int i=0; i < 4; i++) {
			buildingFile >> floorFraction[i];
		}

		buildingFile >> flueShelterFactor;
		buildingFile >> numPipes;
		for(int i=0; i < numPipes; i++) {
			buildingFile >> Pipe[i].wall;
			buildingFile >> Pipe[i].h;
			buildingFile >> Pipe[i].A;
			buildingFile >> Pipe[i].n;
			buildingFile >> Pipe[i].Swf;
			buildingFile >> Pipe[i].Swoff;
		}

		// =========================== Building Inputs =============================
		buildingFile >> Hfloor;
		buildingFile >> rowOrIsolated;
		buildingFile >> houseVolume;
		buildingFile >> floorArea;
		buildingFile >> planArea;
		buildingFile >> storyHeight;
		//buildingFile >> houseLength;
		//buildingFile >> houseWidth;
		buildingFile >> uaWall;
		buildingFile >> uaFloor;
		buildingFile >> uaWindow;

		// ====================== Venting Inputs (Windows/Doors)====================
		buildingFile >> numWinDoor;
		// These are not currently used in the Excel input generating spreadsheet
		// Add them back in if needed
		/*for(int i=0; i < numWinDoor; i++) {
			buildingFile >> winDoor[i].wall;
			buildingFile >> winDoor[i].High;
			buildingFile >> winDoor[i].Wide;
			buildingFile >> winDoor[i].Top;
			buildingFile >> winDoor[i].Bottom;
		}*/

		// ================== Mechanical Venting Inputs (Fans etc) =================
		buildingFile >> numFans;
		for(int i=0;  i < numFans; i++) {
			buildingFile >> fan[i].power;
			buildingFile >> fan[i].q;
			buildingFile >> fan[i].oper;
			fan[i].on = 0;
		}

		buildingFile >> windowWE;
		buildingFile >> windowN;
		buildingFile >> windowS;
		buildingFile >> winShadingCoef;
		buildingFile >> ceilRval_heat;
		buildingFile >> ceilRval_cool;
		buildingFile >> latentLoad;
		buildingFile >> internalGains1;

		// =========================== Attic Inputs ================================
		buildingFile >> atticVolume;
		buildingFile >> atticC;
		buildingFile >> atticPressureExp;
		for(int i=0; i < 5; i++) {
			buildingFile >> soffitFraction[i];
		}
		for(int i=0; i < 4; i++) {
			buildingFile >> soffit[i].h;
		}

		// =========================== Attic Vent Inputs ===========================
		buildingFile >> numAtticVents;
		for(int i=0; i < numAtticVents; i++) {
			buildingFile >> atticVent[i].wall;
			buildingFile >> atticVent[i].h;
			buildingFile >> atticVent[i].A;
			buildingFile >> atticVent[i].n;
		}

		// =========================== Roof Inputs =================================
		buildingFile >> roofPitch;
		buildingFile >> roofPeakOrient;
		buildingFile >> roofPeakHeight;
		buildingFile >> numAtticFans;
		for(int i = 0; i < numAtticFans; i++) {
			buildingFile >> atticFan[i].power;
			buildingFile >> atticFan[i].q;
			buildingFile >> atticFan[i].oper;
		}

		buildingFile >> roofRval;
		buildingFile >> roofType;

		// =========================== Duct Inputs =================================
		buildingFile >> ductLocation;
		buildingFile >> supThickness;
		buildingFile >> retThickness;
		buildingFile >> supRval;
		buildingFile >> retRval;
		buildingFile >> supLF;
		buildingFile >> retLF;
		buildingFile >> supLength;
		buildingFile >> retLength;
		buildingFile >> supDiameter;
		buildingFile >> retDiameter;
		buildingFile >> qAH_cool0;
		buildingFile >> qAH_heat0;
		buildingFile >> supn;
		buildingFile >> retn;
		buildingFile >> supC;
		buildingFile >> retC;
		//buildingFile >> buried; // NOTE: The buried variable is not used but left in code to continue proper file navigation

		// =========================== Equipment Inputs ============================
		buildingFile >> capacityraw;
		buildingFile >> capacityari;
		buildingFile >> EERari;
		buildingFile >> hcapacity;
		buildingFile >> fanPower_heating0;
		buildingFile >> fanPower_cooling0;
		buildingFile >> charge;
		buildingFile >> AFUE;
		//buildingFile >> bathroomSchedule;
		buildingFile >> numBedrooms;
		buildingFile >> numStories;
		buildingFile >> weatherFactor;
		buildingFile >> terrain;
		buildingFile >> Aeq;
		buildingFile >> InfCalc;
		buildingFile >> Crawl;
		buildingFile >> HRV_ASE;
		buildingFile >> ERV_SRE;
		buildingFile >> ERV_TRE;
		// =========================== Filter Inputs ============================
		buildingFile >> filterLoadingFlag;
		buildingFile >> MERV;
		buildingFile >> loadingRate;
		buildingFile >> AHMotorType;
		buildingFile >> rivecFlagInd;
		//The variable from here down were added by Brennan as part of the Smart Ventilation Humidity Control project
		buildingFile >> OccContType; //Now the OccContType for SVC_Occupancy.
		buildingFile >> AuxFanIndex; //Index value that determines if auxiliary fans (dryer, kitchen and bath fans) are counted towards RIVEC relative exposure calculations
// 		if(HumContType > 0) {
// 			buildingFile >> wCutoff;
// 			buildingFile >> wDiffMaxNeg;
// 			buildingFile >> wDiffMaxPos;
// 			for(int i = 0; i < 12; i++) {
// 				buildingFile >> W25[i];
// 			}
// 			for(int i = 0; i < 12; i++) {
// 				buildingFile >> W75[i];
// 			}
// 			buildingFile >> FirstCut;
// 			buildingFile >> SecondCut;
// 			buildingFile >> doseTarget;
// 			buildingFile >> HiDose;
// 			buildingFile >> HiMonths[0] >> HiMonths[1] >> HiMonths[2];
// 			buildingFile >> LowMonths[0] >> LowMonths[1] >> LowMonths[2];
// 			buildingFile >> HiMonthDose;
// 			buildingFile >> LowMonthDose;
// 		}
		buildingFile >> dhCapacity;	
		buildingFile >> dhEnergyFactor;
		buildingFile >> dhSetPoint;	
		
		buildingFile >> endOfFile;
		if(endOfFile != "E_O_F") {
			cout << "Error in input file. Last line: >>" << endOfFile << "<<" << endl;
			return 1;
		}	

		buildingFile.close();

		// [END] Read in Building Inputs ============================================================================================================================================

		// Read in Thermostat Settings ==================================================================
		ifstream tstatFile(tstatFileName); 
		if(!tstatFile) { 
			cout << "Cannot open thermostat file: " << tstatFileName << endl;
			return 1; 
		}
		string header;
		getline(tstatFile,header);
		for(int h = 0; h < 24; h++) {
			double heatT;
			double coolT;
			tstatFile >> heatT >> coolT;
			heatThermostat[h] = 273.15 + (heatT - 32) * 5.0 / 9.0;
			coolThermostat[h] = 273.15 + (coolT - 32) * 5.0 / 9.0;
		}
		tstatFile.close();

		// Read in Occupancy Settings ==================================================================
		ifstream occupancyFile(occupancyFileName); 
		if(!occupancyFile) { 
			cout << "Cannot open occupancy file: " << occupancyFileName << endl;
			return 1; 
		}
		getline(occupancyFile,header);
		for(int h = 0; h < 24; h++) {
			occupancyFile >> occupied[0][h] >> occupied[1][h];
		}
		occupancyFile.close();
		

		// Leakage fractions
		double leakFracCeil = (R + X) / 2 ;					// Fraction of leakage in the ceiling
		double leakFracFloor = (R - X) / 2;					// Fraction of leakage in the floor
		double leakFracWall = 1 - leakFracFloor - leakFracCeil;	// Fraction of leakage in the walls

		rowHouse = (rowOrIsolated.compare("R") == 0);
		roofPeakPerpendicular = (roofPeakOrient.compare("D") == 0);

		double supCp = 753.624;										// Specific heat capacity of steel [j/kg/K]
		double suprho = 16.018 * 2;									// Supply duct density. The factor of two represents the plastic and sprical [kg/m^3]
		double retCp = 753.624;										// Specific heat capacity of steel [j/kg/K]
		double retrho = 16.018 * 2;									// Return duct density [kg/m^3]
		double supArea = (supDiameter + 2 * supThickness) * M_PI * supLength;		// Surface area of supply ducts [m2]
		double retArea = (retDiameter + 2 * retThickness) * M_PI * retLength;		// Surface area of return ducts [m2]
		double supVolume = (pow(supDiameter, 2) * M_PI / 4) * supLength;			// Volume of supply ducts [m3]
		double retVolume = (pow(retDiameter, 2) * M_PI / 4) * retLength;			// Volume of return ducts [m3]
		hcapacity = hcapacity * .29307107 * 1000 * AFUE;			// Heating capacity of furnace converted from kBtu/hr to Watts and with AFUE adjustment
		double MWha = .5 * floorArea / 186;							// MWha is the moisture transport coefficient that scales with floor area (to scale with surface area of moisture)
					
		// [START] Filter Loading ==================================================================================

		int filterChanges = 0;			// Number of filters used throughout the year
		double qAH_low = 0;				// Lowest speed of AH (set in sub_filterLoading)

		// Set fan power and airflow rates to initial (0) input values for filter loading
		double qAH_heat = qAH_heat0;
		double qAH_cool = qAH_cool0;
		double fanPower_heating = fanPower_heating0;
		double fanPower_cooling = fanPower_cooling0;

		double massFilter_cumulative = 0;	// Cumulative mass that has flown through the filter
		double massAH_cumulative = 0;			// Cumulative mass that has flown through the AH

		// Filter loading coefficients (initialization- attributed values in sub_filterLoading)
		double A_qAH_heat = 0;		// Initial AH airflow change (in heating mode) from installing filter [%]
		double A_qAH_cool = 0;		// Initial AH airflow change (in cooling mode) from installing filter [%]
		double A_wAH_heat = 0;		// Initial AH power change (in heating mode) from installing filter [%]
		double A_wAH_cool = 0;		// Initial AH power change (in cooling mode) from installing filter [%]
		double A_DL = 0;				// Intial return duct leakage change from installing filter [%]
		
		double k_qAH = 0;				// Gradual change in AH airflow from filter loading [% per 10^6kg of air mass through filter]
		double k_wAH = 0;				// Gradual change in AH power from filter loading [% per 10^6kg of air mass through filter]
		double k_DL = 0;				// Gradual change in return duct leakage from filter loading [% per 10^6kg of air mass through filter]
		
		// Open filter loading file
		ofstream filterFile;
		if(printFilterFile) {
			filterFile.open(filterFileName);
			if(!filterFile) { 
				cout << "Cannot open filter file: " << filterFileName << endl;
				return 1; 
			}
			filterFile << "mAH_cumu\tqAH\twAH\tretLF" << endl;
		}
		
		// Filter loading coefficients are in the sub_filterLoading sub routine
		if(filterLoadingFlag == 1)
			sub_filterLoading(MERV, loadingRate, AHMotorType, A_qAH_heat, A_qAH_cool, A_wAH_heat, A_wAH_cool, A_DL, k_qAH, k_wAH, k_DL, qAH_heat0, qAH_cool0, qAH_low);
		// [END] Filter Loading ====================================================================================

		// Cooling capacity air flow correction term (using CFM)
		double qAH_cfm = qAH_cool / .0004719;
		double qAHcorr = 1.62 - .62 * qAH_cfm / (400 * capacityraw) + .647 * log(qAH_cfm / (400 * capacityraw));	

		// For RIVEC calculations
		int peakFlag = 0;						// Peak flag (0 or 1) prevents two periods during same day
		int rivecFlag = 0;					// Dose controlled ventilation (RIVEC) flag 0 = off, 1 = use rivec (mainly for HRV/ERV control)
		double qRivec = 1.0;					// Airflow rate of RIVEC fan for max allowed dose and exposure [L/s]
		// int numFluesActual = numFlues;	// numFluesActual used to remember the number of flues for when RIVEC blocks them off as part of a hybrid ventilation strategy

		// For Economizer calculations
		int economizerUsed = 0;			// 1 = use economizer and change C (house leakage) for pressure relief while economizer is running (changes automatically)
		double Coriginal = 0;			// House leakage before pressure relief opens for economizer
		double ELAeconomizer = 0;		// Size of leakage increase for economizer pressure relief
		double Ceconomizer = 0;			// Total leakage of house inncluding extra economizer pressure relief (only while economizer is running)

		for(int i=0; i < numFans; i++) {
			if(fan[i].oper == 21 || fan[i].oper == 22)
				economizerUsed = 1;		// Lets REGCAP know that an economizer is being used
				//economizerUsed = 0;		// Force economizer to be off
			if(fan[i].oper == 30 || fan[i].oper == 31 || fan[i].oper == 50 || fan[i].oper == 51) { // RIVEC fans
				qRivec = -1 * fan[i].q	* 1000;
				rivecFlag = 1;
			}
		}

		for(int i=0; i < numFans; i++) {
			if(fan[i].oper == 5 || fan[i].oper == 16) {	// If using an HRV (5 or 16) with RIVEC. Avoids a divide by zero
				qRivec = -1 * fan[i].q * 1000;
			}
			if((fan[i].oper == 13 && rivecFlagInd == 1) || (fan[i].oper == 17 && rivecFlagInd == 1)){ //If using CFIS or ERV+AHU with RIVEC control.
					rivecFlag = 1;
			}
		}

		// Presure Relief for Economizer, increase envelope leakage C while economizer is operating
		if(economizerUsed == 1) {
			Coriginal = envC;
			if(qAH_cool >= qAH_heat) {			// Dependent on the largest of the heating/cooling AH fan power
				// sized to 2Pa of pressure while economizer running
				ELAeconomizer = qAH_cool * (sqrt(airDensityRef / (2 * 2)) / 1);
				Ceconomizer = 1 * ELAeconomizer * sqrt(2 / airDensityRef) * pow(4, (.5 - .65));
			} else {
				ELAeconomizer = qAH_heat * (sqrt(airDensityRef / (2 * 2)) / 1);
				Ceconomizer = 1 * ELAeconomizer * sqrt(2 / airDensityRef) * pow(4, (.5 - .65));
			}
		}

		// [START] Read in Shelter Values ================================================================================================
		// (reading urban shelter values from a data file: Bshelter.dat)
		// Computed for the houses at AHHRF for every degree of wind angle

		ifstream shelterFile(shelterFileName); 
		if(!shelterFile) { 
			cout << "Cannot open shelter file: " << shelterFileName << endl;
			return 1; 
		}

		// FF: This angle variable is being overwritten to read Swinit in proper ductLocation per FOR iteration. Not used in code.
		double angle;
		for(int i=0; i < 361; i++) {
			shelterFile >> angle >> Swinit[0][i] >> Swinit[1][i] >> Swinit[2][i] >> Swinit[3][i];
		}

		shelterFile.close();
		// [END] Read in Shelter Values ================================================================================================

		// [START] Terrain ============================================================================================
		// Terrain where the house is located (for wind shelter etc.) See ASHRAE Fundamentals 2009 F24.3
		double windPressureExp;		// Power law exponent of the wind speed profile at the building site
		double layerThickness;		// Atmospheric boundary layer thickness [m]

		switch (terrain) {
		case 1:
			// 1. Large city centres, at least 50% of buildings are higher than 25m
			windPressureExp = 0.33;
			layerThickness = 460;
			break;

		case 2:
			// 2. Urban and suburban areas, wooded areas
			windPressureExp = 0.22;
			layerThickness = 370;
			break;

		case 3:
			// 3. Open terrain with scattered obsructions e.g. meteorlogical stations
			windPressureExp = 0.14;
			layerThickness = 270;
			break;
		case 4:
			// 4. Completely flat e.g. open water
			windPressureExp = 0.10;
			layerThickness = 210;
			break;
		}

		// Wind speed correction to adjust met wind speed to that at building eaves height
		const double metExponent = 0.14;		// Power law exponent of the wind speed profile at the met station
		const double metThickness = 270;		// Atmospheric boundary layer thickness at met station [m]
		const double metHeight = 10;			// Height of met station wind measurements [m]
		double windSpeedCorrection = pow((metThickness / metHeight), metExponent) * pow((eaveHeight / layerThickness), windPressureExp);
		//double windSpeedCorrection = pow((80/ 10), .15) * pow((h / 80), .3);		// The old wind speed correction
		// [END] Terrain ============================================================================================

		// AL4 is used to estimate flow velocities in the attic
		double AL4 = atticC * sqrt(airDensityRef / 2) * pow(4, (atticPressureExp - .5));

		// New char velocity for unvented attics
		if(AL4 == 0)
			AL4 = envC * leakFracCeil * sqrt(airDensityRef / 2) * pow(4, (atticPressureExp - .5));

		// ================= CREATE OUTPUT FILE =================================================
		ofstream outputFile;
		if(printOutputFile) {
			outputFile.open(outputFileName); 
			if(!outputFile) { 
				cout << "Cannot open output file: " << outputFileName << endl;
				return 1; 
			}
			outputFile << "Time\tMin\twindSpeed\ttempOut\ttempHouse\tsetpoint\ttempAttic\ttempSupply\ttempReturn\tAHflag\tAHpower\tcompressPower\tmechVentPower\tHR\tSHR\tMcoil\thousePress\tQhouse\tACH\tACHflue\tventSum\tnonRivecVentSum\tfan1\tfan2\tfan3\tfan4\tfan5\tfan6\tfan7\trivecOn\trelExp\trelDose\toccupied\tHROUT\tHRattic\tHRreturn\tHRsupply\tHRhouse\tHRmaterials\tRHhouse\tRHind60\tRHind70\tHumidityIndex\tDHcondensate\tPollutantConc" << endl; 
		}
		
		// 5 HR nodes
		// Node 1 is attic air (Node HR[0])
		// Node 2 is return air (Node HR[1])
		// Node 3 is supply air (Node HR[2])
		// Node 4 is house air (Node HR[3])
		// Node 5 is house materials that interact with house air only (Node HR[4])

		// ================== OPEN WEATHER FILE FOR INPUT ========================================
		ifstream weatherFile(weatherFileName);
		if(!weatherFile) { 
			cout << "Cannot open weather file: " << weatherFileName << endl;
			return 1; 
		}

		// Determine weather file type and read in header
		CSVRow row;					// row of data to read from TMY3 csv
		weatherData begin, end;	// hourly weather data
		int weatherFileType;
		weatherFile >> row;
		if(row.size() > 2 && row.size() < 10) {			//TMY3, was >2
			//double siteID = row[0];
			//string siteName = row[1];
			//string State = row[2];
			timeZone = row[3];
			latitude = row[4];
			longitude = row[5];
			altitude = row[6];
			cout << "ID=" << row[0] << " TZ=" << timeZone << " lat=" << latitude << " long=" << longitude << " alt=" << altitude << endl;
			cout << "Using TMY3 weather data file." << endl;
			weatherFile >> row;					// drop header
			begin = readTMY3(weatherFile);	// duplicate first hour for interpolation of 0->1
			end = begin;
			weatherFileType = 1;
		}
		
		else if(row.size() == 10) {			//EnergyPlus weather file EPW
			//double siteID = row[0];
			//string siteName = row[1];
			//string State = row[2];
			timeZone = row[8];
			latitude = row[6];
			longitude = row[7];
			altitude = row[9];
			cout << "ID=" << row[5] << " TZ=" << timeZone << " lat=" << latitude << " long=" << longitude << " alt=" << altitude << endl;
			cout << "Using Energy Plus weather data file." << endl;
 			for(int i = 0; i < 7; ++i) {	// drop rest of header lines
			   weatherFile >> row;
			   }
			begin = readEPW(weatherFile);	// duplicate first hour for interpolation of 0->1
			end = begin;
			weatherFileType = 2;
		}
		
		else {							// 1 minute data
			weatherFile.close();
			weatherFile.open(weatherFileName);
			weatherFile >> latitude >> longitude >> timeZone >> altitude;
			weatherFileType = 0;
			cout << "1 minute:" << "Lat:" << latitude << " Alt:" << altitude << endl;
			cout << "Using REGCAP weather data file." << endl;
		}
		latitude = M_PI * latitude / 180.0;					// convert to radians
		// set 1 = air handler off, and house air temp below setpoint minus 0.5 degrees
		// set 0 = air handler off, but house air temp above setpoint minus 0.5 degrees
		int set = 0;

		int AHflag = 0;			// Air Handler Flag (0/1/2 = OFF/HEATING MODE/COOLING MODE, 100 = fan on for venting, 102 = heating cool down for 1 minute at end of heating cycle)
		int AHflagPrev = 0;
		int econoFlag = 0;		// Economizer Flag (0/1 = OFF/ON)

		// ---------------- Attic and duct system heat transfer nodes --------------------
		for(int k=0; k < ATTIC_NODES; k++) {
			tempOld[k] = airTempRef;
		}
		tempOld[2] = 278;
		tempOld[4] = 278;

		double tempAttic = tempOld[0];
		double tempReturn = tempOld[11];
		double tempSupply = tempOld[14];
		double tempHouse = tempOld[15];

		// Setting initial values of air mass for moisture balance:
		double M1 = atticVolume * airDensityRef;		// Mass of attic air
		double M12 = retVolume * airDensityRef;		// Mass of return air
		double M15 = supVolume * airDensityRef;		// Mass of supply air
		double M16 = houseVolume * airDensityRef;		// Mass of house air
		double Mw5 = 60 * floorArea;						// Active mass containing moisture in the house (empirical)

		// Like the mass transport coefficient, the active mass for moisture scales with floor area

		// The following are defined for RIVEC ==============================================================================================================
		// Start and end times for the RIVEC base, peak and recovery periods [h]
		int baseStart = 0;
		int baseEnd = 0;
		int peakStart = 0;
		int peakEnd = 0;
		int recoveryStart = 0;
		int recoveryEnd = 0;

		double ELA = envC * sqrt(airDensityRef / 2) * pow(4, (envPressureExp - .5));			// Effective Leakage Area
		double NL = 1000 * (ELA / floorArea) * pow((eaveHeight-Hfloor)/2.5, 0.4);				// Normalized Leakage Calculation 62.2-2016 Equation 4.4. 
		double wInfil = (weatherFactor * NL * floorArea) / 1.44;								// Effective annual average infiltration rate, 62.2-2016 Equation 4.5b, L/s.
// 		double defaultInfil = .001 * ((floorArea / 100) * 10) * 3600 / houseVolume;	// Default infiltration credit [ACH] (ASHRAE 62.2, 4.1.3 p.4). This NO LONGER exists in 62.2-2013
// 		double rivecX = (.05 * floorArea + 3.5 * (numBedrooms + 1)) / qRivec;
// 		double rivecY = 0;
//		double Aeq = 0;		// Equivalent air change rate of house to meet 62.2 minimum for RIVEC calculations. Choose one of the following:

// 		switch (AeqCalcs) {
// 			// 1. 62.2-2013 ventilation rate with no infiltration credit [ACH]. Brennan
// 		case 1:
// 			Aeq = .001 * (.15 * floorArea + 3.5 * (numBedrooms + 1)) * (3600 / houseVolume); //Brennan, equals Qtot from 62.2-2013 in air changes per hour
// 			rivecY = 0;
// 			break;
// 			// 2. 62.2 ventilation rate + infiltration credit from weather factors (w x NL)  [ACH]
// 		case 2:
// 			Aeq = .001 * (.05 * floorArea + 3.5 * (numBedrooms + 1)) * (3600 / houseVolume) + wInfil;
// 			rivecY = (wInfil * houseVolume / 3.6)/qRivec;
// 			break;
// 			// 3. 62.2 ventilation rate + 62.2 default infiltration credit (62.2, 4.1.3 p.4)  [ACH]
// 		case 3:
// 			Aeq = .001 * (.05 * floorArea + 3.5 * (numBedrooms + 1) + (floorArea / 100) * 10) * (3600 / houseVolume);
// 			rivecY = ((floorArea / 100) * 10)/qRivec;
// 			break;
// 			// 4. 62.2-2013 ventilation rate with no infiltration credit [ACH], with Existing home deficit of 65 l/s. Brennan
// 		case 4:
// 			Aeq = .001 * (.15 * floorArea + 3.5 * (numBedrooms + 1)+(65/4)) * (3600 / houseVolume); //Brennan, in VentTempControl simulations, this equals an AER 0.319 hr-1 and equals Qtot PLUS an existing home flow deficit of 65 L/s ((65/4).
// 			rivecY = 0;
// 			break;
// 		}

		double expLimit = 5; //Maximum relative exposure limit, 62.2-2016. Old Max way: 1 + 4 * (1 - rivecX) / (1 + rivecY);	// Exposure limit for RIVEC algorithm. This is the Max Sherman way of calculating the maximum limit. 
		//We have moved away from using this, and instead just use a value of 2.5, based on ratios of chronic to acute pollutant exposure limits.

		long int rivecMinutes = 0;					// Counts how many minutes of the year that RIVEC is on
		long int occupiedMinCount = 0;			// Counts the number of minutes in a year that the house is occupied
		
		double relDose = 1;							// Initial value for relative dose used in the RIVEC algorithm
		//double occupiedDose = 1;					// When occupied, this equals relDose, when unoccupied, this equals 0. 
		double totalRelDose = 1;
		//double totalOccupiedDose = 0;				//Cumulative sum for occupiedDose
		//double meanOccupiedDose = 1;				// Mean occupied relative dose over the year
		double meanRelDose = 1;
		double relDoseOld = 1;							// Relative Dose from the previous time step
		
		double relExp = 1;							// Initial value for relative exposure used in the RIVEC algorithm
		double totalRelExp = 1;
		double meanRelExp = 1;
// 		double occupiedExp = 1;						// When occupied, this equals relExp, when unoccupied, this equals 0. 
// 		double totalOccupiedExp = 0;				//Cumulative sum for occupiedExp
// 		double meanOccupiedExp = 1;					// Mean occupied relative exposure over the year
		double relExpOld = 1;						//Initial value for prior minute's relative exposure.
		
		//double turnover = 1 / Aeq;					// Initial value for turnover time (hrs) used in the RIVEC algorithm. Turnover is calculated the same for occupied and unoccupied minutes.					
		double Q_total = 0;							//Total airflow (infiltration + mechanical) (L/s), 62.2-2016
		double Q_wind = 0; 							//Wind airflow (L/s), 62.2-2016 
		double Q_stack = 0; 						//Stack airflow (L/s), 62.2-2016 
		double Q_infiltration = 0; 					//Infiltration airflow (L/s), 62.2-2016 
		
		double outdoorConc = 3.68; //Outdoor average formaldehyde concentration, average is roughly 3 ppb ~ 3.68 ug/m3
		double indoorConc = 20; //Initialized value at typical indoor concentration, ug/m3
		double indoorSource = 23 * floorArea / 3600; //Current default source term (23 ug/m2-hr) is based on the median from low-emitting homes in Hult et al. (2015). Converted here to, ug/s
		//40 ug/m2-hr is a high estimate based on the median from Offermann (2009) winter homes. 
		double DepositionRate = 0.0; // deposition rate in air changes per hour, typical expression in the literature. 
		double qDeposition = DepositionRate * houseVolume / 3600; //deposition rate converted to an equivalent airflow for this particular house, m3/s.
		double penetrationFactor = 1.0; //Penetration factor is the fraction of infiltrating outside mass that gets inside. Default to 1 for no losses. 
		double filterEfficiency = 0.0; //Filtration efficiency in the central HVAC system. Default to 0 for no filter. 

		//For calculating the "real" exposure and dose, based on the actual air change of the house predicted by the mass balance. Standard exposure and dose use the sum of annual average infiltration and current total fan airflow.
// 		double relDoseReal = 1;						// Initial value for relative dose using ACH of house, i.e. the real rel dose not based on ventSum
// 		double relExpReal = 1;						// Initial value for relative exposure using ACH of house i.e. the real rel exposure not based on ventSum
// 		double turnoverReal = 1 / Aeq;			// Initial value for real turnover using ACH of house

		// for calculating the dose and exposure based on the hours of occupancy
		
		
		
// 		double meanRelExp = 0;						// Mean relative exposure over the year, used as a cumulative sum and then ultimately annual average value	
// 		double meanRelDose = 0;						// Mean relative dose over the year, used as a cumulative sum and then ultimately annual average value
// 
// 		double meanRelExpReal = 0;					//Mean "real" relative exposure over the year, used as a cumulative sum and then ultimately annual average value
// 		double meanRelDoseReal = 0;				//Mean "real" relative dose over the year, used as a cumulative sum and then ultimately annual average value
// 
// 		double totalOccupiedExpReal = 0;			//Cumulative sum 
// 		double meanOccupiedExpReal = 0;			//Annual average
// 
// 		double totalOccupiedDoseReal = 0;		// Cumulative sum for "real" calculations. total dose and exp over the occupied time period
// 		double meanOccupiedDoseReal = 0;			// Mean for "real" calculations. total dose and exp over the occupied time period





// 		double turnoverRealOld = 0;
// 		double relDoseRealOld = 0;
		
		double relExpTarget = 1;		//This is a relative expsoure value target used for humidity control simulations. It varies between 0 and 2.5, depending on magnitude of indoor-outdoor humidity ratio difference.

// 		double envC; //Envelope leakage coefficient, m3/s/Pa^n
// 		double envPressureExp; //Envelope pressure exponent
// 		double G; //Wind speed multiplier
// 		double s; //Shelter Factor
// 		double Cs; //Stack coefficient
// 		double Cw; //Wind coefficient
// 		double windSpeed; //Wind speed, corrected for site conditions in main.cpp, m/s
// 		double dryBulb; //Outside temp, K
// 		double ventSum; //Sum of the larger of the mechanical inflows and outflows, ACH
// 		double houseVolume; //House volume, m3
// 		double windSpeedCorrection; //Correction factor used in main.cpp for windspeed
// 	
// 		double Aeq; //Qtot calculated according to 62.2-2016 without infiltration factor, ACH. 
// 		double Q_total; //Total airflow combined infiltration and mechanical, L/s
// 		double relExp_old; //Relative exposure from the prior time-step.
// 		double rivecdt; //RIVCEC timestep, currently defaults to 60/3600, sec. 
// 		double houseVolume; //House volume, m3


		// [START] Fan Schedule Inputs =========================================================================================
		// Read in fan schedule (lists of 1s and 0s, 1 = fan ON, 0 = fan OFF, for every minute of the year)
		// Different schedule file depending on number of bathrooms
		ifstream fanScheduleFile;
		fanScheduleFile.open(fanScheduleFileName); 
		if(!fanScheduleFile) { 
			cout << "Cannot open fan schedule: " << fanScheduleFileName << endl;
			return 1; 		
		}
		// [END] Fan Schedule Inputs =======================================================================================

		int AHminutes;
		int target;
		int dryerFan = 0;			// Dynamic schedule flag for dryer fan (0 or 1)
		int kitchenFan = 0;		// Dynamic schedule flag for kitchen fan (0 or 1)
		int bathOneFan = 0;		// Dynamic schedule flag for first bathroom fan (0 or 1)
		int bathTwoFan = 0;		// Dynamic schedule flag for second bathroom fan (0 or 1)
		int bathThreeFan = 0;	// Dynamic schedule flag for third bathroom fan (0 or 1)
		int weekend;
		int compTime = 0;
		int compTimeCount = 0;
		int rivecOn = 0;		   // 0 (off) or 1 (on) for RIVEC devices
		int mainIterations;
		int ERRCODE = 0;
		int hcFlag = 1;         // Start with HEATING (simulations start in January)

		weatherData weather;		// current minute of weather data
		double mFanCycler;
		double hcap;				// Heating capacity of furnace (or gas burned by furnace)
		double mHRV =0;			// Mass flow of stand-alone HRV unit
		double mHRV_AH = 0;		// Mass flow of HRV unit integrated with the Air Handler
		double mERV_AH = 0;		// Mass flow of ERV unit integrated with the Air Handler
		double fanHeat;
		double ventSumIN;				// Sum of all ventilation flows into house
		double ventSumOUT;			// Sum of all ventilation flows out from house
		//double turnoverOld;			// Turnover from the pervious time step
		double nonRivecVentSumIN;	// Ventilation flows into house not including the RIVEC device flow
		double nonRivecVentSumOUT;	// Ventilation flows out from house not including the RIVEC device flow
		double nonRivecVentSum;		// Equals the larger of nonRivecVentSumIN or nonRivecVentSumOUT
		double qAH;						// Air flowrate of the Air Handler (m^3/s)
		double uaSolAir;
		double uaTOut;
		double rceil;
		double econodt = 0;			// Temperature difference between inside and outside at which the economizer operates
		double supVelAH;			// Air velocity in the supply ducts
		double retVelAH;			// Air velocity in the return ducts
		double qSupReg;				// Airflow rate in the supply registers
		double qRetReg;				// Airflow rate in the return registers
		double qRetLeak;			// Leakage airflow rate in the return ducts
		double qSupLeak;			// Leakage airflow rate in the supply ducts
		double mSupLeak1;
		double mRetLeak1;
		double mSupReg1;
		double mRetReg1;
		double mAH1;
		double mSupReg;
		double mAH;					// Mass flow of air through air handler
		double mRetLeak;			// Mass flow of the air that leaks from the attic into the return ducts
		double mSupLeak;			// Mass flow of the air that leaks from the supply ducts into the attic
		double mRetReg;
		double supVel;				// Supply air velocity
		double retVel;				// Return air velocity
		double AHfanPower;
		double AHfanHeat;
		double mSupAHoff=0;
		double mRetAHoff=0;
		double evapcap = 0;
		double latcap = 0;
		double capacity = 0;
		double capacityh = 0;
		double compressorPower = 0;
		double capacityc = 0;
		double chargecapd = 0;
		double chargeeerd = 0;
		double EER = 0;
		double chargecapw = 0;
		double chargeeerw = 0;
		double Mcoilprevious = 0;
		double mCeilingOld;
		double limit;
		double matticenvout = 0;
		double mCeiling = 0;
		//double mretahaoff;
		double matticenvin = 0;
		double mHouseIN = 0;
		double mCeilingIN = 0; //Ceiling mass flows, not including register flows. Brennan added for ventilation load calculations.
		double mHouseOUT = 0;
		//double RHOATTIC;
		double internalGains;
		double mIN = 0;
		double mOUT = 0;
		double Pint = 0;
		double mFlue = 0;
		double dPflue = 0;
		double Patticint = 0;
		double mAtticIN = 0;
		double mAtticOUT = 0;
		double TSKY = 0;
		double mHouse = 0;
		double qHouse = 0;
		double houseACH = 0;
		double flueACH = 0;
		double ventSum = 0;
		//double ventsumD = 0;
		//double ventSumOUTD = 0;
		//double ventSumIND = 0;
		double Dhda = 0; //Dry-air enthalpy difference, non-ceiling flows [kJ/kg]. Brennan added these elements for calculating ventilation loads.
		double Dhma = 0; //Moist-air enthalpy difference, non-ceiling flows [kJ/kg]. Brennan added these elements for calculating ventilation loads.
		double ceilingDhda = 0; //Dry-air enthalpy difference, ceiling flows [kJ/kg]. Brennan added these elements for calculating ventilation loads.
		double ceilingDhma = 0; //Moist-air enthalpy difference, ceiling flows [kJ/kg]. Brennan added these elements for calculating ventilation loads.
		double DAventLoad = 0; //Dry-air enthalpy load [kJ/min]. Brennan added these elements for calculating ventilation loads.
		double MAventLoad = 0; //Moist-air enthaply load [kJ/min]. Brennan added these elements for calculating ventilation loads.
		double TotalDAventLoad = 0; //Brennan added these elements for calculating ventilation loads [kJ].
		double TotalMAventLoad = 0; //Brennan added these elements for calculating ventilation loads [kJ].
		double HumidityIndex = 0;
		double HumidityIndex_Sum = 0;
		double HumidityIndex_Avg = 0;
		//double coolingLoad = 0;
		//double latLoad = 0;
		//double heatingLoad = 0;
		
		//Mold Index variables
		double moldIndex_South = 0;
		double moldIndex_North = 0;
		double moldIndex_BulkFraming = 0;
		int Time_decl_South = 0;
		int Time_decl_North = 0;
		int Time_decl_Bulk = 0;

		vector<double> averageTemp (0);	// Array to track the 7 day running average
		double HRAttic = 0.008, HRReturn = 0.008, HRHouse = 0.008, HRSupply = 0.008;

		//Variables Brennan added for Temperature Controlled Smart Ventilation.
		double dailyCumulativeTemp = 0;
		double dailyAverageTemp = 0;
		double runningAverageTemp = 0;
		//double AIM2 = 0; //Brennan added
		//double AEQaim2FlowDiff = 0; //Brennan added
		//double FlowDiff = 0; //Brennan added
		//double qFanFlowRatio = 0; //Brennan added
		//double FanQ = fan[0].q; //Brennan's attempt to fix the airflow outside of the if() structures in the fan.oper section. fan[0].q was always the whole house exhaust fan, to be operated continuously or controlled by temperature controls.
		//double FanP = fan[0].power; //Brennan's attempt to fix the fan power outside of the if() structures in the fan.oper section.

		// Call class constructors
		Dehumidifier dh(dhCapacity, dhEnergyFactor, dhSetPoint, dhDeadBand);	// Initialize Dehumidifier 
		Moisture moisture_nodes(atticVolume, retVolume, supVolume, houseVolume, floorArea, planArea, roofPitch, atticMCInit);   // initialize moisture model

		cout << endl;
		cout << "Simulation: " << simNum << endl;
		cout << "Batch File:\t " << batchFileName << endl;
		cout << "Input File:\t " << inputFileName << endl;
		cout << "Output File:\t " << outputFileName << endl;
		cout << "Weather File:\t " << weatherFileName << endl;

		
		// =================================================================
		// ||				 THE SIMULATION LOOPS START HERE:					   ||
		// =================================================================
		//hcFlag = 1;					// Start with HEATING (simulations start in January)
		for(int day = 1; day <= 365; day++) {
			cout << "\rDay = " << day << flush;
			peakFlag = 0;

			dailyAverageTemp = dailyCumulativeTemp / 1440;
			averageTemp.push_back (dailyAverageTemp); //provides the prior day's average temperature...need to do something for day one
			dailyCumulativeTemp = 0;			// Resets daily average outdoor temperature to 0 at beginning of new day
			// For 7 day moving average heating/cooling thermostat decision
			if(day > 7) {
				averageTemp.erase(averageTemp.begin());
				runningAverageTemp = 0;
				for(int i = 0; i < 7 ; i++) {
					runningAverageTemp = runningAverageTemp + averageTemp[i];
				}
				runningAverageTemp = runningAverageTemp/7;
			}

			// Day 1 of simulation is a Sunday then weekend = 1 every Saturday and Sunday, equals 0 rest of the time
			if(day % 7 <= 1)
				weekend = 1;
			else
				weekend = 0;

			// Solar declination and equation of time from ASHRAE HOF 2009 SI ch14
			double dec = 23.45 * sin(360 * (day + 284) / 365.0 * M_PI / 180) * M_PI / 180;
			double gamma = 360.0 * (day - 1) / 365.0 * M_PI / 180;
			double equationOfTime = 2.2918 * (0.0075 + 0.1868 * cos(gamma) - 3.2077 * sin(gamma)
										 - 1.4615 * cos(2 * gamma) - 4.089 * sin(2 * gamma));
			double timeCorrection = equationOfTime / 60 + (longitude - 15 * timeZone) / 15.0;
			// month used for humidity control.
			int month;
			if(day <= 31) {
				month = 1;
			}
			if(day > 31 && day <= 59) {
				month = 2;
			}
			if(day > 60 && day <= 90) {
				month = 3;
			}
			if(day > 91 && day <= 120) {
				month = 4;
			}
			if(day > 121 && day <= 151) {
				month = 5;
			}
			if(day > 152 && day <= 181) {
				month = 6;
			}
			if(day > 182 && day <= 212) {
				month = 7;
			}
			if(day > 213 && day <= 243) {
				month = 8;
			}
			if(day > 244 && day <= 273) {
				month = 9;
			}
			if(day > 274 && day <= 304) {
				month = 10;
			}
			if(day > 305 && day <= 334) {
				month = 11;
			}
			if(day > 335) {
				month = 12;
			}
			// =================================== HOUR LOOP ================================	
			for(int hour = 0; hour < 24; hour++) {
				AHminutes = 0;					// Resetting air handler operation minutes for this hour
				mFanCycler = 0;					// Fan cycler?
				if (hour == peakEnd)
					peakFlag = 1;			// Prevents two peak periods in the same day when there is heating and cooling

				// ============================== MINUTE LOOP ================================	
				for(int minute = 0; minute < 60; minute++) {
					double hourAngle;
					double sinBeta;
					double beta = 0;
					double phi = 0;
					double diffuse = 0;
					double ssolrad = 0;
					double nsolrad = 0;
					double solgain = 0;
					double tsolair;
					double H2, H4, H6;
					double airDensityOUT,airDensityIN,airDensityATTIC,airDensitySUP,airDensityRET;

					target = minute - 39;			// Target is used for fan cycler operation currently set for 20 minutes operation, in the last 20 minutes of the hour.
					if(target < 0)
						target = 0;						// For other time periods, replace the "40" with 60 - operating minutes

					hcap = 0;							// Heating Capacity
					mechVentPower = 0;				// Mechanical Vent Power
					mHRV = 0;							// Mass flow of HRV
					mHRV_AH = 0;						// Mass flow of HRV synced to Air Handler
					fanHeat = 0;
	
					// Resets leakage in case economizer has operated previously (increased leakage for pressure relief)
					if(economizerUsed == 1) {
						envC = Coriginal;
					}
						
					// [START] Filter loading calculations ===============================================================
					if(filterLoadingFlag == 1) {	// Perform filter loading calculations if flag set to 1
						// Air handler airflow rate alterations due to filter loading
						qAH_heat = qAH_heat0 * A_qAH_heat + (k_qAH/100) * qAH_heat0 * massFilter_cumulative;
						qAH_cool = qAH_cool0 * A_qAH_cool + (k_qAH/100) * qAH_cool0 * massFilter_cumulative;
				
						// AH power alterations due to filter loading
						fanPower_heating = fanPower_heating0 * A_wAH_heat + (k_wAH/100) * fanPower_heating0 * massFilter_cumulative;
						fanPower_cooling = fanPower_cooling0 * A_wAH_cool + (k_wAH/100) * fanPower_cooling0 * massFilter_cumulative;
				
						// Return duct leakage alterations due to filter loading
						retLF = retLF * A_DL + (k_DL/100) * retLF * massFilter_cumulative;
				
						// Change filter when airflow rate drops to half lowest speed (reset AH airflow, power and duct leakage)
						if((qAH_heat <= (0.5 * qAH_low)) || (qAH_cool <= (0.5 * qAH_low))) {
								qAH_heat = qAH_heat0 * A_qAH_heat;
								qAH_cool = qAH_cool0 * A_qAH_cool;
								fanPower_heating = fanPower_heating0 * A_wAH_heat;
								fanPower_cooling = fanPower_cooling0 * A_wAH_cool;
								retLF = retLF * A_DL;
								massFilter_cumulative = 0;
								filterChanges++;
							}
						// Cooling capacity air flow correction term (in CFM)
						qAH_cfm = qAH_cool / .0004719;
						qAHcorr = 1.62 - .62 * qAH_cfm / (400 * capacityraw) + .647 * log(qAH_cfm / (400 * capacityraw));
					}
					// [END] Filter loading calculations =================================================================

					ventSumIN = 0;						// Setting sum of supply mechanical ventilation to zero
					ventSumOUT = 0;						// Setting sum of exhaust mechanical ventilation to zero

					// RIVEC dose and exposure calculations
					//turnoverOld = turnover;
					

					// Dose and Exposure calculations based on ACH of house rather than just mech vent
// 					turnoverRealOld = turnoverReal;
// 					relDoseRealOld = relDoseReal;

					nonRivecVentSumIN = 0;				// Setting sum of non-RIVEC supply mechanical ventilation to zero
					nonRivecVentSumOUT = 0;				// Setting sum of non-RIVEC exhaust mechanical ventilation to zero

					// Read in or interpolate weather data
					if(weatherFileType == 0) {  // minute
						weather = readOneMinuteWeather(weatherFile);
					}
					else {
						weather = interpWeather(begin, end, minute);
					}
					weather.windSpeed *= windSpeedCorrection;		// Correct met wind speed to speed at building eaves height [m/s]
					if(weather.windSpeed < 1)					// Minimum wind velocity allowed is 1 m/s to account for non-zero start up velocity of anenometers
						weather.windSpeed = 1;					// Wind speed is never zero
					dailyCumulativeTemp = dailyCumulativeTemp + weather.dryBulb - C_TO_K;
					tsolair = weather.dryBulb;

					for(int k=0; k < 4; k++)			// Wind direction as a compass direction?
						Sw[k] = pow(Swinit[k][weather.windDirection],2);		// store square of Sw to pass to attic_leak and house_leak

					// Fan Schedule Inputs
					// Assumes operation of dryer and kitchen fans, then 1 - 3 bathroom fans
					fanScheduleFile >> dryerFan >> kitchenFan >> bathOneFan >> bathTwoFan >> bathThreeFan;

					if(minuteYear == 1) {				// Setting initial humidity conditions
						for(int i = 0; i < 5; i++) {
							HR[i] = weather.humidityRatio;
						}
					}

					// Calculate air densities
					airDensityOUT = airDensityRef * airTempRef / weather.dryBulb;		// Outside Air Density
					airDensityIN = airDensityRef * airTempRef / tempHouse;		// Inside Air Density
					airDensityATTIC = airDensityRef * airTempRef / tempAttic;	// Attic Air Density
					airDensitySUP = airDensityRef * airTempRef / tempSupply;		// Supply Duct Air Density
					airDensityRET = airDensityRef * airTempRef / tempReturn;		// Return Duct Air Density

					// Solar calculations
					hourAngle = 15 * (hour + minute / 60.0 + timeCorrection - 12) * M_PI / 180;
					sinBeta = cos(latitude) * cos(dec) * cos(hourAngle) + sin(latitude) * sin(dec);
					if(sinBeta > 0) {			// sun is up
						beta = asin(sinBeta);
						phi = copysign(acos((sinBeta * sin(latitude) - sin(dec)) / (cos(beta) * cos(latitude))),hourAngle);
						double hdirect = sinBeta * weather.directNormal;
						diffuse = max(0.0,weather.globalHorizontal - hdirect);

						// Roof solar gain
						double sigma = M_PI * roofPitch / 180;				// roof tilt in radians
						ssolrad = surfaceInsolation(weather.directNormal, diffuse, beta, sigma, phi - 0);
						nsolrad = surfaceInsolation(weather.directNormal, diffuse, beta, sigma, phi - M_PI);
					
						// Wall and window solar gain
						double totalSolar = 0;
						double incsolar[4];
						double psi[4] = {0, M_PI/2, M_PI, -M_PI/2};
						for(int i=0; i < 4; i++) {
							incsolar[i] = surfaceInsolation(weather.directNormal, diffuse, beta, M_PI/2, phi - psi[i]);
							totalSolar += incsolar[i];
						}
						solgain = winShadingCoef * (windowS * incsolar[0] + windowWE / 2 * incsolar[1] + windowN * incsolar[2] + windowWE / 2 * incsolar[3]);
						tsolair = totalSolar / 4 * .03 + weather.dryBulb;		// the .03 is from 1993 AHSRAE Fund. SI 26.5
					}
					// 7 day outdoor temperature running average. If <= 60F we're heating. If > 60F we're cooling
					if(runningAverageTemp <= (60 - 32) * 5.0 / 9.0)
						hcFlag = 1;					// Heating
					else
						hcFlag = 2;					// Cooling
			
					// ====================== HEATING THERMOSTAT CALCULATIONS ============================
					if(hcFlag == 1) {
						qAH = qAH_heat;            // Heating Air Flow Rate [m3/s]
						uaSolAir = uaWall;
						uaTOut = uaFloor + uaWindow;
						rceil = ceilRval_heat;

						if(tempOld[15] > (heatThermostat[hour] + .5))
							AHflag = 0;					// Heat off if building air temp above setpoint

						if(AHflag == 0 || AHflag == 100) {
							if(tempOld[15] <= (heatThermostat[hour] - .5))
								set = 1;					// Air handler off, and tin below setpoint - 0.5 degrees
							else
								set = 0;					// Air handler off, but tin above setpoint - 0.5 degrees
						}

						if(tempOld[15] < (heatThermostat[hour] - .5))
							AHflag = 1;					// turn on air handler/heat

						if(tempOld[15] >= (heatThermostat[hour] - .5) && tempOld[15] <= (heatThermostat[hour] + .5)) {
							if(set == 1)
								AHflag = 1;
							else
								AHflag = 0;
						}

						if(AHflag == 0 && AHflag != AHflagPrev)
							endrunon = minuteYear + 1;		// Adding 1 minute runon during which heat from beginning of cycle is put into air stream

						if(AHflag == 1)
							hcap = hcapacity / AFUE;	// hcap is gas consumed so needs to divide output (hcapacity) by AFUE

						// ====================== COOLING THERMOSTAT CALCULATIONS ========================
					} else {
						hcap = 0;
						qAH = qAH_cool;						// Cooling Air Flow Rate
						uaSolAir = uaWall;
						uaTOut = uaWindow;
						rceil = ceilRval_cool;
						endrunon = 0;

						if(tempOld[15] < (coolThermostat[hour] - .5))
							AHflag = 0;

						if(AHflag == 0 || AHflag == 100) {
							if(tempOld[15] >= (coolThermostat[hour] + .5))
								set = 1;				// 1 = AH OFF and house temp below thermostat setpoint - 0.5
							else
								set = 0;				// 0 = AH OFF and house temp above thermostat setpoint - 0.5
						}

						if(tempOld[15] >= (coolThermostat[hour] + .5)) {
							AHflag = 2;
						}

						if(tempOld[15] < (coolThermostat[hour] + .5) && tempOld[15] >= (coolThermostat[hour] - .5)) {
							if(set == 1)
								AHflag = 2;
							else
								AHflag = 0;
						}

						if(AHflag != AHflagPrev && AHflag == 2) {				// First minute of operation
							compTime = 1;
							compTimeCount = 1;
						} else if(AHflag == 2 && compTimeCount == 1) {			// Second minute of operation
							compTime = 2;
							compTimeCount = 0;
						} else if(AHflag == 2 && compTimeCount == 0) {
							compTime = 0;
						} else if(AHflag != AHflagPrev && AHflag == 0) {		// End of cooling cycle
							compTime = 0;
						}

						// [START] ====================== ECONOMIZER RATIONALE ===============================
						econodt = tempHouse - weather.dryBulb;			// 3.333K = 6F

						// No cooling, 6F dt (internal/external) and tempHouse > 21C for economizer operation. 294.15 = 21C
							//if(AHflag == 0 && econodt >= 3.333 && tempHouse > 294.15 && economizerUsed == 1) {
						if(AHflag == 0 && econodt >= 3.33 && tempHouse > 294.15 && economizerUsed == 1 && hcFlag == 2) {
								econoFlag = 1;
								envC = Ceconomizer;
							} else {
								econoFlag = 0;
							}
					}	// [END] ========================== END ECONOMIZER ====================================



					// [START] RIVEC Decision ==================================================================================================================

					// Choose RIVEC time periods depending on heating or cooling
					switch (hcFlag) {

					case 1:	// HEATING TIMES
						baseStart		= 12;
						baseEnd			= 4;
						peakStart		= 4;
						peakEnd			= 8;
						recoveryStart	= 8;
						recoveryEnd		= 12;
						break;

					case 2: // COOLING TIMES
						baseStart		= 22;
						baseEnd			= 14;
						peakStart		= 14;
						peakEnd			= 18;
						recoveryStart	= 18;
						recoveryEnd		= 22;
						break;
					}


					if(minute == 0 || minute == 10 || minute == 20 || minute == 30 || minute == 40 || minute == 50) {
// 						for(int i=0; i < numFans; i++) {
// 							// [START] ---------------------FAN 50---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v6
// 							// v6 of the algorithm only uses the peakStart and peakEnd variables, no more base or recovery periods.(
// 							if(fan[i].oper == 50 || fan[i].oper == 13 || fan[i].oper == 17) { //traditional (50), cfis (13) or erv+ahu (17) fans for rivec control
// 								// rivecOn = 1 or 0: 1 = whole-house fan ON, 0 = whole-house fan OFF
// 								rivecOn = 0;
// 								//// RIVEC fan operation after algorithm decision
// 								//if(rivecOn) {	  												// RIVEC has turned ON this fan
// 								//	fan[i].on = 1;
// 								//	mechVentPower = mechVentPower + fan[i].power;				// vent fan power
// 								//	if(fan[i].q > 0) { 											// supply fan - its heat needs to be added to the internal gains of the house
// 								//		fanHeat = fan[i].power * .84;							// 16% efficiency for this fan
// 								//		ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
// 								//	} else { 													// exhaust fan
// 								//		ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
// 								//	}
// 								//	rivecMinutes++;
// 								//}
// 								//else {
// 								//	fan[i].on = 0;
// 								//}
// 							}
// 						//}
						if(OccContType == 2){ //Aux fan control only.
							if(relExp >= 1.0 || relDose > 1.0){
								rivecOn = 1;
							} else {
								rivecOn = 0;
							}
						}
					
						else if(OccContType == 3 || OccContType == 4){ //Occupancy control only (3) or Aux Fans + Occupancy control (4). Aux fans are accounted for with the AuxFanIndex variable. 
							if(occupied[weekend][hour] == 1){ //occupied
								if(relExp >= 1.0 || relDose > 1.0){
									rivecOn = 1;
								} else {
									rivecOn = 0;
								}
								
							} else{ //unoccupied
								if(relExp > expLimit){ //expLimit defaults to 5, based on ASHRAE 62.2-2016 Addendum C. 
									rivecOn = 1;
								} else{
									rivecOn = 0;
								}
							}
						}				
					}
						// [END] ========================== END RIVEC Decision ====================================

	// 					========================== Start Humidity Control Logic ================================
// 						Cooling system tie-in.
// 						if(HumContType == 1){			   				
// 							if(hcFlag == 2){ //test if we're in cooling season
// 								if(AHflag == 2){
// 									rivecOn = 1;
// 								} else
// 									if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;					
// 									}
// 							} else //if NOT in cooling season
// 								if(relExp >= 0.95 || relDose > 1.0){ //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;					
// 									}
// 							}
// 
// 						Fixed control. 	   				
// 						Indoor and Outdoor sensor based control. 
// 						if(HumContType == 2){
// 							if(RHhouse >= 55){ //Engage increased or decreased ventilation only if house RH is >60 (or 55% maybe?). 
// 								if(weather.humidityRatio > HR[3]){ //do not want to vent. Add some "by what amount" deadband value. 
// 									if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relExp >= 0.50 || relDose > 1.0){ //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(relExp >= 0.95 || relDose > 1.0){ 
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Fixed control + cooling system tie-in.		   				
// 						Indoor and Outdoor sensor based control. 
// 						if(HumContType == 3){							
// 							if(RHhouse >= 55){ //Engage increased or decreased ventilation only if house RH is >60 (or 55% maybe?). 
// 								if(weather.humidityRatio > HR[3]){ //do not want to vent. Add some "by what amount" deadband value. 
// 									if(AHflag == 2){
// 										rivecOn = 1;
// 									} else
// 										if(relExp >= 2.5 || relDose > 1.0){ //have to with high exp
// 											rivecOn = 1;
// 										} else { //otherwise off
// 											rivecOn = 0;
// 										}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relExp >= 0.50 || relDose > 1.0){ //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(relExp >= 0.95 || relDose > 1.0){ 
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Proportional control.  				   				
// 						Indoor and Outdoor sensor based control. 
// 						if(HumContType == 4) {
// 							if(RHhouse >= 55) {
// 								if(weather.humidityRatio > HR[3]){ //More humid outside than inside, want to under-vent.  
// 									relExpTarget = 1 + (2.5-1) * abs((HR[3]-weather.humidityRatio) / (wDiffMaxNeg)); //wDiffMax has to be an avergaed value, because in a real-world controller you would not know this. 
// 									if(relExpTarget > 2.5){
// 										relExpTarget = 2.5;
// 									}
// 									if(relExp >= relExpTarget || relDose > 1){ //relDose may be over a 1-week or 2-week time span...
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { // More humid inside than outside, want to over-vent
// 									relExpTarget = 1 - abs((HR[3]- weather.humidityRatio) / (wDiffMaxPos));
// 									if(relExpTarget < 0){
// 										relExpTarget = 0;
// 									}
// 									if(relExp >= relExpTarget || relDose > 1){ //
// 										rivecOn = 1;
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(relExp >= 0.95 || relDose > 1){ 
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Proportional control + cooling system tie-in. 				   				
// 						Indoor and Outdoor sensor based control.
// 						if(HumContType == 5) { 
// 							if(RHhouse >= 55) {
// 								if(weather.humidityRatio > HR[3]) { //More humid outside than inside, want to under-vent.  
// 									relExpTarget = 1 + (2.5-1) * abs((HR[3]-weather.humidityRatio) / (wDiffMaxNeg)); //wDiffMax has to be an avergaed value, because in a real-world controller you would not know this. 
// 									if(relExpTarget > 2.5) {
// 										relExpTarget = 2.5;
// 									}
// 									if(AHflag == 2) {
// 										rivecOn = 1;
// 									} else if(relExp >= relExpTarget || relDose > 1) { //relDose may be over a 1-week or 2-week time span...
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { // More humid inside than outside, want to over-vent
// 									relExpTarget = 1 - abs((HR[3]- weather.humidityRatio) / (wDiffMaxPos));
// 									if(relExpTarget < 0) {
// 										relExpTarget = 0;
// 									}
// 									if(relExp >= relExpTarget || relDose > 1) {
// 										rivecOn = 1;
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(relExp >= 0.95 || relDose > 1) { 
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Monthly Seasonal Control
// 						Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
// 						if(HumContType == 6) {
// 							if(month <= FirstCut || month >= SecondCut) { //High ventilation months with net-humidity transport from inside to outside.
// 								if(relDose > doseTarget) {
// 									rivecOn = 1;
// 								} else {						
// 									rivecOn = 0;
// 								} 
// 							} else { //Low ventilation months with net-humidity transport from outside to inside.
// 								if(relExp >= 2.5 || relDose > HiDose) { //Brennan changed from fixed 1.5 to variable HiDose value.
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Fixed control + cooling system tie-in + Monthly Seasonal Control.		   				
// 						Indoor and Outdoor sensor based control.
// 						if(HumContType == 7) { 
// 							if(weather.humidityRatio > HR[3]) { //do not want to vent. 
// 								if(AHflag == 2) {
// 									rivecOn = 1;
// 								} else if(relExp >= 2.5 || relDose > HiDose) { 
// 									rivecOn = 1;
// 								} else { //otherwise off
// 									rivecOn = 0;
// 								}
// 							} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 								if(relDose > doseTarget) { 
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Fixed control + cooling system tie-in + Monthly Seasonal Control.		   				
// 						Indoor and Outdoor sensor based control.
// 						if(HumContType == 7) { 
// 							if(weather.humidityRatio > HR[3]) { //do not want to vent. 
// 								if(RHhouse >= 55){
// 								if(AHflag == 2) {
// 									rivecOn = 1;
// 								} else if(relExp >= 2.5 || relDose > HiDose) { 
// 									rivecOn = 1;
// 								} else { //otherwise off
// 									rivecOn = 0;
// 								}
// 								}
// 								else if(relExp >= 0.95 || relDose > 1){ 
// 									rivecOn = 1;
// 								} 
// 								else{
// 									rivecOn = 0;
// 								}		
// 							}		
// 							else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 								if(relDose > doseTarget) { 
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Monthly Seasonal controller + time of day
// 						Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average
// 						targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
// 						if(HumContType == 8) {
// 							if(month <= FirstCut || month >= SecondCut) { //High ventilation months with net-humidity transport from inside to outside.
// 								if(hour >= 3 && hour <= 7) { //Was 3 and 7. Maybe change this to the warmest hours of the day.
// 									rivecOn = 1; //Relatively dry time of day, vent more.
// 								} else {
// 									if(relDose > doseTarget) {
// 										rivecOn = 1;
// 									} else {						
// 										rivecOn = 0;
// 									} 
// 								}
// 							} else { //Low ventilation months with net-humidity transport from outside to inside.
// 								if(hour >= 14 && hour <= 18) { //Brennan changed from 12 and 6, to 14 to 18. Always vent during peak cooling period.
// 									rivecOn = 1;
// 								} else if(hour >= 7 && hour <= 11){
// 									rivecOn = 0;
// 								} else {
// 									if(relExp >= 2.5 || relDose > HiDose) {
// 										rivecOn = 1;
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							}
// 						}
// 
// 						Fixed control + Monthly Seasonal Control.		   				
// 						Indoor and Outdoor sensor based control.
// 						if(HumContType == 9) { 
// 							if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 								if(relExp >= 2.5 || relDose > HiDose) { //have to with high exp 0.61
// 									rivecOn = 1;
// 								} else { //otherwise off
// 									rivecOn = 0;
// 								}
// 							} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 								if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
// 									rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 					
// 						Fixed control + Monthly Seasonal Control.		   				
// 						Indoor and Outdoor sensor based control.
// 						if(HumContType == 9) { 
// 							if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 								if(RHhouse >= 55) {
// 								if(relExp >= 2.5 || relDose > HiDose) { //have to with high exp 0.61
// 									rivecOn = 1;
// 								} 
// 								else { //otherwise off
// 									rivecOn = 0;
// 								}
// 								}	
// 								else if(relExp >= 0.95 || relDose > 1){ 
// 									rivecOn = 1;
// 								} 
// 								else{
// 									rivecOn = 0;
// 								}	
// 							}			
// 							else if(relDose > doseTarget){ //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 								rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 							} 
// 							else {
// 								rivecOn = 0;
// 							}
// 						}						
// 					
// 
// 						The real opportunities for control based on outside are when the outside value is changing rapidly. Sharp increases, decrease ventilation. Sharp decreases, increase ventilation. 
// 						Need to undervent at above the 75th percentile and overvent below the 25th percentile based on a per month basis. 
// 						Outdoor-only sensor based control
// 						if(HumContType == 10) {
// 							if(weather.humidityRatio > 0.012) { //If humid outside, reduce ventilation. 
// 								if(relExp >= 2.5 || relDose > 1) {
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							} else { //If dry outside, increase vnetilation.
// 								if(relExp >= 0.95 || relDose > 1) { //but we do if exp is high
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0; //otherwise don't vent under high humidity condition.
// 								}
// 							}
// 						}
// 
// 						Monthly Advanced Seasonal Control
// 						Monthly timer-based control, based on mean HRdiff by month. 
// 						doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
// 						Setting the appropriate Dose Target based on the month						
// 						if(HumContType == 11) {
// 							double doseTargetTmp;
// 							if(month == HiMonths[0] || month == HiMonths[1] || month == HiMonths[2]){
// 								doseTargetTmp = HiMonthDose;
// 							} else if(month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
// 								doseTargetTmp = LowMonthDose;
// 							} else if(month > FirstCut && month < SecondCut) {
// 								doseTargetTmp = HiDose; //Brennan changed from fixed 1.5 to HiDose.
// 							} else {
// 								doseTargetTmp = doseTarget;
// 							}
// 							if(doseTargetTmp >= HiDose) {
// 								if(relExp >= 2.5 || relDose > doseTargetTmp) {
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 							else {
// 								if(relDose > doseTargetTmp) {
// 									rivecOn = 1;
// 								} else {						
// 									rivecOn = 0;
// 								} 
// 							}
// 						}
// 
// Consider combining 8 and 12 - Monthly + Time of day + Cooling tie-in.
// 						Monthly Seasonal Control + Cooling system tie-in.		   				
// 						Indoor and Outdoor sensor based control. 
// 						if(HumContType == 12) {
// 							double doseTargetTmp;
// 							if(month <= FirstCut || month >= SecondCut) {
// 								doseTargetTmp = doseTarget;
// 							} else {
// 								doseTargetTmp = HiDose;
// 							}
// 							if(AHflag == 2) {
// 									rivecOn = 1;
// 							} else if(doseTargetTmp >= HiDose) {
// 								if(relExp >= 2.5 || relDose > doseTargetTmp) {
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							} else {
// 								if(relDose > doseTargetTmp) {
// 									rivecOn = 1;
// 								} else {						
// 									rivecOn = 0;
// 								} 
// 							}
// 						}
// 
// 						Outdoor-only sensor based control, with variable dose targets
// 						if(HumContType == 13) { 
// 							if(weather.humidityRatio >= wCutoff) { //If humid outside, reduce ventilation. 
// 								if(relExp >= 2.5 || relDose > 1.45) {
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							} else { //If dry outside, increase ventilation.
// 								if(relDose > 0.45) { //but we do if exp is high
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0; //otherwise don't vent under high humidity condition.
// 								}
// 							}
// 						}
// 
// 						Outdoor-only sensor based control, control based on 25th and 75th percentile monthly values for each month and climate zone.
// 						if(HumContType == 14) { 
// 							if(weather.humidityRatio >= W75[month-1]) { //If humid outside, reduce ventilation. 
// 								if(relExp >= 2.5 || relDose > 1.5) {
// 									rivecOn = 1; 
// 								} else {
// 									rivecOn = 0;
// 								}
// 							} else if (weather.humidityRatio <= W25[month-1]) { //If dry outside, increase vnetilation.
// 								if(relDose > 0.5) { //but we do if exp is high
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0; //otherwise don't vent under high humidity condition.
// 								}
// 							} else {
// 								if(relExp >= 0.95 || relDose > 1.0){ //Need to reduce this target to make equivalence work out...
// 									rivecOn = 1;
// 								} else {
// 									rivecOn = 0;
// 								}
// 							}
// 						}
// 
// 						Monthly Advanced Seasonal Control + Fixed Control + Cooling System Tie_in
// 						Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
// 						Setting the appropriate Dose Target based on the month
// 						if(HumContType == 15) {
// 							if(month == HiMonths[0] || month == HiMonths[1]  || month == HiMonths[2] ||
// 								month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
// 								doseTargetTmp = HiMonthDose;
// 								if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 									if(AHflag == 2) {
// 										rivecOn = 1;
// 									} else if(relExp >= 2.5 || relDose > LowMonthDose) { //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relDose > HiMonthDose) { //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 									if(AHflag == 2) {
// 										rivecOn = 1;
// 									} else if(relExp >= 2.5 || relDose > 1.5) { //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							}
// 						}
// 
// 						Monthly Advanced Seasonal Control + Fixed Control
// 						Monthly timer-based control, based on mean HRdiff by month. doseTargets are based on weighted average targeting dose = 1.5 during low-ventilation months, targeting annual dose of 0.98. 
// 						Setting the appropriate Dose Target based on the month
// 						if(HumContType == 16) {
// 							if(month == HiMonths[0] || month == HiMonths[1]  || month == HiMonths[2] || 
// 								month == LowMonths[0] || month == LowMonths[1]  || month == LowMonths[2]) {
// 								doseTargetTmp = HiMonthDose;
// 								if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 									if(relExp >= 2.5 || relDose > LowMonthDose) { //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relDose > HiMonthDose) { //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							} else {
// 								if(weather.humidityRatio > HR[3]) { //do not want to vent. Add some "by what amount" deadband value. 
// 									if(relExp >= 2.5 || relDose > 1.5) { //have to with high exp
// 										rivecOn = 1;
// 									} else { //otherwise off
// 										rivecOn = 0;
// 									}
// 								} else { //want to vent due to high indoor humidity, so maybe we just let it run, without relExp control?
// 									if(relDose > doseTarget) { //control to exp = 0.5. Need to change this value based on weighted avg results.
// 										rivecOn = 1; //OR we can change the does calculation based on expected periods of contol function (i.e., 1-week,1-month, etc.)
// 									} else {
// 										rivecOn = 0;
// 									}
// 								}
// 							}
// 						}
// 
// 						Rivec control of number 50 fan type, algorithm v6					     
// 						if(HumContType == 0) {
// 							if(occupied[weekend][hour]) {				            	// Base occupied
// 								if(relExp >= 0.95 || relDose >= 1.0)
// 									rivecOn = 1;
// 							} else {						                	// Base unoccupied
// 								if(relExp >= expLimit)
// 									rivecOn = 1;
// 							}
// 							if(hour >= peakStart && hour < peakEnd && peakFlag == 0) {		// PEAK Time Period
// 								rivecOn = 0;												// Always off
// 								if(relExp >= expLimit)
// 									rivecOn = 1;
// 							}
// 						}
// 						================== End Humidity Control Logic ================================
// 					}
				
					AHflagPrev = AHflag;
					if(AHflag != 0)
						AHminutes++;			// counting number of air handler operation minutes in the hour

					// the following air flows depend on if we are heating or cooling
					supVelAH = qAH / (pow(supDiameter,2) * M_PI / 4);
					retVelAH = qAH / (pow(retDiameter,2) * M_PI / 4);
					qSupReg = -qAH * supLF + qAH;
					qRetReg = qAH * retLF - qAH;
					qRetLeak = -qAH * retLF;
					qSupLeak = qAH * supLF;
					mSupLeak1 = qSupLeak * airDensityIN;
					mRetLeak1 = qRetLeak * airDensityIN;
					mSupReg1 = qSupReg * airDensityIN;
					mRetReg1 = qRetReg * airDensityIN;
					mAH1 = qAH * airDensityIN;			// Initial mass flow through air handler
					mFanCycler = 0;
					mechVentPower = 0;					// Total vent fan power consumption

					if(AHflag == 0) {			// AH OFF
						// the 0 at the end of these terms means that the AH is off
						if(minuteYear >= endrunon) {  // endrunon is the end of the heating cycle + 1 minutes
							mSupReg = 0;
							mAH = 0;
							mRetLeak = 0;
							mSupLeak = 0;
							mRetReg = 0;
							supVel = abs(mSupAHoff) / airDensitySUP / (pow(supDiameter,2) * M_PI / 4);
							retVel = abs(mRetAHoff) / airDensityRET / (pow(retDiameter,2) * M_PI / 4);
							AHfanPower = 0;		// Fan power consumption [W]
							AHfanHeat = 0;		// Fan heat into air stream [W]
							hcap = 0;			// Gas burned by furnace NOT heat output (that is hcapacity)
						} else {
							AHflag = 102;		// Note: need to have capacity changed to one quarter of burner capacity during cool down
							mAH = mAH1;
							mRetLeak = mRetLeak1;
							mSupLeak = mSupLeak1;
							mRetReg = mRetReg1;
							mSupReg = mSupReg1;
							supVel = supVelAH;
							retVel = retVelAH;
							AHfanHeat = fanPower_heating * 0.85;	//0.85		// Heating fan power multiplied by an efficiency
							AHfanPower = fanPower_heating;
							hcap = hcapacity / AFUE;

							for(int i = 0; i < numFans; i++) {     // For outside air into return cycle operation
								if(fan[i].oper == 6 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 10) {			// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 14) {			// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 18) {			// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 15 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 13) {
									if((rivecFlag == 0 && AHminutes <= 20) || (rivecFlag == 1 && rivecOn == 1)){ //This allows RIVEC to operate this fan until it reaches its Exposure setpoint?
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
									}
								}
								if(fan[i].oper == 22 && AHminutes <= 20 && econoFlag == 0) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
							}
						}
					} else {
						if(AHflag == 2) {	//	cooling 
							hcap = 0;
							mAH = mAH1;
							mRetLeak = mRetLeak1;
							mSupLeak = mSupLeak1;
							mRetReg = mRetReg1;
							mSupReg = mSupReg1;
							supVel = supVelAH;
							retVel = retVelAH;
							AHfanHeat = fanPower_cooling * 0.85;	//0.85		// Cooling fan power multiplied by an efficiency (15% efficient fan)
							AHfanPower = fanPower_cooling;
							for(int i = 0; i < numFans; i++) {			// for outside air into return cycle operation
								if(fan[i].oper == 6 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN * qAH_cool / qAH_heat;	// correcting this return intake flow so it remains a constant fraction of air handler flow
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 10) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 14) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 18) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 15 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 13) {
									if((rivecFlag == 0 && AHminutes <= 20) || (rivecFlag == 1 && rivecOn == 1)){ //This allows RIVEC to operate this fan until it reaches its Exposure setpoint?
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;								
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
									}
								}
								if(fan[i].oper == 22 && AHminutes <= 20 && econoFlag == 0) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
							}
						} else {
							mAH = mAH1;
							mRetLeak = mRetLeak1;
							mSupLeak = mSupLeak1;
							mRetReg = mRetReg1;
							mSupReg = mSupReg1;
							supVel = supVelAH;
							retVel = retVelAH;
							AHfanHeat = fanPower_heating * 0.85;	//0.85
							AHfanPower = fanPower_heating;
							hcap = hcapacity / AFUE;

							for(int i = 0; i < numFans; i++) {			// for outside air into return cycle operation
								if(fan[i].oper == 6 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 10) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 14) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 18) {				// vent always open during any air handler operation
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 15 && AHminutes <= 20) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
								if(fan[i].oper == 13) { 
									if((rivecFlag == 0 && AHminutes <= 20) || (rivecFlag == 1 && rivecOn == 1)){ //This allow RIVEC to operate this fan until it reaches its Exposure setpoint?
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;								
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
									}
								}
								if(fan[i].oper == 22 && AHminutes <= 20 && econoFlag == 0) {
									mFanCycler = fan[i].q * airDensityIN;
									mRetReg = mRetReg1 - mFanCycler;
								}
							}
						}
					}

					for(int i=0; i < numFans; i++) {
						// ***********************  20 minute minimum central fan integrated supply (CFIS) system. 20 minute minimum operation per hour. 
						if(fan[i].oper == 13) {			// fan cycler operation without exhaust fan
							if((rivecFlag == 0 && AHminutes < target && AHflag != 1 && AHflag != 2 && AHflag != 102) || (rivecFlag == 1 && rivecOn == 1 && AHflag != 1 && AHflag != 2 && AHflag != 102)){		// conditions for turning on air handler only for fan cycler operation
								AHflag = 100;
								AHminutes = AHminutes + 1;
								qAH = qAH_cool;
								supVelAH = qAH / (pow(supDiameter,2) * M_PI / 4);
								retVelAH = qAH / (pow(retDiameter,2) * M_PI / 4);
								qSupReg = -qAH * supLF + qAH;
								qRetReg = qAH * retLF - qAH;
								qRetLeak = -qAH * retLF;
								qSupLeak = qAH * supLF;
								mSupLeak = qSupLeak * airDensityIN;
								mRetLeak = qRetLeak * airDensityIN;
								mFanCycler = fan[i].q * airDensityIN;
								mSupReg = qSupReg * airDensityIN;
								mRetReg = qRetReg * airDensityIN - mFanCycler;							
								mAH = qAH * airDensityIN;
								mechVentPower = mechVentPower + fanPower_cooling;			// note ah operates at cooling speed
								AHfanHeat = fanPower_cooling * .85;							// for heat
								ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
								//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
							}
						}
					//}

						if(fan[i].oper == 22) {						// Fan cycler operation only when economizer is off
							if(econoFlag == 0) {					// Only open CFI if economizer is off
								if(AHminutes < target && AHflag != 1 && AHflag != 2 && AHflag != 102) {		// conditions for turning on air handler only for fan cycler operation
									AHflag = 100;
									AHminutes = AHminutes + 1;
									qAH = qAH_cool;
									supVelAH = qAH / (pow(supDiameter,2) * M_PI / 4);
									retVelAH = qAH / (pow(retDiameter,2) * M_PI / 4);
									qSupReg = -qAH * supLF + qAH;
									qRetReg = qAH * retLF - qAH;
									qRetLeak = -qAH * retLF;
									qSupLeak = qAH * supLF;
									mSupLeak = qSupLeak * airDensityIN;
									mRetLeak = qRetLeak * airDensityIN;
									mFanCycler = fan[i].q * airDensityIN;
									mSupReg = qSupReg * airDensityIN;
									mRetReg = qRetReg * airDensityIN - mFanCycler;
									mAH = qAH * airDensityIN;
									mechVentPower = mechVentPower + fanPower_cooling;    // note ah operates at cooling speed
									AHfanHeat = fanPower_cooling * .85;  // for heat
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
								}
							}
						}
					}


					// [START] HRV's and ERV's===================================================================================================================================
					// (Can be RIVEC controlled so after RIVEC algorithm decision)
					for(int i=0; i < numFans; i++) {

						//Brennan. Stand-alone ERV MUST be entered two times in the input files, with half the fan power attributed to each. Same airflow, but one positive and one negative.

						if(fan[i].oper == 5) {		// Standalone HRV not synched to air handler

							if((rivecFlag == 0 && minute >= 30) || (rivecFlag == 1 && rivecOn == 1)) {	// on for last 30 minutes of every hour or RIVEC controlled
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								if(fan[i].q > 0) {
									mHRV = fan[i].q * airDensityOUT;;
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;	// use if HRV not to be included in RIVEC whole house ventilation calculations
								} else {
									ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
									//nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								}
							} else {
								fan[i].on = 0;
							}
						}

						if(fan[i].oper == 16) {		// HRV + Air Handler

							if((rivecFlag == 0 && minute >= 30) || (rivecFlag == 1 && rivecOn == 1)) {	// on for last 30 minutes of every hour or RIVEC controlled
								fan[i].on = 1;
								ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
								ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;

								mechVentPower = mechVentPower + fan[i].power;
								if(AHflag != 1 && AHflag != 2 && AHflag != 102) {			// ah is off for heat/cool so turn it on for venting
									AHflag = 100;
									AHminutes = AHminutes + 1;
									qAH = qAH_cool;
									supVelAH = qAH / (pow(supDiameter,2) * M_PI / 4);
									retVelAH = qAH / (pow(retDiameter,2) * M_PI / 4);
									qSupReg = -qAH * supLF + qAH;
									qRetReg = qAH * retLF - qAH;
									qRetLeak = -qAH * retLF;
									qSupLeak = qAH * supLF;
									mSupLeak = qSupLeak * airDensityIN;
									mAH = qAH * airDensityIN;
									mHRV_AH = abs(fan[i].q * airDensityIN) * -1.0;			// mHRV_AH needs to be a negative number - so does fan(i).q to get right mRetReg and then fan(i).q is also the exhaust from house flow
									mSupReg = qSupReg * airDensityIN;
									mRetReg = qRetReg * airDensityIN - mHRV_AH;
									mechVentPower = mechVentPower + fanPower_cooling;		// note AH operates at cooling speed
									AHfanHeat = fanPower_cooling * .85;						// for heat
								} else {													// open the outside air vent
									mHRV_AH = fan[i].q * airDensityIN;
									mRetReg = qRetReg * airDensityIN - mHRV_AH;
								}
							} else {
								fan[i].on = 0;
								mHRV_AH = 0;
							}
						}

						//From the HVI directory (as of 5/2015) of all ERV units, average SRE=0.63, LR/MT=0.51, TRE=0.52, cfm/w=1.1

						// ERV + Air Handler Copied over from ARTI code
						if(fan[i].oper == 17) {												
					
							if((rivecFlag == 0 && minute >= 40) || (rivecFlag == 1 && rivecOn == 1)) {	// on for last 20 minutes of every hour or RIVEC controlled
								fan[i].on = 1;
								ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
								ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								mechVentPower = mechVentPower + fan[i].power;
								if(AHflag != 1 && AHflag != 2 && AHflag != 102) {			// AH is off for heat/cool so turn it on for venting
									AHflag = 100;
									AHminutes = AHminutes + 1;
									qAH = qAH_cool;
									supVelAH = qAH / (pow(supDiameter,2) * M_PI / 4);
									retVelAH = qAH / (pow(retDiameter,2) * M_PI / 4);
									qSupReg = -qAH * supLF + qAH;
									qRetReg = qAH * retLF - qAH;
									qRetLeak = -qAH * retLF;
									qSupLeak = qAH * supLF;
									mSupLeak = qSupLeak * airDensityIN;
									mRetLeak = qRetLeak * airDensityIN;
									mAH = qAH * airDensityIN;
									mERV_AH = abs(fan[i].q * airDensityIN) * -1.0;
									mSupReg = qSupReg * airDensityIN;
									mRetReg = qRetReg * airDensityIN - mERV_AH;
									mechVentPower = mechVentPower + fanPower_cooling;		// note ah operates at cooling speed
									AHfanHeat = fanPower_cooling * .85;						// for heat
								} else {													// open the outside air vent
									mERV_AH = fan[i].q * airDensityIN;						// mass flow of ERV unit synced to Air Handler
									mRetReg = qRetReg * airDensityIN - mERV_AH;
								}
							} else {
								fan[i].on = 0;
								mERV_AH = 0;
							}
						}
					}
					// [END] HRV's and ERV's =================================================================================================================================


					// [START] Auxiliary Fan Controls ============================================================================================================

					for(int i = 0; i < numFans; i++) {

						if(fan[i].oper == 1 && OccContType > 1) {		// RIVEC controlled exhaust fan.
							if(rivecOn == 1){
								fan[i].on = 1;
								//fan[i].q = FanQ;
								//fan[i].on = 0;	// Use this to disable the whole-house exhaust fan
								mechVentPower = mechVentPower + fan[i].power;				// vent fan power
								if(fan[i].q > 0) {							// supply fan - its heat needs to be added to the internal gains of the house
									fanHeat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
									ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
								}
								else											// exhasut fan
									ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;	
							} else{
								fan[i].on = 0;
								//fan[i].q = 0;
							}	
						
						} else if(fan[i].oper == 1 && OccContType <= 1){ //Continuous exhaust fan simulations.
							fan[i].on = 1;
							//fan[i].q = FanQ;
							mechVentPower = mechVentPower + fan[i].power;
							if(fan[i].q > 0) {							// supply fan - its heat needs to be added to the internal gains of the house
								fanHeat = fan[i].power * .84;			// 16% efficiency for the particular fan used in this study.
								ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
							}
							else											// exhasut fan
								ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;		
						
						} else if(fan[i].oper == 2) {					// FIXED SCHEDULE BATHROOM
							if(hour == 7 && minute >= 30) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 3) {					// FIXED SCHEDULE KITCHEN
							if(hour > 17 && hour <= 18) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;			// vent fan power
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 4) {					// FIXED SCHEDULE EXHAUST FAN
							if(hcFlag == 1) {  // if heating
								if(hour > 1 && hour <= 5)
									fan[i].on = 0;
								else {
									fan[i].on = 1;
									mechVentPower = mechVentPower + fan[i].power;		// vent fan power
									//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
									nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								}
							}
							if(hcFlag == 2) {  // if cooling
								if(hour > 15 && hour <= 19) {
									fan[i].on = 0;
								} else {
									fan[i].on = 1;
									mechVentPower = mechVentPower + fan[i].power;		// vent fan power
									//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
									nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								}
							}

						} else if(fan[i].oper == 19) {					// DRYER FAN
							if(weekend == 1) {						// Three hours of operation, two days per week
								if(hour > 12 && hour <= 15) {
									fan[i].on = 1;
									mechVentPower = mechVentPower + fan[i].power;
									//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
									nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								} else
									fan[i].on = 0;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 21) {					// Economizer
							if(econoFlag == 1) {
								fan[i].on = 1;
								AHfanPower = AHfanPower + fan[i].power;			// vent fan power
								//ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumIN = nonRivecVentSumIN + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 23) {					// DYNAMIC SCHEDULE DRYER
							if(dryerFan == 1) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 24) {					// DYNAMIC SCHEDULE KITCHEN FAN
							if(kitchenFan == 1) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 25) {					// DYNAMIC SCHEDULE BATHROOM 1
							if(bathOneFan == 1) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 26) {					// DYNAMIC SCHEDULE BATHROOM 2
							if(bathTwoFan == 1) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;

						} else if(fan[i].oper == 27) {					// DYNAMIC SCHEDULE BATHROOM 3. Brennan, we can redirect this fan.oper=27 to look at the bathOneFan value to increaes usage for 6 occupants.
							if(bathThreeFan == 1) {
							//if(bathOneFan == 1) {
								fan[i].on = 1;
								mechVentPower = mechVentPower + fan[i].power;
								//ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
								nonRivecVentSumOUT = nonRivecVentSumOUT + abs(fan[i].q) * 3600 / houseVolume;
							} else
								fan[i].on = 0;
							
								
//sum the non-Rivec ventilation airflows (prior section) and set nonRivecVentSum to the largest of inflow/outflow.								
								
						//} else if(fan[i].oper == 30) {// Continuous, Infiltration-Balanced Ventilation Temperature Controller ("Gold Standard"). Brennan.
						//	if(numStories == 1){
						//		AIM2 = pow(abs((tempHouse - 273.15) - weatherTemp) , 0.67) * 0.069 * C; //estimated stack infiltration (m3/s). Alwasy positive.
						//		AEQaim2FlowDiff = (Aeq * houseVolume / 3600) - AIM2; //difference between 62.2-2013 AEQ flow rate (m3/s) and the estimated infiltration. Sometimes negative.
						//		qFanFlowRatio = AEQaim2FlowDiff / (Aeq * houseVolume / 3600); //dimensionless ratio of flows. Sometimes negative.
						//			if(AIM2 < (Aeq * houseVolume / 3600)){ //tests if infiltration exceeds Aeq
						//				fan[i].on = 1;
						//				fan[i].q = FanQ * qFanFlowRatio; //adjusts fan flow based on ratio of AIM2 / Aeq
						//				fan[i].power = FanP * qFanFlowRatio; //adjusts fan power based on ratio of AIM2 / Aeq
						//				mechVentPower = mechVentPower + fan[i].power ;
						//				ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume ;  
						//			} else
						//				fan[i].on = 0;
						//	} else 
						//		if(numStories == 2){
						//			AIM2 = pow(abs((tempHouse - 273.15) - weatherTemp) , 0.67) * 0.089 * C;
						//			AEQaim2FlowDiff = (Aeq * houseVolume / 3600) - AIM2;
						//			qFanFlowRatio = AEQaim2FlowDiff / (Aeq * houseVolume / 3600);
						//				if(AIM2 < (Aeq * houseVolume / 3600)){
						//					fan[i].on = 1;
						//					fan[i].q = FanQ * qFanFlowRatio;
						//					fan[i].power = FanP * qFanFlowRatio;
						//					mechVentPower = mechVentPower + fan[i].power ;
						//					ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume ; 
						//				} else 
						//					fan[i].on = 0;
						//	} 


						//} else if(fan[i].oper == 31) {// Continuous, Infiltration-Balanced Ventilation Temperature Controller ("Gold Standard"). Brennan. Quadrature!
						//	if(numStories == 1){
						//		AIM2 = 0.069 * envC * pow(abs((tempHouse - 273.15) - weatherTemp), n); //estimated stack infiltration (m3/s).
						//		AEQaim2FlowDiff = pow((Aeq * houseVolume / 3600),2) - pow(AIM2,2); //difference between squares of 62.2-2013 AEQ flow rate (m3/s) and the estimated infiltration. Tells us if AIM2 inf > Aeq
						//		if(AEQaim2FlowDiff > 0){ 
						//			fan[i].on = 1;
						//			qFanFlowRatio = sqrt(AEQaim2FlowDiff) / (Aeq * houseVolume / 3600); //sqrt of the squared differences solves for the required fan flow, this is the ratio of the req fan flow to the Aeq
						//			fan[i].q = FanQ * qFanFlowRatio;
						//			fan[i].power = FanP * qFanFlowRatio;
						//			mechVentPower = mechVentPower + fan[i].power ;
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//		} else
						//			fan[i].on = 0;
						//	} else 
						//			if(numStories == 2){
						//			AIM2 = 0.089 * envC * pow(abs((tempHouse - 273.15) - weatherTemp), n); 
						//			AEQaim2FlowDiff = pow((Aeq * houseVolume / 3600),2) - pow(AIM2,2); 
						//				if(AEQaim2FlowDiff > 0){
						//					fan[i].on = 1;
						//					qFanFlowRatio = sqrt(AEQaim2FlowDiff) / (Aeq * houseVolume / 3600); //
						//					fan[i].q = FanQ * qFanFlowRatio;
						//					fan[i].power = FanP * qFanFlowRatio;
						//					mechVentPower = mechVentPower + fan[i].power ;
						//					ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//				} else
						//					fan[i].on = 0;
								


							} 
							
						if(nonRivecVentSumIN > nonRivecVentSumOUT){						//ventSum based on largest of inflow or outflow
							nonRivecVentSum = nonRivecVentSumIN;
						} else{
							nonRivecVentSum = nonRivecVentSumOUT;			
						}	
						
						}
					//}
					// [END] Auxiliary Fan Controls ========================================================================================================

			
						// [END] ---------------------FAN 50---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v6

						// ---------------------FAN 30---------------------- OLD RIVEC OPERATION BASED ON CONTROL ALGORITHM v1
						////if(fan[i].oper == 30) {
						////	if(minute == 1 || minute == 11 || minute == 21 || minute == 31 || minute == 41 || minute == 51) {
						////		rivecOn = 1;		// start with controlled fan on as default					
						////		if(hcFlag == 1) {									// HEATING TIMES
						////			if(hour >= 10 && hour < 22) {					// base time period
						////				if(relDose <= 1 && relExp <= 0.8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum < Aeq)
						////					rivecOn = 1;
						////			} else if(hour >= 22 || hour < 2) {				// pre-peak time period
						////				if(relDose <= 1 && relExp <= 0.8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum > 1.25 * Aeq)
						////					rivecOn = 0;
						////			} else if(hour >= 2 && hour < 6) {				// peak time period
						////				rivecOn = 0;								// always off
						////			} else {                                        // post-peak time period
						////				if(relDose <= 1 && relExp <= 0.8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum > 1.25 * Aeq)
						////					rivecOn = 0;
						////			}
						////		} else {											// COOLING TIMES
						////			if(hour >= 22 || hour < 10) {					// base time
						////				if(relDose <= 1 && relExp <= .8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum < Aeq)
						////					rivecOn = 1;
						////			} else if(hour >= 10 && hour < 14) { 			// prepeak
						////				if(relDose <= 1 && relExp <= 0.8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum > 1.25 * Aeq)
						////					rivecOn = 0;
						////			} else if(hour >= 14 && hour < 18) {			// peak
						////				rivecOn = 0;								// always off
						////			} else {										// post peak
						////				if(relDose <= 1 && relExp <= 0.8)
						////					rivecOn = 0;
						////				if(nonRivecVentSum > 1.25 * Aeq)
						////					rivecOn = 0;
						////			}

						////		}
						////	}

						//	//RIVEC decision after all other fans added up
						//	if(rivecOn == 0)										//RIVEC has turned off this fan
						//		fan[i].on = 0;
						//	else {													
						//		fan[i].on = 1;
						//		rivecMinutes++;
						//		mechVentPower = mechVentPower + fan[i].power;		// vent fan power
						//		if(fan[i].q > 0) {									// supply fan - its heat needs to be added to the internal gains of the house
						//			fanHeat = fan[i].power * .84;					// 16% efficiency for the particular fan used in this study.
						//			ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
						//		} else												// ex fan
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//	}
						//}
						//---------------------FAN 30 end----------------------

						//// ---------------------FAN 31---------------------- OLD RIVEC OPERATION BASED ON CONTROL ALGORITHM v2
						//if(fan[i].oper == 31) {
						//	if(minute == 1 || minute == 11 || minute == 21 || minute == 31 || minute == 41 || minute == 51) {
						//		rivecOn = 1;		// start with controlled fan on as default					
						//		if(hcFlag == 1) {									// HEATING TIMES
						//			if(hour >= 10 && hour < 22) {					// base time period
						//				if(relDose <= 1 && relExp <= .8)
						//					rivecOn = 0;
						//				//if(nonRivecVentSum < Aeq)
						//				//rivecOn = 1; (removed from revised algorithm)
						//			} else if(hour >= 22 || hour < 2) {				// pre-peak time period
						//				if(relDose <= 1 && relExp <= 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			} else if(hour >= 2 && hour < 6) {				// peak time period
						//				rivecOn = 0;								// always off
						//			} else {                                        // post-peak time period
						//				if(relDose <= 1 && relExp <= 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}
						//		} else {											// COOLING TIMES
						//			if(hour >= 22 || hour < 10) {					// base time
						//				if(relDose <= 1 && relExp <= .8)
						//					rivecOn = 0;
						//				//IF nonRivecVentSum < Aeq THEN rivecOn = 1 (removed from revised algorithm)
						//			} else if(hour >= 10 && hour < 14) { 			// prepeak
						//				if(relDose <= 1 && relExp <= 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			} else if(hour >= 14 && hour < 18) {			// peak
						//				rivecOn = 0;								// always off
						//			} else {										// post peak
						//				if(relDose <= 1 && relExp <= 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}

						//		}
						//	}

						//	//RIVEC decision after all other fans added up
						//	if(rivecOn == 0)										//RIVEC has turned off this fan
						//		fan[i].on = 0;
						//	else {													
						//		fan[i].on = 1;
						//		rivecMinutes++;
						//		mechVentPower = mechVentPower + fan[i].power;		// vent fan power
						//		if(fan[i].q > 0) {									// supply fan - its heat needs to be added to the internal gains of the house
						//			fanHeat = fan[i].power * .84;					// 16% efficiency for the particular fan used in this study.
						//			ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
						//		} else												// ex fan
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//	}
						//}
						//---------------------FAN 31 end----------------------

						//---------------------FAN 50---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v5
						//if(fan[i].oper == 50) {
						//	if(minute == 1 || minute == 11 || minute == 21 || minute == 31 || minute == 41 || minute == 51) {
						//	// rivecOn = 1 or 0: 1 = whole-house fan ON, 0 = whole-house fan OFF
						//	rivecOn = 0;
						//		if(hour >= baseStart || hour < baseEnd ) {        		// BASE Time Period
						//			
						//			if(occupied[hour]) {				            	// Base occupied
						//				if(relExp > 0.8 && nonRivecVentSum < Aeq) 
						//					rivecOn = 1;

						//			} else {						                	// Base unoccupied
						//				if(relExp > expLimit && nonRivecVentSum < 1.25 * Aeq)
						//					rivecOn = 1;
						//			}

						//		} else if(hour >= peakStart && hour < peakEnd) {		// PEAK Time Period
						//			rivecOn = 0;										// Always off

						//		} else if(hour >= recoveryStart && hour < recoveryEnd) {// RECOVERY Time Period
						//	
						//			if (occupied[hour]) {								// Recovery occupied
						//				if(relExp > 1 && nonRivecVentSum < 1.25 * Aeq)
						//					rivecOn = 1;

						//			} else {											// Recovery unoccupied
						//				if(relExp > expLimit && nonRivecVentSum < 1.25 * Aeq)
						//					rivecOn = 1;
						//			}
						//		}
						//	}

						//	// RIVEC fan operation after algorithm decision
						//	if(rivecOn) {	  												// RIVEC has turned ON this fan
						//		fan[i].on = 1;
						//		mechVentPower = mechVentPower + fan[i].power;				// vent fan power
						//		if(fan[i].q > 0) { 											// supply fan - its heat needs to be added to the internal gains of the house
						//			fanHeat = fan[i].power * .84;							// 16% efficiency for this fan
						//			ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
						//		} else { 													// exhaust fan
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//		}
						//		rivecMinutes++;
						//	}
						//	else {  														// limited set of fans controlled by RIVEC. Only supply, exhaust and HRV are controlled
						//		fan[i].on = 0;

						//	}
						//}
						// --------------------FAN 50 end----------------------


						////---------------------FAN 50---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v4 WJNT May 2011
						//if(fan[i].oper == 50) {
						//	// if(minute == 1 || minute == 11 || minute == 21 || minute == 31 || minute == 41 || minute == 51) {
						//	if(minute == 1 || minute == 16 || minute == 31 || minute == 46) {

						//		// rivecOn = 1 or 0: 1 = whole-house fan ON, 0 = whole-house fan OFF
						//		if(hcFlag == 1) {      										// HEATING TIMES
						//			if((hour >= 0 && hour < 4) || (hour >= 12)) {        	// BASE Time Period
						//				if(occupied[hour] == 1) {                        	// Occupied
						//					rivecOn = 1;									// Turn rivec fan ON
						//					if(relDose < 1 && relExp < .8) {
						//						if(nonRivecVentSum > Aeq) {
						//							rivecOn = 0;
						//						}											
						//					}									
						//				} else if(occupied[hour] == 0) {                	// Unoccupied
						//					//if(relDose > 1.2 || relExp > expLimit)
						//					if(relDose > 1.2 || relExp > 2)
						//						rivecOn = 1;
						//					//if(relDose < 1 && relExp < expLimit)
						//					if(relDose < 1 && relExp < 2)
						//						rivecOn = 0;
						//					if(nonRivecVentSum > Aeq)
						//						rivecOn = 0;									
						//				}
						//			} else if(hour >= 4 && hour < 8) {						// PEAK Time Period
						//				rivecOn = 0;										// Always off
						//			} else if(hour >= 8 && hour < 12) {         			// RECOVERY Time Period
						//				if(occupied[hour] == 0)
						//					rivecOn = 0;
						//				if(nonRivecVentSum < 1.25 * Aeq) {
						//					if(relExp > expLimit)
						//					//if(relExp > 1)
						//						rivecOn = 1;
						//				}
						//				if(relExp < 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}
						//		} else {    												// COOLING TIMES
						//			if(hour >= 22 || (hour >= 0 && hour < 14)) {           	// BASE Time Period
						//				if(occupied[hour] == 1) {                           // Occupied
						//					rivecOn = 1;									// Turn rivec fan ON
						//					if(relDose < 1 && relExp < .8) {
						//						if(nonRivecVentSum > Aeq)
						//							rivecOn = 0;
						//					}
						//				} else if(occupied[hour] == 0) {	             	// Unoccupied
						//					if(relDose > 1.2 || relExp > expLimit)
						//					//if(relDose > 1.2 || relExp > 2)
						//						rivecOn = 1;
						//					if(relDose < 1 && relExp < expLimit)
						//					//if(relDose < 1 && relExp < 2)
						//						rivecOn = 0;
						//					if(nonRivecVentSum > Aeq)
						//						rivecOn = 0;
						//				}
						//			} else if(hour >= 14 && hour < 18) {                 	// PEAK Time Period
						//				rivecOn = 0;                                		// Always off
						//			} else if(hour >= 18 && hour < 22) {         			// RECOVERY Time Period
						//				if(occupied[hour] == 0)
						//					rivecOn = 0;
						//				if(nonRivecVentSum < 1.25 * Aeq) {
						//					if(relDose > 1 || relExp > expLimit)
						//					//if(relDose > 1 || relExp > 1)
						//						rivecOn = 1;
						//				}
						//				if(relExp < 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}
						//		}
						//	}
						//	// RIVEC decision after all other fans added up
						//	if(rivecOn == 0)	  											// RIVEC has turned off this fan
						//		fan[i].on = 0;
						//	else {  														// limited set of fans controlled by RIVEC. Only supply, exhaust and HRV are controlled
						//		rivecMinutes++;
						//		fan[i].on = 1;
						//		mechVentPower = mechVentPower + fan[i].power;				// vent fan power
						//		if(fan[i].q > 0) { 											// supply fan - its heat needs to be added to the internal gains of the house
						//			fanHeat = fan[i].power * .84;							// 16% efficiency for the particular fan used in this study.
						//			ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
						//		} else { 													// exhaust fan
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//		}
						//	}
						//}
						//// --------------------FAN 50 end----------------------

						//---------------------FAN 51---------- RIVEC OPERATION BASED ON CONTROL ALGORITHM v3 WJNT May 2011
						//if(fan[i].oper == 51) {
						//	if(minute == 1 || minute == 11 || minute == 21 || minute == 31 || minute == 41|| minute == 51) {
						//		if(hcFlag == 1) {      										// HEATING TIMES
						//			if((hour >= 0 && hour < 4) || (hour >= 12)) {        	// BASE Time Period
						//				if(occupied[hour] == 1) {                        	// Occupied
						//					rivecOn = 1;									// Turn rivec fan ON
						//					if(relDose < 1 && relExp < .8) {
						//						if(nonRivecVentSum > Aeq) {
						//							rivecOn = 0;
						//						}											
						//					}									
						//				} else if(occupied[hour] == 0) {                	// Unoccupied
						//					if(relDose > 1.2 || relExp > 2)
						//						rivecOn = 1;
						//					if(relDose < 1 && relExp < 2)
						//						rivecOn = 0;
						//					if(nonRivecVentSum > Aeq)
						//						rivecOn = 0;									
						//				}
						//			} else if(hour >= 4 && hour < 8) {						// PEAK Time Period
						//				rivecOn = 0;										// Always off
						//			} else if(hour >= 8 && hour < 12) {         			// RECOVERY Time Period
						//				if(nonRivecVentSum < 1.25 * Aeq) {
						//					if(relDose > 1 || relExp > 1)
						//						rivecOn = 1;
						//				}
						//				if(relDose < 1 && relExp < 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}
						//		} else {    												// COOLING TIMES
						//			if(hour >= 22 || (hour >= 0 && hour < 14)) {           	// BASE Time Period
						//				if(occupied[hour] == 1) {                           // Occupied
						//					rivecOn = 1;									// Turn rivec fan ON
						//					if(relDose < 1 && relExp < .8) {
						//						if(nonRivecVentSum > Aeq)
						//							rivecOn = 0;
						//					}
						//				} else if(occupied[hour] == 0) {					// Unoccupied
						//					if(relDose > 1.2 || relExp > 2)
						//						rivecOn = 1;
						//					if(relDose < 1 && relExp < 2)
						//						rivecOn = 0;
						//					if(nonRivecVentSum > Aeq)
						//						rivecOn = 0;
						//				}
						//			} else if(hour >= 14 && hour < 18) {                 	// PEAK Time Period
						//				rivecOn = 0;                                		// Always off
						//			} else if(hour >= 18 && hour < 22) {         			// RECOVERY Time Period
						//				if(nonRivecVentSum < 1.25 * Aeq) {
						//					if(relDose > 1 || relExp > 1)
						//						rivecOn = 1;
						//				}
						//				if(relDose < 1 && relExp < 1)
						//					rivecOn = 0;
						//				if(nonRivecVentSum > 1.25 * Aeq)
						//					rivecOn = 0;
						//			}
						//		}
						//	}
						//	// RIVEC decision after all other fans added up
						//	if(rivecOn == 0)	  											// RIVEC has turned off this fan
						//		fan[i].on = 0;
						//	else {  														// limited set of fans controlled by RIVEC. Only supply, exhaust and HRV are controlled
						//		rivecMinutes++;
						//		fan[i].on = 1;
						//		mechVentPower = mechVentPower + fan[i].power;								// vent fan power
						//		if(fan[i].q > 0) { 											// supply fan - its heat needs to be added to the internal gains of the house
						//			fanHeat = fan[i].power * .84;							// 16% efficiency for the particular fan used in this study.
						//			ventSumIN = ventSumIN + abs(fan[i].q) * 3600 / houseVolume;
						//		} else { 													// exhaust fan
						//			ventSumOUT = ventSumOUT + abs(fan[i].q) * 3600 / houseVolume;
						//		}
						//	}
						//}
						// ---------------------FAN 51 end----------------------


					//} // [End RIVEC Fans] =======================================================================================================================================



					// [START] Hybrid Systems ================================================================================================================================
					// Block off passive stacks/flues when the RIVEC fan is running for hybrid ventilation strategy
			
					//if(rivecOn)
					//	numFlues = 0;
					//else
					//	numFlues = numFluesActual;

					// [END] Hybrid Systems ==================================================================================================================================

					// [START] Equipment Model ===============================================================================================================================
					
					evapcap = 0;
					latcap = 0;
					capacity = 0;
					capacityh = 0;
					compressorPower = 0;
					

					if(AHflag == 0) {										// AH OFF
						capacityc = 0;
						capacityh = 0;
						compressorPower = 0;
					} else if(AHflag == 100) {								// Fan on no heat/cool
						capacityh = AHfanHeat;								// Include fan heat
						capacityc = 0;
					} else if(AHflag == 102) {								// End if heat cycle
						capacityh = .25 * hcapacity + AHfanHeat;			// NOTE 0.25 burner capacity at 1/4 during cool down
						capacityc = 0;
					} else if(AHflag == 1) {								// Heating
						capacityh = hcapacity + AHfanHeat;					// include fan heat in furnace capacity (converted to W in main program)
						capacityc = 0;
					} else {												// we have cooling
						double Toutf = KtoF(weather.dryBulb);

						// the following corrections are for TXV only
						// test for wet/dry coil using SHR calculation
						SHR = 1 - 50 * (HRReturn - .005);										// using humidity ratio - note only because we have small indoor dry bulb range

						if(SHR < .25)
							SHR = .25;  //setting low limit
						if(SHR > 1)
							SHR = 1;
						if(SHR == 1) { // dry coil
							// charge correction for capacity (dry coil)
							if(charge < .725) {
								chargecapd = 1.2 + (charge - 1);
							} else if(charge >= .725) {
								chargecapd = .925;
							}
							// charge correction for EER (dry coil)
							if(charge < 1) {
								chargeeerd = 1.04 + (charge - 1) * .65;
							} else if(charge >= 1) {
								chargeeerd = 1.04 - (charge - 1) * .35;
							}
							capacity = (capacityari * 1000) * .91 * qAHcorr * chargecapd * ((-.00007) * pow((Toutf - 82),2) - .0067 * (Toutf - 82) + 1);
							// EER for dry coil
							EER = EERari * qAHcorr * chargeeerd * ((-.00007) * pow((Toutf - 82),2) - .0085 * (Toutf - 82) + 1);
						} else {															// wet coil
							// charge correction for capacity (wet coil)
							if(charge < .85) {
								chargecapw = 1 + (charge - .85);
							} else if(charge >= .85) {
								chargecapw = 1;
							}
							// charge correction for EER (wet coil)
							if(charge < .85) {
								chargeeerw = 1 + (charge - .85) * .9;
							} else if(charge >= .85 && charge <= 1) {
								chargeeerw = 1;
							} else if(charge > 1) {
								chargeeerw = 1 - (charge - 1) * .35;
							}
							double hret = .24 * ((tempReturn - 273.15) * 9 / 5 + 32) + HRReturn * (1061 + .444 * ((tempReturn - 273.15) * 9 / 5 + 32));   // return air enthalpy
							hret += AHfanHeat / mAH / 2326.;										// add fan heat, converted to Btu/lb
							capacity = capacityari * 1000 * (1 + (hret - 30) * .025) * qAHcorr * chargecapw * ((-.00007) * pow((Toutf - 95),2) - .0067 * (Toutf - 95) + 1);
							// EER for wet coil
							EER = EERari * qAHcorr * chargeeerw * ((-.00007) * pow((Toutf - 95),2) - .0085 * (Toutf - 95) + 1);
						}

						// split sensible and latent capacities and convert capacity Btu/h to W
						// calculation of power used and total cacapity for air conditioning
						compressorPower = capacity / EER; //compressor power

						// SHR ramps up for first three minutes of operation
						// compTime tracks the number of minutes the compressor is on for
						if(compTime == 1) {
							SHR = SHR + 2.0 / 3.0 * (1 - SHR);
						} else if(compTime == 2) {
							SHR = SHR + 1.0 / 3.0 * (1 - SHR);
						}
						
						// Reduce capacity over cCapAdjustTime minutes (to reduce sensible spike)
						if(compTime > 0 && compTime < cCapAdjustTime) {
							capacity *= (compTime / cCapAdjustTime);
							}

						capacityc = SHR * capacity / 3.413;									// sensible (equals total if SHR=1)
						// correct the sensible cooling for fan power heating
						capacityc = capacityc - AHfanHeat;
						//cout << "minute" << minute << " comptime=" << compTime << " compTimeCount=" << compTimeCount << "capc=" << capacityc << " SHR=" << SHR << endl;
						// tracking mass on coil (Mcoil) including condensation and evaporation - do this in main program
						evapcap = 0;
						latcap = 0;
					}

					// Tracking Coil Moisture
					if(AHflag == 2) {
						if(SHR < 1) {
							latcap = (1 - SHR) * capacity / 3.413;							// latent capacity J/s
							Mcoil = Mcoilprevious + latcap / 2501000 * dtau;				// condenstion in timestep kg/s * time
						} else {
							if(Mcoil > 0) {													// if moisture on coil but no latcap, then we evaporate coil moisture until Mcoil = zero
								latcap = -.3 * capacityraw / 1800 * 2501000;				// evaporation capacity J/s - this is negative latent capacity
								evapcap = latcap;
								Mcoil = Mcoilprevious - .3 * capacityraw / 1800 * dtau;		// evaporation in timestep kg/s * time
							} else {														// no latcap and dry coil
								latcap = 0;
							}
						}
					} else if(AHflag == 100) {												// then we have no cooling, but the fan is on
						if(Mcoil > 0) {														// if moisture on coil but no latcap, then we evaporate coil moisture until Mcoil = zero
							latcap = -.3 * capacityraw / 1800 * 2501000;					// evaporation capacity J/s - this is negative latent capacity
							evapcap = latcap;
							Mcoil = Mcoilprevious - .3 * capacityraw / 1800 * dtau;			// evaporation in timestep kg/s * time
						} else {															// no latcap and dry coil
							latcap = 0;
						}
					}

					if(Mcoil < 0)
						Mcoil = 0;
					if(Mcoil > .3 * capacityraw)
						Mcoil = .3 * capacityraw;											// maximum mass on coil is 0.3 kg per ton of cooling
					Mcoilprevious = Mcoil;													// maybe put this at top of the hour

					if(dhCapacity > 0) {
						dh.run(RHHouse, tempHouse);	// run dehumidifier using house air node conditions
						}
					
					// [END] Equipment Model ======================================================================================================================================

					// [START] Heat and Mass Transport ==============================================================================================================================
					mCeilingOld = -1000;														// inital guess
					mainIterations = 0;
					limit = envC / 10;
					if(limit < .00001)
						limit = .00001;

					// Ventilation and heat transfer calculations
					while(1) {
						mainIterations = mainIterations + 1;	// counting # of temperature/ventilation iterations
						int flag = 0;

						while(1) {
							// Call houseleak subroutine to calculate air flow. Brennan added the variable mCeilingIN to be passed to the subroutine. Re-add between mHouseIN and mHouseOUT
							sub_houseLeak(AHflag, flag, weather.windSpeed, weather.windDirection, tempHouse, tempAttic, weather.dryBulb, envC, envPressureExp, eaveHeight,
								leakFracCeil, leakFracFloor, leakFracWall, numFlues, flue, wallFraction, floorFraction, Sw, flueShelterFactor, numWinDoor, winDoor, numFans, fan, numPipes,
								Pipe, mIN, mOUT, Pint, mFlue, mCeiling, mFloor, atticC, dPflue, Crawl,
								Hfloor, rowHouse, soffitFraction, Patticint, wallCp, mSupReg, mAH, mRetLeak, mSupLeak,
								mRetReg, mHouseIN, mHouseOUT, supC, supn, retC, retn, mSupAHoff, mRetAHoff, airDensityIN, airDensityOUT, airDensityATTIC, houseVolume, windPressureExp);
							//Yihuan : put the mCeilingIN on comment 
							flag = flag + 1;

							if(abs(mCeilingOld - mCeiling) < limit || flag > 5)
								break;
							else
								mCeilingOld = mCeiling;

							// call atticleak subroutine to calculate air flow to/from the attic
							sub_atticLeak(flag, weather.windSpeed, weather.windDirection, tempHouse, weather.dryBulb, tempAttic, atticC, atticPressureExp, eaveHeight, roofPeakHeight,
								flueShelterFactor, Sw, numAtticVents, atticVent, soffit, mAtticIN, mAtticOUT, Patticint, mCeiling, rowHouse,
								soffitFraction, roofPitch, roofPeakPerpendicular, numAtticFans, atticFan, mSupReg, mRetLeak, mSupLeak, matticenvin,
								matticenvout, mSupAHoff, mRetAHoff, airDensityIN, airDensityOUT, airDensityATTIC);
						}

						// adding fan heat for supply fans, internalGains1 is from input file, fanHeat reset to zero each minute, internalGains is common
						internalGains = internalGains1 + fanHeat;

						//bsize = sizeof(b)/sizeof(b[0]);

						// Call heat subroutine to calculate heat exchange
						sub_heat(weather.dryBulb, mCeiling, AL4, weather.windSpeed, ssolrad, nsolrad, tempOld, atticVolume, houseVolume, weather.skyCover, b, ERRCODE, TSKY,
							floorArea, roofPitch, ductLocation, mSupReg, mRetReg, mRetLeak, mSupLeak, mAH, supRval, retRval, supDiameter,
							retDiameter, supArea, retArea, supThickness, retThickness, supVolume, retVolume, supCp, retCp, supVel, retVel, suprho,
							retrho, weather.pressure, weather.humidityRatio, uaSolAir, uaTOut, matticenvin, matticenvout, mHouseIN, mHouseOUT, planArea, mSupAHoff,
							mRetAHoff, solgain, tsolair, mFanCycler, roofPeakHeight, eaveHeight, retLength, supLength,
							roofType, M1, M12, M15, M16, roofRval, rceil, AHflag, mERV_AH, ERV_SRE, mHRV, HRV_ASE, mHRV_AH,
							capacityc, capacityh, evapcap, internalGains, airDensityIN, airDensityOUT, airDensityATTIC, airDensitySUP, airDensityRET, numStories, storyHeight,
							dh.sensible, H2, H4, H6);

						if((abs(b[0] - tempAttic) < .2) || (mainIterations > 10)) {	// Testing for convergence
							tempAttic        = b[0];					
							tempReturn       = b[11];
							tempSupply       = b[14];
							tempHouse        = b[15];
							break;
						}

						tempAttic = b[0];
						tempHouse = b[15];
					}

					// setting "old" temps for next timestep to be current temps:
					// [START] Moisture Balance ===================================================================================================================================

					if(moistureModel == 1) {
						// Call moisture balance
						double mRetOut = mFanCycler + mHRV_AH + mERV_AH * (1 - ERV_TRE);
						moisture_nodes.mass_cond_bal(b, weather.dryBulb, weather.relativeHumidity,
						airDensityOUT, airDensityATTIC, airDensityIN, airDensitySUP, airDensityRET,
						weather.pressure, H4, H2, H6, matticenvin, matticenvout, mCeiling, mHouseIN, mHouseOUT,
						mAH, mRetAHoff, mRetLeak, mRetReg, mRetOut, mERV_AH * ERV_TRE, mSupAHoff, mSupLeak, mSupReg,
						latcap, dh.condensate, latentLoad);

						HRAttic = calcHumidityRatio(moisture_nodes.PW[6],weather.pressure);
						HRReturn = calcHumidityRatio(moisture_nodes.PW[7],weather.pressure);  
						HRSupply = calcHumidityRatio(moisture_nodes.PW[8],weather.pressure);
						HRHouse = calcHumidityRatio(moisture_nodes.PW[9],weather.pressure);  
						RHHouse = moisture_nodes.moistureContent[9];
						RHAttic = moisture_nodes.moistureContent[6];
						}
					else {
						// Call moisture subroutine
						sub_moisture(HR, M1, M12, M15, M16, Mw5, matticenvout, mCeiling, mSupAHoff, mRetAHoff,
							matticenvin, weather.humidityRatio, mSupLeak, mAH, mRetReg, mRetLeak, mSupReg, latcap, mHouseIN, mHouseOUT,
							latentLoad, mFanCycler, mHRV_AH, mERV_AH, ERV_TRE, MWha, airDensityIN, airDensityOUT, dh.condensate);

						//Calculate Saturation Humidity Ratio, Equation 23 in ASHRAE HoF
						double SatVaporPressure = saturationVaporPressure(tempHouse);
						double HRsaturation = calcHumidityRatio(SatVaporPressure, weather.pressure); 
						if(HR[1] == 0)
							HR[1] = weather.humidityRatio;
						if(HR[3] > HRsaturation) // Previously set by Iain to 0.02. Here we've replaced it with the saturation humidity ratio (Ws).
							HR[3] = HRsaturation; //consider adding calculate saturation humidity ratio by indoor T and Pressure and lmit HR[3] to that.
						RHHouse = 100 * ((weather.pressure*(HR[3]/0.621945))/(1+(HR[3]/0.621945)) / SatVaporPressure);
						SatVaporPressure = saturationVaporPressure(tempAttic);
						RHAttic = 100 * ((weather.pressure*(HR[0]/0.621945))/(1+(HR[0]/0.621945)) / SatVaporPressure);
						HRAttic = HR[0];
						HRReturn = HR[1];
						HRSupply = HR[2];
						HRHouse = HR[3];
						}

					// [END] Moisture Balance =======================================================================================================================================

					for(int i = 0; i < ATTIC_NODES; i++) {
						tempOld[i]  = b[i];
						}
						
						
		

					// ************** house ventilation rate  - what would be measured with a tracer gas i.e., not just envelope and vent fan flows
					// mIN has msupreg added in mass balance calculations and mRetLeak contributes to house ventilation rate

					mHouse = mIN - mRetLeak * (1 - supLF) - mSupReg; //Brennan, I think this is what we should use. It's already calculated! For this to equal 0, either all values evaluate to 0, or mIN exactly equals the sum of the other values, first seems likely.
					mHouse = mHouse - mFanCycler - mERV_AH - mHRV_AH;

					qHouse = abs(mHouse / airDensityIN);			// Air flow rate through the house [m^3/s]. Brennan. These suddenly compute to 0, at minute 228,044, which means mHouse went to 0. I think. We get very werid flucations in house pressure and qHouse leading up to the error, where pressure cycles minute-by-minute between 4 and almost 0, and then eventually evaluates to actual 0.
					houseACH = qHouse / houseVolume * 3600;			// Air Changes per Hour for the whole house [h-1]. Brennan. These suddnely compute to 0, at minute 228,044

					if(mFlue < 0)												// mFlue is the mass airflow through all of the flues combined (if any)
						flueACH = (mFlue / airDensityIN) / houseVolume * 3600;	// Flow from house air to outside
					else
						flueACH = (mFlue / airDensityOUT) / houseVolume * 3600;	// Flow from outside to inside


					//// ************** Brennan Less. Calculating the dry and moist air loads due to air exchange.  
						if(mCeiling >= 0) { // flow from attic to house. Brennan, create new variable. Include envelope, ceiling and air handler off flows, but NOT mSupReg or mRetReg
					//mHouseIN = mIN - mCeiling - mSupReg - mSupAHoff - mRetAHoff; //Above, all these things have been added into mIN or mOUT, depending on flow directions.
					mCeilingIN = mCeiling + mSupAHoff + mRetAHoff; //Do we want to include duct leakage in ventilation loads? If so, we need another term, including supply/reutrn temps. 
					//mHouseOUT = mOUT - mRetReg; //Why do we add them in above and then subtract them out here. I DO NOT understand. 
				} else {
					//mHouseIN = mIN - mSupReg;
					mCeilingIN = 0;
					//mHouseOUT = mOUT - mCeiling - mRetReg - mSupAHoff - mRetAHoff;
				} //Yihuan:add these calculation for mCeilingIN calculation 

					//Calculation of dry and moist air ventilation loads. Currently all values are positive, ie don't differentiate between heat gains and losses. Load doens't necessarily mean energy use...
					Dhda = 1.006 * abs(tempHouse - weather.dryBulb); // dry air enthalpy difference across non-ceiling envelope elements [kJ/kg].
					ceilingDhda = 1.006 * abs(tempHouse - tempAttic); //Dry air enthalpy difference across ceiling [kJ/kg]
					Dhma = abs((1.006 * (tempHouse - weather.dryBulb)) + ((HRHouse * (2501 + 1.86 * (tempHouse - C_TO_K)))-(weather.humidityRatio * (2501 + 1.86 * (weather.dryBulb - C_TO_K))))); //moist air enthalpy difference across non-ceiling envelope elements [kJ/kg]
					ceilingDhma = abs((1.006 * (tempHouse - tempAttic)) + ((HRHouse * (2501 + 1.86 * (tempHouse - C_TO_K))) - (HRAttic * (2501 + 1.86 * (tempAttic - C_TO_K))))); //moist air enthalpy difference across ceiling [kJ/kg]


					DAventLoad = (mHouseIN * Dhda) + (mCeilingIN * ceilingDhda); //[kJ/min], I prevoiusly used mHouseIN, but changed to mHouse. This is super simplified, in that it does not account for differing dT across ducts, ceiling, walls, etc. these assume balanced mass flow (mHouseIN=mHouseOUT), I had assumed a need to multiply by 60 to convert kg/s kg/m, but that was not necessary once results were available, so I removed it.
					MAventLoad = (mHouseIN * Dhma) + (mCeilingIN * ceilingDhma); //[kJ/min]. Brennan, mHouseIN needs to be a new variable that we create in the functions.cpp tab.

					TotalDAventLoad = TotalDAventLoad + DAventLoad; //sums the dry air loads
					TotalMAventLoad = TotalMAventLoad + MAventLoad; //sums the moist air loads

					// [END] Heat and Mass Transport ==================================================================================================================================

					// [START] IAQ Calculations =======================================================================================================================================
					
					//Calculate ventSum based on the sum of airflows and the AuxFanIndex value.
					if(AuxFanIndex == 0){ //do NOT count aux fans in relative exposure calculations
					
						if(ventSumIN > ventSumOUT)						//ventSum based on largest of inflow or outflow (includes flue flows as default)
							ventSum = ventSumIN + abs(flueACH);
						else
							ventSum = ventSumOUT + abs(flueACH);
							
					} else if(AuxFanIndex == 1){ //DO count aux fans in relative exposure calculations.
					
						ventSumIN = ventSumIN + nonRivecVentSumIN;
						ventSumOUT = ventSumOUT + nonRivecVentSumOUT;
						
						if(ventSumIN > ventSumOUT)						//ventSum based on largest of inflow or outflow (includes flue flows as default)
							ventSum = ventSumIN + abs(flueACH);
						else
							ventSum = ventSumOUT + abs(flueACH);
							
					}		

					if(ventSum <= 0)
						ventSum = .000001;
						
					sub_infiltrationModel(envC, envPressureExp, windSpeedMultiplier, shelterFactor, stackCoef, windCoef, weather.windSpeed, 	
						weather.dryBulb, ventSum, houseVolume, windSpeedCorrection, Q_wind, 
						Q_stack, Q_infiltration, Q_total, wInfil, InfCalc);
					
					relExpOld = relExp;
					relDoseOld = relDose;
					
					//relExp and relDose are calculated and summed for every minute of the year.
					
					relExp = sub_relativeExposure(Aeq, Q_total, relExpOld, dtau, houseVolume);
					relDose = relExp * (1 - exp(-rivecdt / 24)) + relDoseOld * exp(-rivecdt / 24);
					totalRelExp = totalRelExp + relExp;
					totalRelDose = totalRelDose + relDose;

					if(OccContType > 2 && occupied[weekend][hour] == 0){ //unoccupied period AND occupancy SVC
						totalRelExp = totalRelExp - relExp; //undo the addition.	
						relDose = relDoseOld; //revert to prior time step value of relDose. relExp keeps being calculated in real-time.
						totalRelDose = totalRelDose - relDose; //undo the addition
					}
					
					if(occupied[weekend][hour] == 1) { //house is occupied, then increment counter. 
						occupiedMinCount += 1; // counts number of minutes while house is occupied
					}
					

					indoorConc = sub_Pollutant(outdoorConc, indoorConc, indoorSource, houseVolume, qHouse, qDeposition, penetrationFactor, qAH, AHflag, filterEfficiency);

					
					//occupiedExp and occupiedDose are calculated and running summed only for occupied minutes of the year. Their values are otherwise remain fixed at the prior time step value. 
					//Occupancy-based controls use relExp and occupiedDose.
					
// 					if(occupied[weekend][hour] == 1) { //house is occupied, then calculate new occupiedDose and increment counter. 
// 					//Otherwise leave occupiedDose fixed at value from prior time step.
// 						occupiedMinCount += 1; // counts number of minutes while house is occupied
// 						occupiedDose = relDose;				// For annual average calculation, occupiedDose is unchanged when unoccupied.
// 						occupiedExp = relExp;
// 						totalOccupiedDose = totalOccupiedDose + occupiedDose; // Sums all dose values during occupied mins.
// 						totalOccupiedExp = totalOccupiedExp + occupiedExp;
// 					}

					/* -------dkm: To add flue in relDose/relExp calc the following two IF sentences have been added------
					if(flueACH <= 0)
					ventSumOUTD = ventSumOUT + abs(flueACH);
					else
					ventSumIND = ventSumIN + flueACH;

					if(ventSumIND > ventSumOUTD)	//ventSum based on largest of inflow or outflow
					ventsumD = ventSumIND;
					else
					ventsumD = ventSumOUTD;

					if(ventsumD <= 0)
					ventsumD = .000001;*/

					//Calculation of turnover
					// Automatically chosen based on previous selection of Aeq calculations (above) using AeqCalcs variable. Brennan. Need to change. 
// 					switch (AeqCalcs) {
// 					case 1:
// 						// 1. 62.2 ventilation rate with no infiltration credit
// 						turnover = (1 - exp(-(ventSum + wInfil) * rivecdt)) / (ventSum + wInfil) + turnoverOld * exp(-(ventSum + wInfil) * rivecdt); //This is used to compare Qtot (i.e., AEQ) to ventSum+wInfil
// 						//turnover = (1 - exp(-ventSum * rivecdt)) / ventSum + turnoverOld * exp(-ventSum * rivecdt); //This was the orignal formulation here. Assumes that ONLY fan flow is known for IAQ calculations. I had to change this based on use of AeqCalcs earlier in the code.
// 						break;
// 					case 2:
// 						// 2. 62.2 ventilation rate + infiltration credit from weather factors (w x NL)
// 						turnover = (1 - exp(-(ventSum + wInfil) * rivecdt)) / (ventSum + wInfil) + turnoverOld * exp(-(ventSum + wInfil) * rivecdt);
// 						break;
// 					case 3:
// 						// 3. 62.2 ventilation rate + 62.2 default infiltration credit (62.2, 4.1.3 p.4)
// 						turnover = (1 - exp(-(ventSum + defaultInfil) * rivecdt)) / (ventSum + defaultInfil) + turnoverOld * exp(-(ventSum + defaultInfil) * rivecdt);
// 						//turnover = (1 - exp(-(houseACH) * rivecdt)) / (houseACH) + turnoverOld * exp(-(houseACH) * rivecdt);
// 						break;
// 					case 4:
// 						// 1. 62.2 ventilation rate with no infiltration credit + existing home deficit. This is the same as the original formulation for AeqCalcs =1 (commented out above). Again, this is just deal with the use of AeqCalcs variable earlier in the codes.
// 						turnover = (1 - exp(-ventSum * rivecdt)) / ventSum + turnoverOld * exp(-ventSum * rivecdt);
// 						break;
// 					}

					// The 'real' turnover of the house using ACH of the house (that includes infiltration rather than just mechanical ventilation and flues)
					//turnoverReal = (1 - exp(-(houseACH) * rivecdt)) / (houseACH) + turnoverRealOld * exp(-(houseACH) * rivecdt);
					//relExp = Aeq * turnover; 
					//relExpReal = Aeq * turnoverReal;

					// Relative Dose and Exposure for times when the house is occupied
// 					if(occupied[weekend][hour] == 0) { //unoccupied. fixes dose at 1
// 						relDose = 1 * (1 - exp(-rivecdt / 24)) + relDoseOld * exp(-rivecdt / 24);
// 						relDoseReal = 1 * (1 - exp(-rivecdt / 24)) + relDoseRealOld * exp(-rivecdt / 24);
// 					} else { //occupied. dose calculated by based on turnover.
// 						relDose = relExp * (1 - exp(-rivecdt / 24)) + relDoseOld * exp(-rivecdt / 24);
// 						relDoseReal = relExpReal * (1 - exp(-rivecdt / 24)) + relDoseRealOld * exp(-rivecdt / 24);
// 					}


					
					
			
					//occupiedDose = relDose * occupied[weekend][hour];				// For annual average calculation, occupiedDose = 0 when unoccupied. 
					//occupiedExp = relExp * occupied[weekend][hour];					// For annual average calculation, occupiedExp = 0 when unoccupied. 
			
// 					occupiedDoseReal = relDoseReal * occupied[weekend][hour];		// dose and exposure for when the building is occupied
// 					occupiedExpReal = relExpReal * occupied[weekend][hour];			//If house is always occupied, then relExpReal = occupiedExpReal

// 					totalOccupiedDose = totalOccupiedDose + occupiedDose;   // Sums all dose and exposure values during the simulation (unoccupied periods are 0). 
// 					totalOccupiedExp = totalOccupiedExp + occupiedExp;

// 					totalOccupiedDoseReal = totalOccupiedDoseReal + occupiedDoseReal;   // Brennan. Added these for "real" calculations. total dose and exp over the occupied time period
// 					totalOccupiedExpReal = totalOccupiedExpReal + occupiedExpReal;
			
// 					meanRelDose = meanRelDose + relDose;
// 					meanRelExp = meanRelExp + relExp;

// 					meanRelDoseReal = meanRelDoseReal + relDoseReal; //Brennan, added these to include REAL calculations.
// 					meanRelExpReal = meanRelExpReal + relExpReal;

		




					// [END] IAQ calculations =======================================================================================================================================
					if(RHHouse >= 60){ //RHind60 and 70 count the minutes where high RH occurs
						RHind60 = 1;
						HumidityIndex = (RHHouse - 60) * (1.0 / (100 - 60));
					} else {
						RHind60 = 0;
						HumidityIndex = 0;
					}
					
					if(RHHouse >= 70){
						RHind70 = 1;
					} else {
						RHind70 = 0;
					}	
						

					RHtot60 = RHtot60 + RHind60; //These are summed for the year and summarized in the output file. 
					RHtot70 = RHtot70 + RHind70;
					HumidityIndex_Sum = HumidityIndex_Sum + HumidityIndex;
			

					// Writing results to a csv file
					double setpoint;

					if(hcFlag == 1)
						setpoint = heatThermostat[hour];
					else
						setpoint = coolThermostat[hour];

					// ================================= WRITING RCO DATA FILE =================================

					//File column names, for reference.
					//outputFile << "Time\tMin\twindSpeed\ttempOut\ttempHouse\tsetpoint\ttempAttic\ttempSupply\ttempReturn\tAHflag\tAHpower\tHcap\tcompressPower\tCcap\tmechVentPower\tHR\tSHR\tMcoil\thousePress\tQhouse\tACH\tACHflue\tventSum\tnonRivecVentSum\tfan1\tfan2\tfan3\tfan4\tfan5\tfan6\tfan7\trivecOn\tturnover\trelExpRIVEC\trelDoseRIVEC\toccupiedExpReal\toccupiedDoseReal\toccupied\toccupiedExp\toccupiedDose\tDAventLoad\tMAventLoad\tHROUT\tHRhouse\tRH%house\tRHind60\tRHind70" << endl; 

					// tab separated instead of commas- makes output files smaller
					if(printOutputFile) {
						outputFile << hour << "\t" << minuteYear << "\t" << weather.windSpeed << "\t" << weather.dryBulb << "\t" << tempHouse << "\t" << setpoint << "\t";
						outputFile << tempAttic << "\t" << tempSupply << "\t" << tempReturn << "\t" << AHflag << "\t" << AHfanPower << "\t";
						outputFile << compressorPower << "\t" << mechVentPower << "\t" << HR[3] * 1000 << "\t" << SHR << "\t" << Mcoil << "\t";
						outputFile << Pint << "\t"<< qHouse << "\t" << houseACH << "\t" << flueACH << "\t" << ventSum << "\t" << nonRivecVentSum << "\t";
						outputFile << fan[0].on << "\t" << fan[1].on << "\t" << fan[2].on << "\t" << fan[3].on << "\t" << fan[4].on << "\t" << fan[5].on << "\t" << fan[6].on << "\t";
						outputFile << rivecOn << "\t" << relExp << "\t" << relDose << "\t";
						outputFile << occupied[weekend][hour] << "\t"; 
						outputFile << weather.humidityRatio << "\t" << HR[0] << "\t" << HR[1] << "\t" << HR[2] << "\t" << HR[3] << "\t" << HR[4] << "\t" << RHHouse << "\t" << RHind60 << "\t" << RHind70 << "\t" << HumidityIndex << "\t" << dh.condensate << "\t" << indoorConc << endl;
						//outputFile << mHouse << "\t" << mHouseIN << "\t" << mHouseOUT << "\t" << mIN << "\t" << mOUT << "\t" << mCeiling << "\t" << mSupReg << "\t" << mSupAHoff << "\t" << mRetAHoff << "\t" << mRetReg << "\t" << mFanCycler << "\t" << mFlue << "\t" << mFloor << "\t" << mAH << endl;
						//outputFile << mHouse << "\t" << mHouseIN << "\t" << mHouseOUT << mCeiling << "\t" << mHouseIN << "\t" << mHouseOUT << "\t" << mSupReg << "\t" << mRetReg << "\t" << mSupAHoff << "\t" ;
						//outputFile << mRetAHoff << "\t" << mHouse << "\t"<< flag << "\t"<< AIM2 << "\t" << AEQaim2FlowDiff << "\t" << qFanFlowRatio << "\t" << C << endl; //Breann/Yihuan added these for troubleshooting
					}

					// ================================= WRITING MOISTURE DATA FILE =================================
					//File column names, for reference.
					//moistureFile << "HROUT\tHRattic\tHRreturn\tHRsupply\tHRhouse\tHRmaterials\tRH%house\tRHind60\tRHind70" << endl;
					if(printMoistureFile) {
						if(moistureModel == 1)
							for(int i=0; i<6; i++) {   // new humidity model wood nodes
								moistureFile << moisture_nodes.moistureContent[i] << "\t" << moisture_nodes.mTotal[i] << "\t";
								}
						moistureFile << weather.humidityRatio << "\t" << HRHouse << "\t" << HRAttic << "\t" << HRSupply << "\t" << HRReturn << "\t";
						moistureFile << RHHouse << "\t" << RHAttic << "\t" << tempHouse - C_TO_K << "\t" << tempAttic - C_TO_K << endl;
					}
			
					// ================================= WRITING Filter Loading DATA FILE =================================
			
					massFilter_cumulative = massFilter_cumulative + (mAH * 60);	// Time steps are every minute and mAH is in [kg/s]
					massAH_cumulative = massAH_cumulative + (mAH * 60);			// Time steps are every minute and mAH is in [kg/s]
			

					// Filter loading output file
					if(printFilterFile) {
						filterFile << massAH_cumulative << "\t"  << qAH << "\t" << AHfanPower << "\t" << retLF << endl;
					}
			
					// Calculating sums for electrical and gas energy use
					AH_kWh = AH_kWh + AHfanPower / 60000;						// Total air Handler energy for the simulation in kWh
					compressor_kWh = compressor_kWh + compressorPower / 60000;	// Total cooling/compressor energy for the simulation in kWh
					mechVent_kWh = mechVent_kWh + mechVentPower / 60000;		// Total mechanical ventilation energy for over the simulation in kWh
					gasTherm = gasTherm + hcap * 60 / 1000000 / 105.5;			// Total Heating energy for the simulation in therms
					furnace_kWh = gasTherm * 29.3;								// Total heating/furnace energy for the simulation in kWh
					dehumidifier_kWh += dh.power / 60000;
					//coolingLoad += capacityc / 60000;
					//latLoad += latcap / 60000;
					//heatingLoad += capacityh / 60000;

					// Average temperatures and airflows
					meanOutsideTemp = meanOutsideTemp + (weather.dryBulb - 273.15);		// Average external temperature over the simulation
					meanAtticTemp = meanAtticTemp + (tempAttic -273.15);		// Average attic temperatire over the simulation
					meanHouseTemp = meanHouseTemp + (tempHouse - 273.15);		// Average internal house temperature over the simulation
					meanHouseACH = meanHouseACH + houseACH;						// Average ACH for the house
					meanFlueACH = meanFlueACH + abs(flueACH);					// Average ACH for the flue (-ve and +ve flow considered useful)

					// Time keeping
					minuteYear++;													// Minute count of year

					if(minuteYear > (365 * 1440))
						break;
				}     // end of minute loop
				begin = end;
				switch (weatherFileType) {
				case 1:
					end = readTMY3(weatherFile);
					break;
				case 2:
					end = readEPW(weatherFile);
					break;
				}

				// Mold Index Calculations per ASHRAE 160, BDL 12/2016
				// Sheathing surfaces are Sensitive (1) and bulk wood is Very Sensitive (0)
				// Should these be using minute or avg hourly data?
				/*
				moldIndex_South = sub_moldIndex(1, moldIndex_South, b[3], moisture_nodes.PW[0], Time_decl_South); //South Roof Sheathing Surface Node						
				moldIndex_North = sub_moldIndex(1, moldIndex_North, b[1], moisture_nodes.PW[1], Time_decl_North); //North Roof Sheathing Surface Node
				moldIndex_BulkFraming = sub_moldIndex(0, moldIndex_BulkFraming, b[5], moisture_nodes.PW[2], Time_decl_Bulk); //Bulk Attic Framing Surface Node
				*/
			}        // end of hour loop
		}           // end of day loop
		//} while (weatherFile);			// Run until end of weather file

		//[END] Main Simulation Loop ==============================================================================================================================================

		// Close files
		weatherFile.close();
		if(printOutputFile)
			outputFile.close();
		if(printMoistureFile) 
			moistureFile.close();
		if(printFilterFile)
			filterFile.close();

		fanScheduleFile.close();

		double total_kWh = AH_kWh + furnace_kWh + compressor_kWh + mechVent_kWh + dehumidifier_kWh;

		meanOutsideTemp = meanOutsideTemp / minuteYear;
		meanAtticTemp = meanAtticTemp / minuteYear;
		meanHouseTemp = meanHouseTemp / minuteYear;
		meanHouseACH = meanHouseACH / minuteYear;
		meanFlueACH = meanFlueACH / minuteYear;
// 		meanRelDose = meanRelDose / minuteYear;
// 		meanRelExp = meanRelExp / minuteYear;
// 		meanRelDoseReal = meanRelDoseReal / minuteYear;				//Brennan. Added these and need to define in the definitions area. 
// 		meanRelExpReal = meanRelExpReal / minuteYear;

		meanRelExp = totalRelExp / minuteYear;
		meanRelDose = totalRelDose / minuteYear;
		if(OccContType > 2){
			meanRelExp = totalRelExp / occupiedMinCount;
			meanRelDose = totalRelDose / occupiedMinCount;
		}
		//meanOccupiedDose = totalOccupiedDose / occupiedMinCount; //Annual average relDose
		//meanOccupiedExp = totalOccupiedExp / occupiedMinCount; //Annual average relExp; occupiedMinCount = 525,600 for continuous occupancy.

// 		meanOccupiedDoseReal = totalOccupiedDoseReal / occupiedMinCount;
// 		meanOccupiedExpReal = totalOccupiedExpReal / occupiedMinCount;

		RHexcAnnual60 = RHtot60 / minuteYear;
		RHexcAnnual70 = RHtot70 / minuteYear;
		HumidityIndex_Avg = HumidityIndex_Sum / minuteYear;
		
		// Write summary output file (RC2 file)
		ofstream ou2File(summaryFileName); 
		if(!ou2File) { 
			cout << "Cannot open summary file: " << summaryFileName << endl;
			return 1; 
		}

			// ================================= WRITING ANNUAL SUMMAR (.RC2) DATA FILE =================================

		//Column header names
		ou2File << "Temp_out\tTemp_attic\tTemp_house";
		ou2File << "\tAH_kWh\tfurnace_kWh\tcompressor_kWh\tmechVent_kWh\ttotal_kWh\tmean_ACH\tflue_ACH";
		ou2File << "\tmeanRelExp\tmeanRelDose";
		ou2File << "\toccupiedMinCount\trivecMinutes\tNL\tenvC\tAeq\tfilterChanges\tMERV\tloadingRate\tDryAirVentLoad\tMoistAirVentLoad\tRHexcAnnual60\tRHexcAnnual70\tHumidityIndex_Avg\tdehumidifier_kWh" << endl;
		//Values
		ou2File << meanOutsideTemp << "\t" << meanAtticTemp << "\t" << meanHouseTemp << "\t";
		ou2File << AH_kWh << "\t" << furnace_kWh << "\t" << compressor_kWh << "\t" << mechVent_kWh << "\t" << total_kWh << "\t" << meanHouseACH << "\t" << meanFlueACH << "\t";
		ou2File << meanRelExp << "\t" << meanRelDose << "\t";
		ou2File << occupiedMinCount << "\t" << rivecMinutes << "\t" << NL << "\t" << envC << "\t" << Aeq << "\t" << filterChanges << "\t" << MERV << "\t" << loadingRate << "\t" << TotalDAventLoad << "\t" << TotalMAventLoad;
		ou2File << "\t" << RHexcAnnual60 << "\t" << RHexcAnnual70 << "\t" << HumidityIndex_Avg << "\t" << dehumidifier_kWh << endl;

		ou2File.close();

	}
	batchFile.close();

	//------------Simulation Start and End Times----------
	time(&endTime);
	string runEndTime = ctime(&endTime);

	cout << "\nStart of simulations\t= " << runStartTime;
	cout << "End of simulations\t= " << runEndTime << endl;
	//----------------------------------------------------

	return 0;
}
