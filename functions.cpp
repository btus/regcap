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

void f_flueFlow(double& tempHouse, double& flueShelterFactor, double& dPwind, double& dPtemp, double& h, double& Pint, int& numFlues, flue_struct* flue, double& mFlue,
	double& airDensityOUT, double& airDensityIN, double& dPflue, double& tempOut, double& Aeq, double& houseVolume, double& windPressureExp);

void f_floorFlow3(double& Cfloor, double& Cpfloor, double& dPwind, double& Pint,
	double& n, double& mFloor, double& airDensityOUT, double& airDensityIN, double& Hfloor, double& dPtemp);

void f_ceilingFlow(int& AHflag, double& Patticint, double& h, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& mCeiling, double& atticC, double& airDensityATTIC,
	double& airDensityIN, double& tempAttic, double& tempHouse, double& tempOut, double& airDensityOUT,
	double& mSupAHoff, double& mRetAHoff, double& supC, double& supn, double& retC, double& retn, double CCeiling);

void f_wallFlow3(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& wallCp,
	double& n, double& Cwall, double& h, double& Pint, double& dPtemp, double& dPwind,
	double& mWallIn, double& mWallOut, double& Hfloor);

void f_fanFlow(fan_struct& fan, double& airDensityOUT, double& airDensityIN);

void f_pipeFlow(double& airDensityOUT, double& airDensityIN, double& CP, double& dPwind, double& dPtemp,
	double& Pint, pipe_struct& Pipe, double& tempHouse, double& tempOut);

void f_winDoorFlow(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& h,
	double& wallCp, double& n, double& Pint, double& dPtemp, double& dPwind, winDoor_struct& winDoor);

void f_roofCpTheta(double* Cproof, int& windAngle, double* Cppitch, double& roofPitch);

double f_neutralLevel(double dPtemp, double dPwind, double Pint, double Cpr, double h);

void f_roofFlow(double& tempAttic, double& tempOut, double& airDensityATTIC, double& airDensityOUT, double& Cpr,
	double& atticPressureExp, double& Croof, double& roofPeakHeight, double& Patticint, double& dPtemp, double& dPwind,
	double& mRoofIn, double& mRoofOut, double& H);

void f_atticVentFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	atticVent_struct& atticVent, double& tempAttic, double& tempOut);

void f_soffitFlow(double& airDensityOUT, double& airDensityATTIC, double& CP, double& dPwind, double& dPtemp, double& Patticint,
	soffit_struct& soffit, double& soffitFraction, double& atticC, double& atticPressureExp, double& tempAttic, double& tempOut);

void f_atticFanFlow(fan_struct& atticFan, double& airDensityOUT, double& airDensityATTIC);

double heatTranCoef(double tempi, double tempa, double velocity);

double radTranCoef(double emissivity, double tempi, double tempj, double shapeFactor, double areaRatio);

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
	double& sc, 
	double* x,
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
) {
	
	int rhoSheathing;
	int rhoWood;
	int heatIterations;
	//int asize;
	//int asize2;
	int cpShingles;
	
	//double incsolar[4] = {0,0,0,0};
	vector<double> b(ATTIC_NODES+1,0);
	vector< vector<double> > A(ATTIC_NODES,b);
	//double A[ATTIC_NODES][ATTIC_NODES] = {0};
	double toldcur[ATTIC_NODES];
	double woodThickness;
	//double pws;
	double PW;
	double A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16;
	double denShingles;
	double M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M13, M14;
	double cp1, cp2, cp3, cp4, cp5, cp6, cp7, cp8, cp9, cp10, cp11, cp12, cp13, cp14, cp15, cp16;
	double kWood;
	double kAir;
	double muAir;
	double Rshingles;
	double Rval2, Rval3, Rval4, Rval5, Rval7, Rval8, Rval9, Rval10, Rval11, Rval14;
	double u;
	double H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H13, H14;
	double HI11, HI14;
	double F8t2, F8t4, F8t11, F8t14;
	double F2t4, F2t8, F2t11, F2t14;
	double F4t2, F4t8, F4t11, F4t14;
	double F11t2, F11t4;
	double F14t2, F14t4;
	double HR2t4, HR2t8, HR2t11, HR2t14;
	double HR4t2, HR4t8, HR4t11, HR4t14;
	double HR8t4, HR8t2;
	double HR11t2, HR11t4;
	double HR14t2, HR14t4;
	double EPS1, epsshingles;	
	double hr7;
	double HRG3, HRG5, HRS3, HRS5;
	double FRS, FG;
	double TGROUND;
	double alpha3, alpha5;

	// Node Identification (NOTE: The C++ array indices are one less than the node number e.g. node 1 = 0, node 2 = 1 etc.)
	// Node 1 is the Attic Air
	// Node 2 is the Inner North Sheathing
	// Node 3 is the Outer North Sheathing
	// Node 4 is the Inner South Sheathing
	// Node 5 is the Outer South Sheathing
	// Node 6 is all of the Wood (joists, trusses, etc.) lumped together
	// Node 7 is the Ceiling of the House
	// Node 8 is the Floor of the Attic
	// Node 9 is the Inner Gable Wall (both lumped together)
	// Node 10 is the Outer Gable Wall (both lumped together)
	// Node 11 is the Return Duct Outer Surface
	// Node 12 is the Return Duct Air
	// Node 13 is The Mass of the House
	// Node 14 is the Supply Duct Outer Surface
	// Node 15 is the Supply Duct Air
	// Node 16 is the House Air (all one zone)

	woodThickness = .015;		// thickness of sheathing material
	EPS1 = .9;         			// Emissivity of building materials
	switch(roofType) {
		case 1:			// asphalt shingles
			alpha5 = .92;
			alpha3 = .92;
			epsshingles = .91;
			Rshingles = .078;										// ASHRAE Fundamentals 2011 pg 26.7
			denShingles = 1100 * 2 * .005;					// asphalt shingles (factor of two because they overlap)
			cpShingles = 1260;									// CP asphalt shingles
			break;
		case 2: 			// red clay tile - edited for ConSol to be light brown concrete
			alpha5 = .58; 											// .67
			alpha3 = .58; 											// .67
			epsshingles = .9;
			Rshingles = .5;
			denShingles = 50;										// kg/m2
			cpShingles = 880;										// CP for tile roof
			break;
		case 3:			// low coating clay tile
			alpha5 = .5;
			alpha3 = .5;
			epsshingles = .9;
			Rshingles = .5;
			denShingles = 50;
			cpShingles = 880;										// CP for tile roof
			break;
		case 4:			// asphalt shingles  & white coating
			alpha5 = .15;
			alpha3 = .15;
			epsshingles = .91;
			Rshingles = .078;										// ASHRAE Fundamentals 2011 pg 26.7
			denShingles = 1100 * 2 * .005;					// asphalt shingles (factor of two because they overlap)
			cpShingles = 1260;									// CP asphalt shingles
			break;
	}

	PW = HROUT * pRef / (.621945 + HROUT);						// water vapor partial pressure pg 1.9 ASHRAE fundamentals 2009
	PW = PW / 1000 / 3.38;											// CONVERT TO INCHES OF HG
	TSKY = tempOut * pow((.55 + .33 * sqrt(PW)), .25);		// TSKY DEPENDS ON PW

	// Surface Area of Nodes
	A2 = planArea / 2 / cos(roofPitch * M_PI / 180);		// PITCHED SLOPE AREA
	A4 = A2;
	
	// the following are commented out for ConSOl becasue cement tile is flat and does not have increased surface area
	if(roofType == 2 || roofType == 3) {
        A3 = 1.5 * A2;	   // tile roof has more surface area for convection heat transfer
	} else {
		A3 = A2;																
	}

    A5 = A3;
	A6 = planArea * 1.5;							// Attic wood surface area
	A7 = planArea;									// Ceiling
	A8 = A7;											// Attic floor
	A9 = planArea / 2 * tan(roofPitch * M_PI / 180);	// Total endwall area
	A10 = A9;
	A11 = retArea;
	A12 = M_PI * retLength * retDiameter;
	A14 = supArea;
	A15 = M_PI * supLength * supDiameter;
	
	// surface area of inside of house minus the end walls, roof and ceiling
	// Currently assuming two stories with heights of 2.5m and 3.0m
	//A16 = 3 * pow(floorArea, .5) * 2 + 2.5 * pow(floorArea, .5) * 2 + 2 * floorArea;
	//A16 = numStories * storyHeight * pow(planArea, .5) * 2 + (2 * (numStories -1)) * planArea;
	A16 = 11 * pow(floorArea,0.5) + 2 * floorArea;	// Empirically derived relationship
	A13 = 6 * A16;									//  Surface area of everything in the house
	
	// Material densities
	rhoSheathing = 450;							// Sheathing
	rhoWood = 500;									// Wood

	// masses
	M1 = atticVolume * airDensityATTIC;					// mass of attic air
	M2 = .5 * A2 * rhoSheathing * woodThickness;		// 1/2 OF TOTAL
	M3 = M2 + denShingles * A2;							// OTHER 1/2 OUTSIDE SHEATHING, WOOD TCOND W/MMC
	M4 = M2;
	M5 = M3;
	M6 = 10 * planArea;										// Wild speculation
	M7 = .5 * 4 * planArea;									// MASS OF JOISTS DRYWALL AND INSULATION
	M8 = M7;
	M9 = .5 * rhoWood * woodThickness * A9;
	M10 = M9;
	M11 = retLength * M_PI * (retDiameter + retThickness) * retThickness * retrho;	// retArea * retrho * retThickness
	M12 = retVolume * airDensityRET;

	//  maybe not - 02/2004 need to increase house mass with furnishings and their area: say 5000kg furnishings
	M13 = (storyHeight * pow(floorArea, .5) * 4 * 2000 * .01 + planArea * .05 * 2000);		// mass of walls (5 cm effctive thickness) + mass of slab (aso 5 cm thick)
	M14 = supLength * M_PI * (supDiameter + supThickness) * supThickness * suprho;			// supArea * suprho * supThickness
	M15 = supVolume * airDensitySUP;
	M16 = houseVolume * airDensityIN;

	// Specific heat capacities
	cp1 = CpAir;
	cp2 = 1210;													// CP plywood
	cp4 = cp2;
	cp3 = cpShingles;
	cp5 = cp3;
	cp6 = 1630;													// CP wood
	cp7 = 1150;
	cp8 = cp7;
	cp9 = cp2;
	cp10 = cp9;
	cp11 = retCp;												// input
	cp12 = cp1;
	cp13 = 1300;												// combination of wood and drywall
	cp14 = supCp;
	cp15 = cp1;
	cp16 = cp1;

	// Thermal conductivities (k) [W/mK]and R-values [m2K/W] and the like
	kWood = 0.15;												// check with Iain about this
	//kAir = 0.02624;											// Thermal conductivity of air, now as function of air temperature
	kAir = 1.5207e-11 * pow(tempOld[15],3) - 4.8574e-8 * pow(tempOld[15],2) + 1.0184e-4 * tempOld[15] - 0.00039333;
	muAir = 0.000018462;										// Dynamic viscosity of air (mu) [kg/ms] Make temperature dependent  (this value at 300K)

	// changed to account for cathedralized attics
	if(roofRval == 0) {
		Rval2 = (woodThickness / kWood) + Rshingles;
	} else {
		Rval2 = roofRval + Rshingles;
	}

	Rval3 = Rval2;
	Rval4 = Rval2;
	Rval5 = Rval2;
	Rval7 = rceil;												// EFFECTIVE THERMAL RESISTANCE OF CEILING
	Rval8 = Rval7;

	Rval9 = 2.3;												// rvalue of insulated gable end walls
	//Rval9 = .5;												// rvalue of uninsulated gable end walls

	Rval10 = Rval9;
	Rval11 = retRval;
	Rval14 = supRval;

	/* most of the surfaces in the attic undergo both natural and forced convection
	the overall convection is determined by the forced and natural convection coefficients
	to the THIRD power, adding them, and taking the cubed root.  This is Iain's idea
	and it seemed to work for him*/

	// Characteristic velocity
	u = (matticenvin - matticenvout) / airDensityATTIC / AL4 / 4.0;
	if(u == 0) {
		u = abs(mCeiling) / 2 / airDensityATTIC / AL4 * 2 / 4.0;
		if(u == 0)
			u = .1;
	}

	// ITERATION OF TEMPERATURES WITHIN HEAT SUBROUTINE
	// THIS ITERATES BETWEEN ALL TEMPERATURES BEFORE RETURNING TO MAIN PROGRAM
	heatIterations = 0;
	for(int i=0; i < ATTIC_NODES; i++) {
		toldcur[i] = tempOld[i];
	}

	while(1) {

		heatIterations++;
		
		// reset array A to 0
		if(heatIterations > 1) {
			for(int i=0; i < ATTIC_NODES; i++) {
				for(int j=0; j < ATTIC_NODES+1; j++) {
					A[i][j] = 0;
				}
			}
		}

		// convection heat transfer coefficients
		// inner north sheathing
		H2 = heatTranCoef(tempOld[1], tempOld[0], u);

		// outer north sheathing
		H3 = heatTranCoef(tempOld[2], tempOut, windSpeed);

		// inner south sheathing
		H4 = heatTranCoef(tempOld[3], tempOld[0], u);

		// outer north sheathing
		H5 = heatTranCoef(tempOld[4], tempOut, windSpeed);

		// Wood (joists,truss,etc.)
		H6 = heatTranCoef(tempOld[5], tempOld[0], u);

		// Underside of Ceiling
		// modified to use fixed numbers from ASHRAE Fundamentals ch.3 on 05/18/2000
		if(AHflag != 0)
			H7 = 9;
		else
			H7 = 6;

		// House Mass
		// uses ceiling heat transfer coefficient as rest for house heat transfer coefficient	
		H13 = H7;

		// Attic Floor
		H8 = heatTranCoef(tempOld[7], tempOld[0], u);

		// Inner side of gable endwalls (lumped together)
		H9 = heatTranCoef(tempOld[8], tempOld[0], u);

		// Outer side of gable ends
		//tfilm10 = (tempOld[9] + tempOut) / 2;
		//H10 = (18.192 - .0378 * (tfilm10)) * pow(windSpeed, .8);
		H10 = heatTranCoef(tempOld[9], tempOut, windSpeed);

		if(ductLocation == 1) { //  Ducts in the house
			// Outer Surface of Ducts
			if(AHflag != 0) {
				H11 = 9;
				H14 = H11;
			} else {
				H11 = 6;
				H14 = H11;
			}
		} else {
			H11 = heatTranCoef(tempOld[10], tempOld[0], u);			// Outer Surface of Return Ducts
			H14 = heatTranCoef(tempOld[13], tempOld[0], u);			// Outer Surface of Supply Ducts
		}

		// Inner Surface of Ducts
		// from Holman   Nu(D) = 0.023*Re(D)^0.8*Pr(D)^0.4
		// Note Use of HI notation
		// I think that the following may be an imperical relationship
		// Return Ducts
		HI11 = .023 * kAir / retDiameter * pow((retDiameter * airDensityRET * abs(retVel) / muAir), .8) * pow((CpAir * muAir / kAir), .4);
		if(HI11 <= 0)
			HI11 = H11;

		// Supply Ducts
		HI14 = .023 * kAir / supDiameter * pow((supDiameter * airDensitySUP * supVel / muAir), .8) * pow((CpAir * muAir / kAir), .4);
		if(HI14 <= 0)
			HI14 = H14;

		// Radiation shape factors
		/* Only 5 nodes (2,4,8,11,14) are involved in radiation transfer in the attic
		The endwalls have a very small contribution to radiation exchange and are neglected.
		The wood may or may not contribute to radiation exchange, but their geometry is
		too complex to make any assumptions so it is excluded.
		Assumes that the duct is suspended above the floor, completely out of the insulation
		this will change in the future */

		if(ductLocation == 1) { //  Ducts in the house
			F8t2 = 1 / 2.0;
			F8t4 = F8t2;
			F2t8 = F8t2 * A8 / A2;
			F4t8 = F2t8;
			F2t4 = (1 - F2t8);
			F4t2 = (1 - F4t8);
			F4t11 = 0;
			F4t14 = 0;
			F2t11 = 0;
			F2t14 = 0;
			F14t2 = 0;
			F14t4 = 0;
			F11t4 = 0;
			F11t2 = 0;
		} else {			// ducts in the attic
			// 33.3% of each duct sees each sheathing surface (top third of duct)
			F14t2 = 1 / 2.0;
			F11t2 = F14t2;
			F14t4 = F14t2;
			F11t4 = F14t2;

			// Remaining 50% of each duct surface sees the floor
			// changed, the ducts don't see the floor
			// F11t8 = 0
			// F14t8 = F11t8

			// The ducts don't see each other
			// F11t14 = 0
			// F14t11 = 0

			F8t14 = 0;											// F14t8 * (A14 / 3) / A8
			F8t11 = 0;											// F11t8 * (A11 / 3) / A8
			F2t14 = F14t2 * (A14 / 3) / A2;
			F2t11 = F11t2 * (A11 / 3) / A2;
			F4t14 = F14t4 * (A14 / 3) / A4;
			F4t11 = F11t4 * (A11 / 3) / A4;
			F8t2 = 1 / 2.0;										// (1 - F8t14 - F8t11) / 2
			F8t4 = F8t2;
			F2t8 = F8t2 * A8 / A2;
			F4t8 = F2t8;
			F2t4 = (1 - F2t8 - F2t11 - F2t14);
			F4t2 = (1 - F4t8 - F4t11 - F4t14);

			// North Sheathing
			HR2t11 = radTranCoef(EPS1, tempOld[1], tempOld[10], F2t11, A2/(A11/3));
			HR2t14 = radTranCoef(EPS1, tempOld[1], tempOld[13], F2t14, A2/(A14/3));

			// South Sheathing
			HR4t11 = radTranCoef(EPS1, tempOld[3], tempOld[10], F4t11, A4/(A11/3));
			HR4t14 = radTranCoef(EPS1, tempOld[3], tempOld[13], F4t14, A4/(A14/3));

			// Return Ducts (note, No radiative exchange w/ supply ducts)
			HR11t4 = radTranCoef(EPS1, tempOld[10], tempOld[3], F11t4, A11/A4);
			HR11t2 = radTranCoef(EPS1, tempOld[10], tempOld[1], F11t2, A11/A2);

			// Supply Ducts (note, No radiative exchange w/ return ducts)
			HR14t4 = radTranCoef(EPS1, tempOld[13], tempOld[3], F14t4, A14/A4);
			HR14t2 = radTranCoef(EPS1, tempOld[13], tempOld[1], F14t2, A14/A2);
		}

		// North Sheathing
		HR2t4 = radTranCoef(EPS1, tempOld[1], tempOld[3], F2t4, A2/A4);
		HR2t8 = radTranCoef(EPS1, tempOld[1], tempOld[7], F2t8, A2/A8);

		// South Sheathing
		HR4t2 = radTranCoef(EPS1, tempOld[3], tempOld[1], F4t2, A4/A2);
		HR4t8 = radTranCoef(EPS1, tempOld[3], tempOld[7], F4t8, A4/A8);

		// Attic Floor
		HR8t4 = radTranCoef(EPS1, tempOld[7], tempOld[3], F8t4, A8/A4);
		HR8t2 = radTranCoef(EPS1, tempOld[7], tempOld[1], F8t2, A8/A2);

		// underside of ceiling
		hr7 = radTranCoef(EPS1, tempOld[6], tempOld[12], 1, A7/A13);

		FRS = (1 - sc) * (180 - roofPitch) / 180;      			// ROOF-SKY SHAPE FACTOR
		FG = 1 - FRS;                            					// ROOF-GROUND SHAPE FACTOR
		TGROUND = tempOut;                           			// ASSUMING GROUND AT AIR TEMP
		if(sc < 1) {
			HRS5 = radTranCoef(epsshingles, tempOld[4], TSKY, FRS, 0);
		} else {
			HRS5 = 0;
		}
		HRG5 = radTranCoef(epsshingles, tempOld[4], TGROUND, FG, 0);

		// South Sheathing
		if(sc < 1) {
			HRS3 = radTranCoef(epsshingles, tempOld[2], TSKY, FRS, 0);
		} else {
			HRS3 = 0;
		}
		HRG3 = radTranCoef(epsshingles, tempOld[2], TGROUND, FG, 0);

		// NODE 1 IS ATTIC AIR
		if(mCeiling >= 0) {
			// flow from attic to house
			A[0][0] = M1 * cp1 / dtau + H8 * A8 + H6 * A6 + mCeiling * cp1 + mSupAHoff * cp15 + mRetAHoff * cp12 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mRetLeak * cp1;
			b[0] = M1 * cp1 * tempOld[0] / dtau + matticenvin * cp1 * tempOut + mSupLeak * cp1 * toldcur[14];
		} else {
			// flow from house to attic
			A[0][0] = M1 * cp1 / dtau + H8 * A8 + H6 * A6 + H4 * A4 + H2 * A2 + A9 * H9 - matticenvout * cp1 - mRetLeak * cp1;
			b[0] = M1 * cp1 * tempOld[0] / dtau - mCeiling * cp1 * toldcur[15] - mSupAHoff * cp15 * toldcur[14] - mRetAHoff * cp12 * toldcur[11] + matticenvin * cp1 * tempOut + mSupLeak * cp15 * toldcur[14];
		}
		A[0][1] = -H2 * A2;
		A[0][3] = -H4 * A4;
		A[0][5] = -H6 * A6;
		A[0][7] = -H8 * A8;
		A[0][8] = -H9 * A9;
		if(ductLocation == 0) {		// duct surface conduction loss to attic
			A[0][0] +=  H14 * A14 / 2 + H11 * A11 / 2;
			A[0][10] = -H11 * A11 / 2;
			A[0][13] = -H14 * A14 / 2;
		}

		// NODE 2 IS INSIDE NORTH SHEATHING
		A[1][0] = -H2 * A2;
		A[1][1] = M2 * cp2 / dtau + H2 * A2 + A2 / Rval2 + HR2t4 * A2 + HR2t8 * A2;
		b[1] = M2 * cp2 * tempOld[1] / dtau;
		A[1][2] = -A2 / Rval2;
		A[1][3] = -HR2t4 * A2;
		A[1][7] = -HR2t8 * A2;

		if(ductLocation == 0) {			// duct surface radiation to sheathing
			A[1][1] += HR2t11 * A2 + HR2t14 * A2;
			A[1][10] = -HR2t11 * A2;
			A[1][13] = -HR2t14 * A2;
		}

		// NODE 3 IS OUTSIDE NORTH SHEATHING
		A[2][1] = -A2 / Rval3;
		A[2][2] = M3 * cp3 / dtau + H3 * A3 + A2 / Rval3 + HRS3 * A2 + HRG3 * A2;
		b[2] = M3 * cp3 * tempOld[2] / dtau + H3 * A3 * tempOut + A2 * nsolrad * alpha3 + HRS3 * A2 * TSKY + HRG3 * A2 * TGROUND;

		// NODE 4 IS INSIDE SOUTH SHEATHING
		A[3][0] = -H4 * A4;
		A[3][1] = -HR4t2 * A4;
		A[3][3] = M4 * cp4 / dtau + H4 * A4 + A4 / Rval4 + HR4t2 * A4 + HR4t8 * A4;
		b[3] = M4 * cp4 * tempOld[3] / dtau;
		A[3][4] = -A4 / Rval4;
		A[3][7] = -HR4t8 * A4;

		if(ductLocation == 0) {			// duct surface radiation to sheathing
			A[3][3] += HR4t11 * A4 + HR4t14 * A4;
			A[3][10] = -HR4t11 * A4;
			A[3][13] = -HR4t14 * A4;
		}

		// NODE 5 IS OUTSIDE SOUTH SHEATHING
		A[4][3] = -A4 / Rval5;
		A[4][4] = M5 * cp5 / dtau + H5 * A5 + A4 / Rval5 + HRS5 * A4 + HRG5 * A4;
		b[4] = M5 * cp5 * tempOld[4] / dtau + H5 * A5 * tempOut + A4 * ssolrad * alpha5 + HRS5 * A4 * TSKY + HRG5 * A4 * TGROUND;

		// NODE 6 IS MASS OF WOOD IN ATTIC I.E. JOISTS AND TRUSSES
		A[5][0] = -H6 * A6;
		A[5][5] = M6 * cp6 / dtau + H6 * A6;
		b[5] = M6 * cp6 * tempOld[5] / dtau;

		// NODE  7 ON INSIDE OF CEILING
		A[6][6] = M7 * cp7 / dtau + H7 * A7 + hr7 * A7 + A7 / Rval7;
		b[6] = M7 * cp7 / dtau * tempOld[6];
		A[6][7] = -A7 / Rval7;
		A[6][15] = -H7 * A7;
		A[6][12] = -hr7 * A7;

		//if(ductLocation == 1) {
			// ducts in house
		//	A[6][6] = M7 * cp7 / dtau + H7 * A7 + hr7 * A7 + A7 / Rval7;
		//}

		// NODE 8 ON ATTIC FLOOR
		A[7][0] = -H8 * A8;
		A[7][1] = -HR8t2 * A8;
		A[7][3] = -HR8t4 * A8;
		A[7][6] = -A8 / Rval8;
		A[7][7] = M8 * cp8 / dtau + H8 * A8 + HR8t2 * A8 + HR8t4 * A8 + A8 / Rval8;				// + HR8t11 * A8 + HR8t14 * A8
		b[7] = M8 * cp8 / dtau * tempOld[7];

		// NODE 9 IS INSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[8][0] = -H9 * A9;
		A[8][8] = M9 * cp9 / dtau + H9 * A9 + A9 / Rval9;
		A[8][9] = -A9 / Rval9;
		b[8] = M9 * cp9 * tempOld[8] / dtau;

		// NODE 10 IS OUTSIDE ENDWALLS THAT ARE BOTH LUMPED TOGETHER
		A[9][8] = -A10 / Rval10;
		A[9][9] = M10 * cp10 / dtau + H10 * A10 + A10 / Rval10;
		b[9] = M10 * cp10 * tempOld[9] / dtau + H10 * A10 * tempOut;

		// NODE 11 Exterior Return Duct Surface
		// Remember that the fluid properties are evaluated at a constant temperature
		// therefore, the convection on the inside of the ducts is
		b[10] = M11 * cp11 * tempOld[10] / dtau;
		A[10][11] = -A12 / (Rval11 + 1 / HI11);
		if(ductLocation == 1) {			// ducts in house
			A[10][10] = M11 * cp11 / dtau + H11 * A11 + A12 / (Rval11 + 1 / HI11);
			A[10][15] = -A11 * H11;
		} else {
			A[10][0] = -A11 * H11 / 2;
			A[10][1] = -A11 * HR11t2 / 3;
			A[10][3] = -A11 * HR11t4 / 3;
			A[10][10] += M11 * cp11 / dtau + H11 * A11 / 2 + A12 / (Rval11 + 1 / HI11) + A11 / 3 * HR11t2 + A11 / 3 * HR11t4;
		}
		
		// NODE 12 Air in return duct
		A[11][10] = -A12 / (Rval11 + 1 / HI11);
		if(mCeiling >= 0) {
			// flow from attic to house
			A[11][11] = M12 * cp12 / dtau + A12 / (Rval11 + 1 / HI11) + mAH * cp12 + mRetAHoff * cp12;
			b[11] = M12 * cp12 * tempOld[11] / dtau + mRetAHoff * cp1 * toldcur[0] - mRetLeak * cp1 * toldcur[0] - mRetReg * cp1 * toldcur[15] - mFanCycler * cp1 * tempOut - mHRV_AH * cp16 * ((1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15]) - mERV_AH * cp16 * ((1-ERV_SRE) * tempOut + ERV_SRE * tempOld[15]);
			
		} else {
			// flow from house to attic
			A[11][11] = M12 * cp12 / dtau + A12 / (Rval11 + 1 / HI11) + mAH * cp12 - mRetAHoff * cp12;
			b[11] = M12 * cp12 * tempOld[11] / dtau - mRetAHoff * cp16 * toldcur[15] - mRetLeak * cp1 * toldcur[0] - mRetReg * cp1 * toldcur[15] - mFanCycler * cp1 * tempOut - mHRV_AH * cp16 * ((1 - HRV_ASE) * tempOut + HRV_ASE * tempOld[15]) - mERV_AH * cp16 * ((1-ERV_SRE) * tempOut + ERV_SRE * tempOld[15]);
		}

		// node 13 is the mass of the structure of the house that interacts
		// with the house air to increase its effective thermal mass
		// 95% of solar gain goes to house mass, 5% to house air

		A[12][12] = M13 * cp13 / dtau + H13 * A13 + hr7 * A7;
		A[12][15] = -H13 * A13;
		A[12][6] = -hr7 * A7;
		b[12] = M13 * cp13 * tempOld[12] / dtau + .95 * solgain;

		// NODE 14
		b[13] = M14 * cp14 * tempOld[13] / dtau;
		A[13][14] = -A15 / (Rval14 + 1 / HI14);
		if(ductLocation == 1) {			// ducts in house
			A[13][13] = M14 * cp14 / dtau + H14 * A14 + A15 / (Rval14 + 1 / HI14);
			A[13][15] = -A14 * H14;
		} else {
			A[13][0] = -A14 * H14 / 2;
			A[13][1] = -A14 * HR14t2 / 3;
			A[13][3] = -A14 * HR14t4 / 3;
			A[13][13] = M14 * cp14 / dtau + H14 * A14 / 2 + A15 / (Rval14 + 1 / HI14) + A14 * HR14t2 / 3 + A14 / 3 * HR14t4;
		}

		// NODE 15 Air in SUPPLY duct
		// capacity is AC unit capcity in Watts
		// this is a sensible heat balance, the moisture is balanced in a separate routine.
		A[14][13] = -A15 / (Rval14 + 1 / HI14);
		if(mCeiling >= 0) {
			// flow from attic to house
			A[14][14] = M15 * cp15 / dtau + A15 / (Rval14 + 1 / HI14) + mSupReg * cp15 + mSupLeak * cp15 + mSupAHoff * cp15;
			b[14] = M15 * cp15 * tempOld[14] / dtau - capacityc + capacityh + evapcap + mAH * cp12 * toldcur[11] + mSupAHoff * cp1 * toldcur[0];
		} else {
			// flow from house to attic
			A[14][14] = M15 * cp15 / dtau + A15 / (Rval14 + 1 / HI14) + mSupReg * cp15 + mSupLeak * cp15 - mSupAHoff * cp15;
			b[14] = M15 * cp15 * tempOld[14] / dtau - capacityc + capacityh + evapcap + mAH * cp12 * toldcur[11] - mSupAHoff * cp16 * toldcur[15];
		}

		// NODE 16 AIR IN HOUSE
		// use solair tmeperature for house UA
		if(mCeiling >= 0) {
			// flow from attic to house
			A[15][15] = M16 * cp16 / dtau + H7 * A7 - mRetReg * cp16 - mHouseOUT * cp16 + H13 * A13 + uaSolAir + uaTOut;
			b[15] = M16 * cp16 * tempOld[15] / dtau + (mHouseIN - mHRV) * cp16 * tempOut + mHRV * cp16 * (( 1 - HRV_ASE) 
				* tempOut + HRV_ASE * tempOld[15]) + uaSolAir * tsolair + uaTOut * tempOut + .05 * solgain + mSupReg * cp1 * toldcur[14] 
				+ mCeiling * cp1 * toldcur[0] + mSupAHoff * cp15 * toldcur[14] + mRetAHoff * cp12 * toldcur[11] + internalGains + dhSensibleGain;
		} else {
			// flow from house to attic
			A[15][15] = M16 * cp16 / dtau + H7 * A7 - mCeiling * cp16 - mSupAHoff * cp16 - mRetAHoff * cp16 - mRetReg * cp16 - mHouseOUT * cp16 + H13 * A13 + uaSolAir + uaTOut;
			b[15] = M16 * cp16 * tempOld[15] / dtau + (mHouseIN - mHRV) * cp16 * tempOut + mHRV * cp16 * ((1 - HRV_ASE)
				* tempOut + HRV_ASE * tempOld[15]) + uaSolAir * tsolair + uaTOut * tempOut + .05 * solgain + mSupReg * cp1 * toldcur[14] + internalGains + dhSensibleGain;
		}

		A[15][6] = -H7 * A7;
		A[15][12] = -H13 * A13;

		if(ductLocation == 1) {
			// ducts in house
			A[15][15] += A11 * H11 + A14 * H14;
			A[15][10] = -A11 * H11;
			A[15][13] = -A14 * H14;
		}

		//asize = sizeof(A)/sizeof(A[0]);
		//asize2 = sizeof(A[0])/sizeof(A[0][0]);

		//ERRCODE = MatSEqn(A, b);
		for (int i=0; i<ATTIC_NODES; i++) {
			A[i][ATTIC_NODES] = b[i];
		}
		b = gauss(A);

		if(abs(b[0] - toldcur[0]) < .1) {
			break;
		} else {
			for(int i=0; i < ATTIC_NODES; i++) {
				toldcur[i] = b[i];
			}
		}
	} // END of DO LOOP
	for (int i=0; i<ATTIC_NODES; i++) {
		x[i] = b[i];
		}
}

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
	) {
		double Q;
		double R;
		double hrold[5];

		// 5 nodes
		// Node 1 is attic air (Node HR[0])
		// Node 2 is return air (Node HR[1])
		// Node 3 is supply air (Node HR[2])
		// Node 4 is house air (Node HR[3])
		// Node 5 is house materials that interact with house air only (Node HR[4])
		// this routine takes humidity ratios and mass flows and calculates W in each zone

		for(int i=0; i < 5; i++) {
			hrold[i] = HR[i];
		}

		// Node 1 - Attic Air
		if(mCeiling >= 0) {
			// flow from attic to house
			Q = Mw1 / dtau - matticenvout - mRetLeak + mCeiling + mSupAHoff + mRetAHoff;
			R = Mw1 * hrold[0] / dtau + matticenvin * HROUT + mSupLeak * hrold[2];
		} else {
			// flow from house to attic
			Q = Mw1 / dtau - matticenvout - mRetLeak;
			R = Mw1 * hrold[0] / dtau + matticenvin * HROUT + mSupLeak * hrold[2] - mCeiling * hrold[3] - mSupAHoff * hrold[2] - mRetAHoff * hrold[1];
		}
		HR[0] = R / Q;

		// Node 2 - Return Air
		if(mCeiling >= 0) { // flow from attic to house
			Q = Mw2 / dtau + mAH + mRetAHoff;
			R = Mw2 * hrold[1] / dtau - mRetLeak * hrold[0] + mRetAHoff * hrold[0] - mRetReg * hrold[3] - mFanCycler * HROUT - mERV_AH * ((1-ERV_TRE) * HROUT + ERV_TRE * hrold[3]) - mHRV_AH * HROUT; 
		} else { // flow from house to attic
			Q = Mw2 / dtau - mRetAHoff + mAH;
			R = Mw2 * hrold[1] / dtau - mRetAHoff * hrold[3] - mRetLeak * hrold[0] - mRetReg * hrold[3] - mFanCycler * HROUT - mERV_AH * ((1-ERV_TRE) * HROUT + ERV_TRE * hrold[3]) - mHRV_AH * HROUT;
		}
		HR[1] = R / Q;

		// Node 3 - Supply Air
		if(mCeiling >= 0) { // flow from attic to house
			Q = Mw3 / dtau + mSupLeak + mSupReg + mSupAHoff;
			R = Mw3 * hrold[2] / dtau + mAH * hrold[1] + mSupAHoff * hrold[0] - latcap / 2501000;
		} else { // flow from house to attic
			Q = Mw3 / dtau + mSupLeak + mSupReg - mSupAHoff;
			R = Mw3 * hrold[2] / dtau + mAH * hrold[1] - mSupAHoff * hrold[3] - latcap / 2501000;
		}
		HR[2] = R / Q;

		// Node 4 - House Air
		if(mCeiling >= 0) { // flow from attic to house
			Q = Mw4 / dtau - mHouseOUT - mRetReg + MWha;
			R = Mw4 * hrold[3] / dtau + mHouseIN * HROUT + mSupReg * hrold[2] + mCeiling * hrold[0] + mSupAHoff * hrold[2] + mRetAHoff * hrold[1] + latentLoad + MWha * hrold[4] - dhMoistureRemoved;
		} else { // flow from house to attic
			Q = Mw4 / dtau - mHouseOUT - mRetReg - mCeiling - mSupAHoff - mRetAHoff + MWha;
			R = Mw4 * hrold[3] / dtau + mHouseIN * HROUT + mSupReg * hrold[2] + latentLoad + MWha * hrold[4] - dhMoistureRemoved;
		}
		HR[3] = R / Q;

		// Node 5 - House furnishings/storage
		Q = Mw5 / dtau + MWha;
		R = Mw5 * hrold[4] / dtau + MWha * hrold[3];
		HR[4] = R / Q;
		
}

void sub_houseLeak (
	int& AHflag,
	int& flag, 
	double& windSpeed, 
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

		//if(flag < 1) {        // Yihuan: delete the if condition for the flag
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
				f_flueFlow(tempHouse, flueShelterFactor, dPwind, dPtemp, h, Pint, numFlues, flue, mFlue, airDensityOUT, airDensityIN, dPflue, tempOut, Aeq, houseVolume, windPressureExp);

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
				f_ceilingFlow(AHflag, Patticint, h, dPtemp, dPwind, Pint, envC, n, mCeiling, atticC, airDensityATTIC, airDensityIN, tempAttic, tempHouse, tempOut, airDensityOUT, mSupAHoff, mRetAHoff, supC, supn, retC, retn, CCeiling);

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
					
					f_wallFlow3(tempHouse, tempOut, airDensityIN, airDensityOUT, Cpwallvar, n, Cwall, h, Pint, dPtemp, dPwind, mWallIn, mWallOut, Hfloor);
					
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
					f_winDoorFlow(tempHouse, tempOut, airDensityIN, airDensityOUT, h, Cpwallvar, n, Pint, dPtemp, dPwind, winDoor[i]);
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
	int& flag, 
	double& windSpeed, 
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
	
	if(flag < 2) {
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
		f_roofFlow(tempAttic, tempOut, airDensityATTIC, airDensityOUT, Cpr, atticPressureExp, Croof, roofPeakHeight, Patticint, dPtemp, dPwind, mRoofIn, mRoofOut, h);

		mAtticIN = mAtticIN + mRoofIn;
		mAtticOUT = mAtticOUT + mRoofOut;

		// for second pitched part either back, above wall 2, or side above wall 4
		if(roofPeakPerpendicular) {
			Cpr = Cppitch[3] * Sw[3];
		} else {
			Cpr = Cppitch[1] * Sw[1];
		}

		// developed from wallflow3:
		f_roofFlow(tempAttic, tempOut, airDensityATTIC, airDensityOUT, Cpr, atticPressureExp, Croof, roofPeakHeight, Patticint, dPtemp, dPwind, mRoofIn, mRoofOut, h);

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

void f_flueFlow(double& tempHouse, double& flueShelterFactor, double& dPwind, double& dPtemp, double& h, double& Pint, int& numFlues, flue_struct* flue, double& mFlue,
	double& airDensityOUT, double& airDensityIN, double& dPflue, double& tempOut, double& Aeq, double& houseVolume, double& windPressureExp) {

		// calculates flow through the flue

		// dkm: internal variable to calculate external variable mFlue (because there may be more than one flue)
		double massflue = 0;
		double CpFlue;
		double fluePressureExp = 0.5;
		//double P = 0.14;	// Wind pressure coefficient of flue (0.14 for 

		for(int i=0; i < numFlues; i++) {
			//CpFlue = -.5 * pow(flueShelterFactor,2) * pow((flue[i].flueHeight / h),(2 * P));
			CpFlue = -.5 * pow((flue[i].flueHeight / h),(2 * windPressureExp));

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

void f_ceilingFlow(int& AHflag, double& Patticint, double& h, double& dPtemp,
	double& dPwind, double& Pint, double& C, double& n, double& mCeiling, double& atticC, double& airDensityATTIC,
	double& airDensityIN, double& tempAttic, double& tempHouse, double& tempOut, double& airDensityOUT,
	double& mSupAHoff, double& mRetAHoff, double& supC, double& supn, double& retC, double& retn, double CCeiling) {

		double dPceil = Pint - Patticint - airDensityOUT * g * ((tempHouse - tempOut) / tempHouse - (tempAttic - tempOut) / tempAttic) * h;

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
	double& n, double& Cwall, double& h, double& Pint, double& dPtemp, double& dPwind,
	double& mWallIn, double& mWallOut, double& Hfloor) {
		
		double Hwall = h - Hfloor;
		double dPwalltop = Pint + dPwind * wallCp - dPtemp * h;
		double dPwallbottom = Pint + dPwind * wallCp - dPtemp * Hfloor;
		double Bo = f_neutralLevel(dPtemp, dPwind, Pint, wallCp, h);

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

void f_winDoorFlow(double& tempHouse, double& tempOut, double& airDensityIN, double& airDensityOUT, double& h,
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
		double Bo = f_neutralLevel(dPtemp, dPwind, Pint, wallCp, h);

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
			if(Bo * h <= winDoor.Bottom) {
				Kwindow = .6;
				dummy = dummy2 * sqrt(abs(dummy2)) - dummy1 * sqrt(abs(dummy1));
				if(dT > 0) {
					winDoor.mIN = 0;
					winDoor.mOUT = sqrt(airDensityIN * airDensityOUT) * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
				} else {
					winDoor.mIN = -airDensityOUT * Kwindow * winDoor.Wide * tempHouse / 3 / g / dT * dummy;
					winDoor.mOUT = 0;
				}
			} else if(Bo * h > winDoor.Top) {
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
				if(Bo * h > Topcrit) {
					Pwindow = Topcrit * dPtemp - dPwind * wallCp;
				} else if(Bo * h < Bottomcrit) {
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

				if(Bo * h > Topcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * winDoor.High) * (h * Bo - Topcrit) + Kwindow;
				} else if(Bo * h < Bottomcrit) {
					Kwindow = (.6 - Kwindow) / (.1 * winDoor.High) * (Bottomcrit - h * Bo) + Kwindow;
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
	double& mRoofIn, double& mRoofOut, double& H) {
		
		double Hroof = roofPeakHeight - H;
		double dProoftop = Patticint + dPwind * Cpr - dPtemp * roofPeakHeight;
		double dProofbottom = Patticint + dPwind * Cpr - dPtemp * H;
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
				if(bRoof <= H / roofPeakHeight) {
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
				if(bRoof <= H / roofPeakHeight) {
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
* @param flow velocity - (m/s)
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
* @param shapeFactor - view factor between surfaces Fi-j
* @param areaRatio - ratio of surface areas (Ai/Aj)
* @return heat transfer coefficient (W/m2K)
*/

double radTranCoef(double emissivity, double tempi, double tempj, double shapeFactor, double areaRatio) {
	double rT = (1 - emissivity) / emissivity + 1 / shapeFactor + (1 - emissivity) / emissivity * areaRatio;
	return SIGMA * (tempi + tempj) * (pow(tempi, 2) + pow(tempj, 2)) / rT;
}
