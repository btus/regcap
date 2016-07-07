#pragma once
#ifndef gauss_h
#define gauss_h
#include <vector>

using namespace std;

// ----- MatSEqn forward declarations -----
/*
===================================================================
In theory, a matrix can be exactly singular.  Numerically, exact singularity is a rare occurrence because 
round off error turns exact zeroes into very small numbers. This data corruption makes it necessary to set
a criterion to determine how small a number has to be before we flag it as zero and call the matrix singular.

The following constants set the maximum drop allowed in weighted pivot values without error code -1 (matrix singular) 
being returned.  To increase the singularity sensitivity of the MatSEqn routine, increase the size of feps for 
floating precision calls or deps for double precision calls.  Similarly, decreasing the size of the constants will 
make the routines less sensitive to singularity.
===================================================================
*/

const float feps = .00001f;
const double deps = .00000000001;
const int ArraySize = 16;

int MatSEqn(double A[][ArraySize], double* b);
int matlu(double A[][ArraySize], int* rpvt, int* cpvt, int& continuevar);
int matbs(double A[][ArraySize], double* b, double* x, int* rpvt, int* cpvt);
vector<double> gauss(vector< vector<double> > A);

#endif