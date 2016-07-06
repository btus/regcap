#include "gauss.h"
#ifdef __APPLE__
   #include <cmath>        // needed for mac g++
#endif
#include <iostream>

using namespace std;

// ----- MatSEqn definitions -----

int MatSEqn(double A[][ArraySize], double* b) {
	// Error codes returned:
	//      0  no error                     -1  matrix not invertible
	//     -2  matrix not square            -3  inner dimensions different
	//     -4  matrix dimensions different  -5  result matrix dimensioned incorrectly
	//     any other codes returned are standard BASIC errors
	// 
	// -------------------------------------------------------------------
	
	int errcode = 0;
	//int bserrcode = 0;
	//int ERR = 0;
		
	int continuevar = 0;
	
	double x[ArraySize] = {0};
	int rpvt[ArraySize] = {0};
	int cpvt[ArraySize] = {0};
	
	//for(int i=0; i < ArraySize; i++) {
	//	x[i] = 0;
	//	rpvt[i] = 0;
	//	cpvt[i] = 0;
	//}
	
	errcode = matlu(A, rpvt, cpvt, continuevar);			// Get LU matrix

	if(continuevar != -1) {
		errcode = (errcode + 5) % 200 - 5;
		cout << "\nMatSeqn continue error: " << errcode << endl;
		return errcode;	
	}
	
	// check dimension of b
	//if(asize != bsize) {
	//	ERR = 197;
	//	errcode = (ERR + 5) % 200 - 5;
	//	cout << "\nMatSeqn size error: " << errcode << endl;
	//	return errcode;
	//}
	
	errcode = matbs(A, b, x, rpvt, cpvt);						// Backsolve system
	
	for(int i=0; i < ArraySize; i++) {
	   b[i] = x[i];														// Put solution in b for return
	}

	if(errcode !=0) {
		errcode = (errcode + 5) % 200 - 5;
		cout << "\nMatSeqn error: " << errcode << endl;
		return errcode;
	}
	
	//errcode = (ERR + 5) % 200 - 5;
	return errcode;	
}

int matlu(double A[][ArraySize], int* rpvt, int* cpvt, int& continuevar) {
	int errcode = 0;
	int tempswap;
	//int count;
	int r;
	int c;
	int bestrow;
	int bestcol;
	int rp;
	int cp;

	double rownorm[ArraySize];
	double max;
	double temp;
	double oldmax;

	// Checks if A is square, returns error code if not
	//if(asize != asize2) {
	//	errcode = 198;
	//	continuevar = 0;
	//	cout << "\nMatlu error: " << errcode << endl;
	//	return errcode;
	//}

	//count = 0;														// initialize count, continue
	continuevar = -1;
	
	for(int row = 0; row < ArraySize; row++) {							// initialize rpvt and cpvt
		rpvt[row] = row;
		cpvt[row] = row;
		rownorm[row] = 0;                							// find the row norms of A()

		for(int col = 0; col < ArraySize; col++) {
			rownorm[row] = rownorm[row] + abs(A[row][col]);
		}

		// if any rownorm is zero, the matrix is singular, set error, exit and do not continue
		if(rownorm[row] == 0) {
			errcode = 199;
			continuevar = 0;
			cout << "\nMatlu error: " << errcode << endl;
			return errcode;
		}
	}
	
	for(int pvt = 0; pvt < (ArraySize-1); pvt++) {
		// Find best available pivot
		// checks all values in rows and columns not already used for pivoting
		// and finds the number largest in absolute value relative to its row norm
		max = 0;
		
		for(int row = pvt; row < ArraySize; row++) {
			r = rpvt[row];
			for(int col = pvt; col < ArraySize; col++) {
				c = cpvt[col];
				temp = abs(A[r][c]) / rownorm[r];
				if(temp > max) {
					max = temp;
					bestrow = row;		    						// save the position of new max
					bestcol = col;
				}
			}
		}
		
		// if no nonzero number is found, A is singular, send back error, do not continue
		if(max == 0) {
			errcode = 199;
			continuevar = 0;
			cout << "\nMatlu error: " << errcode << endl;
			return errcode;
		} else if(pvt > 1 && max < (deps * oldmax)) {				// check if drop in pivots is too much
			errcode = 199;
		}
		
		oldmax = max;
		
		// if a row or column pivot is necessary, count it and permute rpvt or cpvt.
		// Note: the rows and columns are not actually switched, only the order in which they are used.		
		if(rpvt[pvt] != rpvt[bestrow]) {
			//count = count + 1;
			
			tempswap = rpvt[pvt];
			rpvt[pvt] = rpvt[bestrow];
			rpvt[bestrow] = tempswap;
		}    
		
		if(cpvt[pvt] != cpvt[bestcol]) {
			//count = count + 1;
			
			tempswap = cpvt[pvt];
			cpvt[pvt] = cpvt[bestcol];
			cpvt[bestcol] = tempswap;
		}
		
		//Eliminate all values below the pivot
		rp = rpvt[pvt];
		cp = cpvt[pvt];

		for(int row = (pvt+1); row < ArraySize; row++) {
			r = rpvt[row];
			A[r][cp] = -A[r][cp] / A[rp][cp];						// save multipliers
			
			for(int col = (pvt+1); col < ArraySize; col++) {
				c = cpvt[col];										// complete row operations
				A[r][c] = A[r][c] + A[r][cp] * A[rp][c];				
			}
		}
		
	}
	
	// if last pivot is zero or pivot drop is too large, A is singular, send back error
	if(A[rpvt[ArraySize-1]][cpvt[ArraySize-1]] == 0) {
		errcode = 199;
		continuevar = 0;
		cout << "\nMatlu error: " << errcode << endl;
		return errcode;
	} else if(((abs(A[rpvt[ArraySize-1]][cpvt[ArraySize-1]])) / rownorm[rpvt[ArraySize-1]]) < (deps * oldmax)) {
		// if pivot is not identically zero then continue remains TRUE
		errcode = 199;
	}
	
	if(errcode != 0 && errcode < 199) {
		continuevar = 0;
	}
	
	return errcode;
}

int matbs(double A[][ArraySize], double* b, double* x, int* rpvt, int* cpvt) {

	int c, r;

	// do row operations on b using the multipliers in L to find Lb
	for(int pvt = 0; pvt < (ArraySize-1); pvt++) {
		c = cpvt[pvt];		
		for(int row = pvt+1 ; row < ArraySize; row++) {
			r = rpvt[row];
			b[r] = b[r] + A[r][c] * b[rpvt[pvt]];
		}
	}

	// backsolve Ux=Lb to find x
	for(int row = ArraySize-1; row >= 0; row--) {
		c = cpvt[row];
		r = rpvt[row];
		x[c] = b[r];
		for(int col = (row+1); col < ArraySize; col++) {
			x[c] = x[c] - A[r][cpvt[col]] * x[cpvt[col]];		
		}
		x[c] = x[c] / A[r][c];
	}

	return 0;
}
