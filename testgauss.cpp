#include <iostream>
#include <cmath>
#include <vector>
#include "gauss.h"

using namespace std;

void print(vector< vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            cout << A[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            }
        }
        cout << "\n";
    }
    cout << endl;
}

int main() {
	 double Array[16][16] = {0};
	 double b[16] = {0};
    int n;
    cin >> n;

    vector<double> line(n+1,0);
    vector< vector<double> > A(n,line);

    // Read input data
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            cin >> A[i][j];
            Array[i][j] = A[i][j];
        }
    }

    for (int i=0; i<n; i++) {
        cin >> A[i][n];
        b[i] = A[i][n];
    }

    // Print input
    print(A);

    // Calculate solution
    vector<double> x(n);
    x = gauss(A);

	 MatSEqn(Array, b);

    // Print result
    cout << "Gauss Result:\t";
    for (int i=0; i<n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;
    cout << "MatSEqn Result:\t";
    for (int i=0; i<n; i++) {
        cout << b[i] << " ";
    }
    cout << endl;
}
