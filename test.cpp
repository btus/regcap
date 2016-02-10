// Code for testing C++ methods
// 2/4/16 LIR
#include <iostream>
#include <cmath>

using namespace std;

int main (int argc, char** argv) {
	
	string batchFile_name = "";
	if ( (argc <= 1) || (argv[argc-1] == NULL) || (argv[argc-1][0] == '-') ) {  // there is NO input...
		cerr << "usage: " << argv[0] << " batch_file" << endl;
      return(1);
   }
   else {
      batchFile_name = argv[argc-1];
   }
	
	cout << "batch file name: " << batchFile_name << endl;
  cout << "Here is PI:";
  cout << M_PI -3.14159 << endl;
}