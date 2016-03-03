#include <iostream>
using namespace std;

#include "config.h"

int main(int argc, char* argv[], char* envp[])
{
	cout << "REGCAP config test program" << endl;
	cout << endl;

	// read config file with environment variable expansion support
	Config config("regcap.cfg", envp);
	
	// File paths
	string inPath = config.pString("inPath");
	string outPath = config.pString("outPath");
	string weatherPath = config.pString("weatherPath");
	string shelterFile_name = config.pString("shelterFile_name");
	
	// output file control
	bool printMoistureFile = config.pBool("printMoistureFile");
	bool printFilterFile = config.pBool("printFilterFile");
	bool printOutputFile = config.pBool("printOutputFile");

	cout << "inPath = " << inPath << endl;
	cout << "outPath = " << outPath << endl;
	cout << "weatherPath = " << weatherPath << endl;
	cout << "shelterFile_name = " << shelterFile_name << endl;
	
	if(printMoistureFile)
		cout << "Moisture" << endl;
	if(printFilterFile)
		cout << "Filter" << endl;
	if(printOutputFile)
		cout << "Output" << endl;
		
	cout << endl;

	return 0;
}

