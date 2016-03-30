// Code for testing C++ methods
// 2/4/16 LIR
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

std::istream& safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
}
int main (int argc, char** argv) {
	
	double A[16][16] = {0};
	string batchFileName = "";
	string test = "";
	if ( (argc <= 1) || (argv[argc-1] == NULL) || (argv[argc-1][0] == '-') ) {  // there is NO input...
		cerr << "usage: " << argv[0] << " batch_file" << endl;
      return(1);
   }
   else {
      batchFileName = argv[argc-1];
   }
	
	cout << "batch file name: " << batchFileName << endl;
	
	ifstream batchFile(batchFileName); 
	if(!batchFile) { 
		cout << "Cannot open: " << batchFileName << endl;
		return 1; 
	} 
	//safeGetline(batchFile,test);
	while(batchFile >> test)
	   cout << "\rTest =" << test << "= end of test";
	
	cout << endl;
	cout << "Here is PI:";
	cout << M_PI -3.14159 << endl;
	
	cout << "asize =" << sizeof(A)/sizeof(A[0]);
	cout << " asize2 =" << sizeof(A[0])/sizeof(A[0][0]) << endl;

}