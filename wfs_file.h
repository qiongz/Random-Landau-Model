#include<fstream>
#include<cstring>
#include<iomanip>
#include<complex>
#if defined DP
#define dtype double
#else
#define dtype float
#endif
using namespace std;
void write_wfs(int theta_1, int theta_2, complex<dtype> * wfs, unsigned long size);
void read_wfs(int theta_1, int theta_2, complex<dtype> * wfs,unsigned long size);

