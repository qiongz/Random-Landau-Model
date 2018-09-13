#include<unistd.h>
#include<chrono>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#if defined DP
#define dtype double
#else
#define dtype float
#endif
void usage(char* target);
int init(int argc, char* argv[], long &n_phi, dtype &quanta_concentration, dtype &impurity_concentration, long &n_mesh,
	long &n_sample, unsigned long &seed, long &num_threads);

