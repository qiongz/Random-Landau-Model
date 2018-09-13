#include<unistd.h>
#include<chrono>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
void usage(char* target);
int init(int argc, char* argv[], long &n_phi, float &quanta_concentration, float &impurity_concentration, long &n_mesh,
	long &n_sample, unsigned long &seed, long &num_threads);

