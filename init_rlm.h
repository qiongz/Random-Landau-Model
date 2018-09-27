#include<unistd.h>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#if __cplusplus > 199711L
#include<chrono>
#endif
#if defined DP
#define dtype double
#else
#define dtype float
#endif
void usage(char* target);
int init(int argc, char* argv[], long &n_phi, dtype &quanta_concentration, dtype &impurity_concentration, long &n_mesh,
	long &n_sample, unsigned long &seed, long &num_threads);

class Timer
{
public:
    Timer() {
        clock_gettime(CLOCK_REALTIME, &beg_);
    }

    double elapsed() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return (end_.tv_sec - beg_.tv_sec) +
               (end_.tv_nsec - beg_.tv_nsec)/1000000000.0;
    }

    unsigned long nanoseconds() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return (end_.tv_nsec - beg_.tv_nsec);
    }

    void reset() {
        clock_gettime(CLOCK_REALTIME, &beg_);
    }

private:
    timespec beg_, end_;
};
