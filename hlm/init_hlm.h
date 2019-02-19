#include<unistd.h>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#if __cplusplus > 199711L
#include<chrono>
#endif
void usage(char* target);
int init(int argc, char* argv[], int &lx,int &ly, int & Q, int &band, int &n_mesh, int &n_sample, 
    unsigned long &seed, int & num_threads) ;

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
