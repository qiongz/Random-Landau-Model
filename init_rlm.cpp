#include"init_rlm.h"
using namespace std;

void usage(char* target) {
    cout << "Usage: " << target << " -[OPTIONS]" << endl;
    cout << "Options:\n";
    cout << "  -n  n_phi\n";
    cout << "  -q  quanta concentration \n";
    cout << "  -i  impurity concentration\n";
    cout << "  -g  theta mesh griding size (n_mesh)\n";
    cout << "  -d  disorder potential K-space cut (offhead)\n";
    cout << "  -c  # of disorder configurations\n";
    cout << "  -s  seed to generate seeds[n_sample]\n";
    cout << "  -t  # of threads\n";
}

int init(int argc, char* argv[], long &n_phi, dtype &quanta_concentration, dtype &impurity_concentration, long &n_mesh,
    long &n_sample,unsigned long &seed, long & num_threads) {
    // default values
    n_phi = 32;
    n_mesh = 32;
    n_sample = 1;
    quanta_concentration = 1;
    impurity_concentration = 16;
    num_threads=1;
    #if __cplusplus > 199711L
    seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #else
    Timer tmr;
    seed=tmr.nanoseconds();
    #endif

    /* Parameters initialization */
    void usage(char *);
    int getopt(int argc, char * const argv[], const char * optstring);
    extern char *optarg;
    int ch, err_flag;
    err_flag = 0;
    while((ch = getopt(argc, argv, "n:q:i:g:c:s:h:t:")) != -1) {
        switch(ch) {
        case 'n':
            n_phi = atoi(optarg);
            break;
        case 'q':
            quanta_concentration = atof(optarg);
            break;
        case 'i':
            impurity_concentration = atof(optarg);
            break;
        case 'g':
            n_mesh = atoi(optarg);
            break;
        case 'c':
            n_sample = atoi(optarg);
            break;
        case 's':
            seed = atoi(optarg);
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'h':
            err_flag++;
            break;
        default:
            err_flag++;
            break;
        }
    }
    if(err_flag) {
        usage(argv[0]);
        exit(2);
    }

    // write the params to log file
    ofstream flog;
    flog.open("log");
    flog<<"########################################################"<<endl;
    flog<<setw(25)<<"seed: "<<seed<<endl;
    flog<<setw(25)<<"n_phi: "<<n_phi<<endl;
    flog<<setw(25)<<"n_mesh: "<<n_mesh<<endl;
    flog<<setw(25)<<"n_sample: "<<n_sample<<endl;
    flog<<setw(25)<<"quanta_concentration: "<<quanta_concentration<<endl;
    flog<<setw(25)<<"impurity_concentration: "<<impurity_concentration<<endl;
    flog<<"########################################################"<<endl;
    flog.close();

    return 1;
}

