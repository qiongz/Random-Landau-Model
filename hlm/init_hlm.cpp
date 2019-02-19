#include"init_hlm.h"
using namespace std;

void usage(char* target) {
    cout << "Usage: " << target << " -[OPTIONS]" << endl;
    cout << "Options:\n";
    cout << "  -x  lx\n";
    cout << "  -y  ly \n";
    cout << "  -q  q [p/q, p set to 1] \n";
    cout << "  -g  theta mesh griding size (n_mesh)\n";
    cout << "  -b  band to be projected onto\n";
    cout << "  -c  # of disorder configurations (n_sample)\n";
    cout << "  -s  seed to generate seeds\n";
    cout << "  -t  # of threads\n";
}

int init(int argc, char* argv[], int &lx,int &ly, int & Q, int &band, int &n_mesh, int &n_sample, 
    unsigned long &seed, int & num_threads) {
    // default values
    lx=3;
    ly=3;
    n_mesh = 32;
    n_sample = 1;
    Q = 3;
    band = 0; 
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
    while((ch = getopt(argc, argv, "x:y:b:q:g:c:s:h:t:")) != -1) {
        switch(ch) {
        case 'x':
            lx = atoi(optarg);
            break;
        case 'y':
            ly = atoi(optarg);
            break;
        case 'b':
            band = atoi(optarg);
            break;
        case 'q':
            Q = atoi(optarg);
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
    flog.open("log_hlm.txt");
    flog<<"########################################################"<<endl;
    flog<<setw(25)<<"seed: "<<seed<<endl;
    flog<<setw(25)<<"lx: "<<lx<<endl;
    flog<<setw(25)<<"ly: "<<ly<<endl;
    flog<<setw(25)<<"Q: "<<Q<<endl;
    flog<<setw(25)<<"n_mesh: "<<n_mesh<<endl;
    flog<<setw(25)<<"n_sample: "<<n_sample<<endl;
    flog<<setw(25)<<"n_thread: "<<num_threads<<endl;
    flog<<"########################################################"<<endl;
    flog.close();

    return 1;
}

