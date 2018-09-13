#include"disorder_potential.h"

void generate_disorder_potential(float * impurity_x, float *impurity_y, float* impurity_intensity, int impurity_num, float L1, float L2, unsigned long seed, int sample_offset){
    VSLStreamStatePtr rndStream;
    vslNewStream(&rndStream, VSL_BRNG_MT19937,seed);
    vslSkipAheadStream(rndStream,impurity_num*3*sample_offset);
    unsigned int *i_buffer = new unsigned int[impurity_num];
    // initialize impurity position r=(x,y)
    vsRngUniform(METHOD_FLOAT, rndStream, impurity_num, impurity_x, 0.0, L1);
    vsRngUniform(METHOD_FLOAT, rndStream, impurity_num, impurity_y, 0.0, L2);
    // initialize impurity intensity with 1 and -1  randomly
    viRngUniformBits(METHOD_INT, rndStream, impurity_num, i_buffer);
    for(int i = 0 ; i < impurity_num; i ++)
        if(i_buffer[i] % 2 == 0)
            impurity_intensity[i] = 1;
        else
            impurity_intensity[i] = -1;

    vslDeleteStream(&rndStream);
    delete [] i_buffer;
}

