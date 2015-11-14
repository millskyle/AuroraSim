#ifndef UTILITYFUNCTIONS
#define UTILITYFUNCTIONS

#include <vector> 
#include <math.h>

using namespace std;



vector<float> gen_random(int N) {
   // """Return a vector of (N/2 + N%2) elements sampled from a normal distribution """
    vector<float> ret;
    float rnd1;
    float rnd2;
    for (int i=0; i<(N/2); i++ ) {
       rnd1 = (float)rand() / RAND_MAX;
       rnd2 = (float)rand() / RAND_MAX;
       ret.push_back( sqrt( -2.0 * log(rnd1)) * cos( 2 * M_PI * rnd2)   );
       ret.push_back( sqrt( -2.0 * log(rnd1)) * sin( 2 * M_PI * rnd2)   );
    }
    return ret;
}

#endif
