#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */

#include "testing.cpp"
#include "utility_functions.h"
#include "data_structure.h"

using namespace std;




int main() {
srand ( time(NULL));

//random_emission_distribution();

vector<float> randoms;
vector<float> box_size;
box_size.push_back(100.00);
box_size.push_back(100.00);
box_size.push_back(200.00);

int N=100;
int dead=0;
float E_lost_interaction=1.0;
float E_low_threshold = 10;

float dt = 0.001;

float velocity_scale = 10.0;
float x[N];
float y[N];
float z[N];

float vx[N];
float vy[N];
float vz[N];

float E[N];

float m_e = 1.0;

for (int i = 0; i < N; i++) {
   x[i] = box_size[0] * (float)rand() / RAND_MAX;
   y[i] = box_size[1] * (float)rand() / RAND_MAX;
   z[i] = box_size[2] * 2;
   randoms = gen_random(3);
  
   vx[i] = 0.001 * velocity_scale * randoms[0];
   vy[i] = 0.001 * velocity_scale * randoms[1];
   vz[i] = velocity_scale * randoms[2];
  
   E[i] = 0.5 * m_e * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );
}

for (int t=0; t< 100; t++ ) { 
//   cout << "test" << endl;i
   for (int i=0; i< N ; i++) {

      E[i] -= E_lost_interaction;

      
      
      x[i] += vx[i] * dt; 
      y[i] += vy[i] * dt;
      z[i] += vz[i] * dt;
     
      E[i] = 0.5 * m_e * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );



      if (z[i] < 0 || E[i] < E_low_threshold) {
         cout << "DEATH! " << dead << " dead" << endl;
         x[i] = box_size[0] * (float)rand() / RAND_MAX;
         y[i] = box_size[1] * (float)rand() / RAND_MAX;
         z[i] = box_size[2] * 2;
         randoms = gen_random(3);
        
         vx[i] = 0.001 * velocity_scale * randoms[0];
         vy[i] = 0.001 * velocity_scale * randoms[1];
         vz[i] = velocity_scale * randoms[2];
        
         E[i] = 0.5 * m_e * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i] );
         
         dead+=1;
      }


   }




}











return 0;
}
