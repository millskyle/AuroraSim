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
   vector<float> randoms;
   Simulation sim;


   for (int i=0; i<sim.N; i++) {
      Electron ee;
      ee.reset(&sim);
      ee.sim = &sim;
      sim.electrons.push_back(ee);
      cout << i <<endl;

   }


   float rnd;
   float p_interaction = 0.1;


   for (int t=0; t< sim.tmax; t++ ) {

      for (int i=0; i<sim.N ; i++) {
   cout << "here " << t<< ":"<<i<<"  " << endl;
         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
         rnd = (float)rand()/RAND_MAX;
         if (rnd < e->get_p_interaction() ) {
            cout << "emitting:";
            rnd = (float)rand()/RAND_MAX;
            if (rnd < e->get_p_emit_red()) {
               //emit red
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_red();
               e->emitting_wavelength = sim.wavelength_red;
            } else if (rnd < (e->get_p_emit_red() + e->get_p_emit_green()) ) {
               //emit green
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_green();
               e->emitting_wavelength = sim.wavelength_green;
            } else {
               //emit blue
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_blue();
               e->emitting_wavelength = sim.wavelength_blue;
               cout << "blue emission, ";
            }

         }

         

         
         
         
         
      }

   }
















return 0;
}
