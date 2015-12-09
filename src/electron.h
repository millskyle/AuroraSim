#ifndef ELECTRON
#define ELECTRON

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>

using namespace std;

typedef vector<float> threevector;

class Electron {
public:
   int interaction_count=0;
   int respawn_count=0;
   float x;
   Simulation* sim;
   float y;
   float z;
   float vx;
   float vy;
   float vz;
   float E=0;
   bool emitting;
   float emitting_time_left;
   float emitting_wavelength;
   int dead_counter;
   float Fx;
   float Fy;
   float Fz;
   float t;
   float ID;
   float tmp;
   float rnd; 

   float p_emit=0;
   float p_emit_r=0;
   float p_emit_g=0;
   float p_emit_b=0;

   //polynomial fits of emission rates
   float get_p_emit_red(float x) {  
      return 
         1.97692801419e-17 * pow(x, 8 ) + 
         -1.24675920355e-14 * pow(x, 7 ) + 
         2.61694799395e-12 * pow(x, 6 ) + 
        -1.20009942187e-10 * pow(x, 5 ) + 
        -3.03037545164e-08 * pow(x, 4 ) + 
        4.58332800758e-06 * pow(x, 3 ) + 
        -0.000224277684039 * pow(x, 2 ) + 
        0.00397540318554 * pow(x, 1 ) + 
       -0.00815824073091 * pow(x, 0 ) ;
   }

   float get_p_emit_green(float x) {
      return 1.22309846233e-16 * pow(x, 8 ) + 
         -1.04247753054e-13 * pow(x, 7 ) + 
         3.62402159771e-11 * pow(x, 6 ) + 
         -6.54860670058e-09 * pow(x, 5 ) + 
         6.39054458804e-07 * pow(x, 4 ) + 
         -2.94617910697e-05 * pow(x, 3 ) + 
         0.000107184265124 * pow(x, 2 ) + 
         0.0318371583695 * pow(x, 1 ) + 
         0.0028840983908 * pow(x, 0 ) ;
   }
   
   float get_p_emit_blue(float x) {
      return 
         7.71746456929e-17 * pow(x, 8 ) + 
        -6.74045102078e-14 * pow(x, 7 ) + 
        2.40765868528e-11 * pow(x, 6 ) + 
        -4.49364430269e-09 * pow(x, 5 ) + 
        4.59098779432e-07 * pow(x, 4 ) + 
        -2.34482985981e-05 * pow(x, 3 ) + 
        0.000317362998942 * pow(x, 2 ) + 
        0.0130756074812 * pow(x, 1 ) + 
        0.00402470459253 * pow(x, 0 );
  
   }

   float calculate_probabilities(float h) {
        p_emit_r = get_p_emit_red(h);
        p_emit_g = get_p_emit_green(h);
        p_emit_b = get_p_emit_blue(h)*0.5;
        p_emit = (p_emit_r + p_emit_g + p_emit_b)*exp((h - sim->box_sizez)/50) ;
        p_emit_r = p_emit_r / p_emit;
        p_emit_g = p_emit_g / p_emit;
        p_emit_b = p_emit_b / p_emit;
   }

   //times for which each emission should be active  
   float get_t_emit_red() {
      float rnd = gen_random();
      return (rnd+1) * sim->timescale_red_emission;
   }

   float get_t_emit_green() {
      float rnd = gen_random();
      return (rnd+1) * sim->timescale_green_emission;
   }

   float get_t_emit_blue() {
      float rnd = gen_random();
      return (rnd+1) * sim->timescale_blue_emission;
   }

   int random_collision() {
      //TO-DO: incorporate scattering here
      return 0;
    }

   int rescale_velocities() {
      //rescales the velocities to be consistent with the current kinetic energy
      float R = 0.5*sim->m_e*(vx*vx + vy*vy + vz*vz) / E;
      vz = vz/sqrt(R);
      return 0; 
   }

   
   int reset() {
     //sets up the electron/ resets it upon death.   
      t = 0;
      
      y = 0.6 * sim->box_sizey* (float)rand() / RAND_MAX + 0.2*sim->box_sizey;
      x = 0.6 * sim->box_sizex* (float)rand() / RAND_MAX + 0.2*sim->box_sizex;

      //if the time is past the equilibration time, reset to the top, else, place randomly
      if (sim->t > 100) {
         z = sim->box_sizez + 10*(float)rand()/RAND_MAX;
      } else {
         z = (float)rand()/RAND_MAX* sim->box_sizez ;
      }
     
      rnd = gen_random();

      E = (100000)*sim->E_mean * (abs(rnd)+1.0);
      vx = (0.001 * sqrt(2*E/sim->m_e) * (float)rand()/RAND_MAX );
      vy = (0.001 * sqrt(2*E/sim->m_e) * (float)rand()/RAND_MAX );
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -(sqrt( 2*E / sim->m_e - tmp )) ;
      emitting = 0;  //not currently emitting
      emitting_time_left = 0;  //not going to be emitting 
      emitting_wavelength = 0; //not emitting
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
      respawn_count +=1;
      dead_counter = 0;
   }


};

#endif
