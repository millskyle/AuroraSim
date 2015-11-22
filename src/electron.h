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
   float Fx;
   float Fy;
   float Fz;
   float t;
   float ID;
   float tmp;
   vector<float>  rands;
   

   float p_emit=0;
   float p_emit_r=0;
   float p_emit_g=0;
   float p_emit_b=0;

    
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
        p_emit_b = get_p_emit_blue(h)/2.0;
        p_emit = p_emit_r + p_emit_g + p_emit_b ;
        p_emit_r = p_emit_r / p_emit;
        p_emit_g = p_emit_g / p_emit;
        p_emit_b = p_emit_b / p_emit;
        //cout << h<<"   "   << p_emit_r << "  " << p_emit_g << "  " << p_emit_b << "  "    << p_emit << endl;

   }
   float get_p_interaction(float h) {
      return 0.001 * exp(-6*h/sim->box_sizez)*(get_p_emit_red(h) + get_p_emit_green(h) +  get_p_emit_blue(h)) ;
   }

   //times that each emission should be active for  
   float get_t_emit_red() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_red_emission;
   }
   float get_t_emit_green() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_green_emission;
   }
   float get_t_emit_blue() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_blue_emission;
   }

   int random_collision() {
      vector<float> rands = gen_random(4);
      float theta = 2*M_PI*abs(rands[0]);
      vx = vx * cos(theta);
      vy = vy * sin(theta);
      vz = -sqrt( 2*E / sim->m_e - (vx*vx + vy*vy) );
    }

   int rescale_velocities() {
      //rescales the velocities to be consistent with the current kinetic energy
//      cout << "BEFORE: " << vx << "  " << vy << "  " << vz << endl;
      float R = 0.5*sim->m_e*(vx*vx + vy*vy + vz*vz) / E;
      vx = vx/sqrt(R);
      vy = vy/sqrt(R);
      vz = vz/sqrt(R);
//      cout << "BEFORE: " << vx << "  " << vy << "  " << vz << endl;
//      cout << " " << endl;
      

   }

   float A;
   float shift;
   float n;
   float hshift;

   int reset() {
      
      
      if (ID==3) {
         cout << "Respawning" << endl;
         cout << "  at respawn: Fx=" << Fx << "   Fy=" << Fy << "   Fz="<<Fz << endl;
         cout << "  at respawn:  x=" <<  x << "    y=" <<  y << "    z="<< z << endl;
      }
      t = 0;
      
      
      y = 0.333 * sim->box_sizey* (float)rand() / RAND_MAX + 0.333*sim->box_sizey;
      x = 0.333 * sim->box_sizex* (float)rand() / RAND_MAX + 0.333*sim->box_sizex;



/*      if ((float)rand()/RAND_MAX > 0.4) {
         A = sim->box_sizex/2.0 - 1.0;
         shift = sim->box_sizex/2.0;
         n = 6.0;
         hshift = sim->t/3000.0;
      } else {
         A = sim->box_sizex/2.0 - 1.0;
         shift = sim->box_sizex/2.0;
         n = 4.23455;
         hshift = -sim->t / 1000.00;
      }
      x = A * sin(n*M_PI*(y-hshift)/sim->box_sizex) + shift;
*/




//      x = 3.0;
//      y = 50.0;


      if (sim->t > 10000) {
         z = sim->box_sizez + 10*(float)rand()/RAND_MAX;
      } else {
         z = (float)rand()/RAND_MAX* sim->box_sizez ;
      }
      randoms = gen_random(3);
     
//     if ((float)rand()/RAND_MAX > 0.98) {
//         y = 75.0 + (float)rand()/RAND_MAX;
//      } else {
//         y = 30.0 + (float)rand()/RAND_MAX;
//      }
      
      E = 100000*sim->E_mean * (abs(randoms[0])+1.0);
      vx = 0.00; //(0.001 * sqrt(2*E/sim->m_e) * randoms[1]);
      vy = 0.00; //-(0.001 * sqrt(2*E/sim->m_e) * randoms[2]);
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -(sqrt( 2*E / sim->m_e - tmp )) ;
      emitting = 0;
      emitting_time_left = 0;
      emitting_wavelength = 0;
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
      respawn_count +=1;

   }


};

#endif
