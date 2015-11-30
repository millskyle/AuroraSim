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

    
   float get_p_emit_blue(float x) {
return 4.81871384251e-17 * pow(x, 8 ) + 
-5.14246386909e-14 * pow(x, 7 ) + 
2.14910606454e-11 * pow(x, 6 ) + 
-4.54662002061e-09 * pow(x, 5 ) + 
5.20401501824e-07 * pow(x, 4 ) + 
-3.18407383944e-05 * pow(x, 3 ) + 
0.000985761182095 * pow(x, 2 ) + 
-0.0133320014768 * pow(x, 1 ) + 
0.0612772920524 * pow(x, 0 ) ;
     
  }

  
   float get_p_emit_green(float x) {
return 3.81278117806e-16 * pow(x, 8 ) + 
-3.10001714635e-13 * pow(x, 7 ) + 
1.02888787226e-10 * pow(x, 6 ) + 
-1.78415951306e-08 * pow(x, 5 ) + 
1.71129354466e-06 * pow(x, 4 ) + 
-8.80251099314e-05 * pow(x, 3 ) + 
0.00219479356264 * pow(x, 2 ) + 
-0.0294122754696 * pow(x, 1 ) + 
0.666529089285 * pow(x, 0 );
}
   
   float get_p_emit_red(float x) {
      return 
 -1.11444871653e-15 * pow(x, 8 ) + 
8.60576875019e-13 * pow(x, 7 ) + 
-2.68171320654e-10 * pow(x, 6 ) + 
4.29637206481e-08 * pow(x, 5 ) + 
-3.71855716986e-06 * pow(x, 4 ) + 
0.000166259525058 * pow(x, 3 ) + 
-0.00332817581428 * pow(x, 2 ) + 
0.0249578808359 * pow(x, 1 ) + 
0.0359694822408 * pow(x, 0 ); 
   }


   float calculate_probabilities(float h) {
        p_emit_r = get_p_emit_red(h);
        p_emit_g = get_p_emit_green(h);
        p_emit_b = get_p_emit_blue(h);
        p_emit = p_emit_r + p_emit_g + p_emit_b ;
        p_emit_r = p_emit_r / p_emit;
        p_emit_g = p_emit_g / p_emit;
        p_emit_b = p_emit_b / p_emit;
        //cout << h<<"   "   << p_emit_r << "  " << p_emit_g << "  " << p_emit_b << "  "    << p_emit << endl;

   }

   //times that each emission should be active for  
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
      float rnd = gen_random();


      float theta = 2*M_PI*abs(rnd);
      vx = vx * cos(theta);
      vy = vy * sin(theta);
      vz = -sqrt( 2*E / sim->m_e - (vx*vx + vy*vy) );
    }

   int rescale_velocities() {
      //rescales the velocities to be consistent with the current kinetic energy
//      cout << "BEFORE: " << vx << "  " << vy << "  " << vz << endl;
      float R = 0.5*sim->m_e*(vx*vx + vy*vy + vz*vz) / E;
      //vx = vx/sqrt(R);
      //vy = vy/sqrt(R);
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
      
      
      y = 0.6 * sim->box_sizey* (float)rand() / RAND_MAX + 0.2*sim->box_sizey;
      x = 0.6 * sim->box_sizex* (float)rand() / RAND_MAX + 0.2*sim->box_sizex;

      if (sim->t > 100) {
         z = sim->box_sizez + 10*(float)rand()/RAND_MAX;
      } else {
         z = (float)rand()/RAND_MAX* sim->box_sizez ;
      }
     
//     if ((float)rand()/RAND_MAX > 0.98) {
//         y = 75.0 + (float)rand()/RAND_MAX;
//      } else {
//         y = 30.0 + (float)rand()/RAND_MAX;
//      }
 
      rnd = gen_random();

      E = (1000000)*sim->E_mean * (abs(rnd)+1.0);
      vx = (0.001 * sqrt(2*E/sim->m_e) * (float)rand()/RAND_MAX );
      vy = (0.001 * sqrt(2*E/sim->m_e) * (float)rand()/RAND_MAX );
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -(sqrt( 2*E / sim->m_e - tmp )) ;
      emitting = 0;
      emitting_time_left = 0;
      emitting_wavelength = 0;
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
      respawn_count +=1;
      dead_counter = 0;
   }


};

#endif
