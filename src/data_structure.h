#ifndef DATASTRUCTURE
#define DATASTRUCTURE

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>

using namespace std;

typedef vector<float> threevector;

vector<float> randoms;

class Electron ; 

class Simulation {
   public: 
      int N = 1000;
      int tmax = 1000;

      /* box size, in km */
      float box_sizex = 100.0; //.push_back(100.00);
      float box_sizey = 100.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);
      
      float dt=0.001;
   
      /* Mean energy of electrons when they begin, in keV */
      float E_mean = 10.0;
      
      float E_low_threshold = 10; 
      
      int timescale_red_emission = 50.; //110/dt;
      int timescale_green_emission = 50.;// 0.7/dt;
      int timescale_blue_emission = 50.;// 0.001/dt;


      float wavelength_red = 500.0;
      float wavelength_green = 600.0;
      float wavelength_blue = 800.0;
      float e_chg = 1e3; // elementary charge
      float epsilon_naught = 1e-9;
      
      vector<Electron> electrons;

      const float m_e = 1.0;
      float hc = 1.23984197 *10000; //*plank's constant in keV*nm 

      float number_density(float h) {
         return 7.86778379e6*exp(-h*0.176299767);
      }
      float bethe_energy_loss(float h,float v) {
         // dE/dx due to atmospheric 'drag'
         float effective_Z = 7.2; // 
         float Ionization = 92.0;
         float coeff = 4*M_PI*number_density(h)*pow(e_chg,4)*effective_Z / (m_e * pow(v,2));
         float parenthesis= log(  2*pow(v,2) / (m_e*Ionization  )) - 1.2329;
         return coeff*parenthesis;
      }

      bool in_sheet1(float x, float y) {
         int N = 4;
         float shift = 0.1*(float)rand()/RAND_MAX + 0.95;
         if (x < box_sizex*sin(shift * N * M_PI * (y - box_sizey/2) / box_sizey) 
             && x > box_sizex*sin(shift * N * M_PI * (y - box_sizey/2) / box_sizey) ) {
          return 1  ;
         } else {
            return 0;
         }

      }

 };




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
        p_emit_b = get_p_emit_blue(h);
        p_emit = p_emit_r + p_emit_g + p_emit_b ;
        p_emit_r = p_emit_r / p_emit;
        p_emit_g = p_emit_g / p_emit;
        p_emit_b = p_emit_b / p_emit;
        //cout << h<<"   "   << p_emit_r << "  " << p_emit_g << "  " << p_emit_b << "  "    << p_emit << endl;

   }
   float get_p_interaction(float h) {
      return 0.4*(get_p_emit_red(h) + get_p_emit_green(h) +  get_p_emit_blue(h)) ;
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

   
   int reset() {
      if (ID==3) {
         cout << "Respawning" << endl;
         cout << "  at respawn: Fx=" << Fx << "   Fy=" << Fy << "   Fz="<<Fz << endl;
         cout << "  at respawn:  x=" <<  x << "    y=" <<  y << "    z="<< z << endl;
      }
      t = 0;
      x = sim->box_sizex *(float)rand() / RAND_MAX;
      y = sim->box_sizey* (float)rand() / RAND_MAX;

       
     
      randoms = gen_random(6);
      
      x = (randoms[3])*75.0 + sim->box_sizex/2.0;
      y = abs(randoms[4])*75.0;

/*      if ((float)rand()/RAND_MAX > 0.9) {
         y = 75.0;
      } else {
         y = 30.0;
      }
      x = 50.0;

      */
      

      //
      //z = 0.1 * sim->box_sizez* (float)rand() / RAND_MAX + 0.9*sim->box_sizez;
      z = 1.1*sim->box_sizez -0.0001;

  //    z = 180.0 ; //sim->box_sizez*(float)rand()/RAND_MAX ;  

      E = 100000*sim->E_mean * (abs(randoms[0])+1.0);
      vx = (0.0001 * sqrt(2*E/sim->m_e) * randoms[1]);
      vy = (0.0001 * sqrt(2*E/sim->m_e) * randoms[2]);
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -sqrt( 2*E / sim->m_e - tmp ) ;
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
