#ifndef SIMULATION
#define SIMULATION

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>

using namespace std;

typedef vector<float> threevector;

vector<float> randoms;

class Electron ; 

class Simulation {
   public: 
      int k = 20;  //number of particles will be N = 2^k 
      int tmax = 100;
      int t = 0;

      /* box size, in km */
      float box_sizex = 200.0; //.push_back(100.00);
      float box_sizey = 200.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);

      int E_field_recalc = 1;  //how many timesteps between E_field recalculations.

      float dt=0.001;
   
      /* Mean energy of electrons when they begin, in keV */
      float E_mean = 10.0;
      float E_loss_factor = 100.0;
      
      float E_low_threshold = 10; 
      
//      int timescale_red_emission = 1100.; //110/dt;
//      int timescale_green_emission = 50.;// 0.7/dt;
//      int timescale_blue_emission = 2.;// 0.001/dt;
      
      int timescale_red_emission = 500.; //110/dt;
      int timescale_green_emission = 50.;// 0.7/dt;
      int timescale_blue_emission = 5.;// 0.001/dt;


      float wavelength_red = 500.0;
      float wavelength_green = 600.0;
      float wavelength_blue = 800.0;
      float e_chg = 1e3; // elementary charge
      float epsilon_naught = 1e-9;
      
      vector<Electron> electrons;

      const float m_e = 1.0;
      float hc = 1.23984197 *10000; //*plank's constant in keV*nm 


      int N;

      void init() {
         N = pow(2,k);
//         for (int i=0; i<10000; i++) {
//            sheet_seeds.push_back((float)rand()/RAND_MAX) ; // = gen_random(100);
//         }

      }

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

 };





#endif
