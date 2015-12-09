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
      int tmax = 10000;
      int t = 0;

      /* box size, in 'kilometers' */
      float box_sizex = 200.0; 
      float box_sizey = 200.0; 
      float box_sizez = 200.0; 

      int E_field_recalc = 1;  //how many timesteps between E_field recalculations.

      float dt=0.001; //timestep
   
      /* Mean energy of electrons when they begin, in keV */
      float E_mean = 8.0;
      float E_loss_factor = 100.0; //how much energy (factor) a photon emission costs.
      
      float E_low_threshold = 10; 
      
      //relative timescales for excited states
      int timescale_red_emission = 500.; 
      int timescale_green_emission = 100.;
      int timescale_blue_emission = 50.;

      //wavelength of the red emisison. Used to scale energy emission values (E=hc/lambda)
      float wavelength_red = 500.0;
      float wavelength_green = 600.0;
      float wavelength_blue = 800.0;


      //physical constants
      float e_chg = 1e3; // elementary charge
      float epsilon_naught = 1e-9;
      const float m_e = 1.0;
      float hc = 1.23984197 *10000; //*plank's constant * speed of light 
      
      vector<Electron> electrons;

      int N;

      void init() {
         N = pow(2,k);
      }

      float number_density(float h) {
         return 7.86778379e6*exp(-h*0.176299767);
      }


      float bethe_energy_loss(float h,float v) {
         // dE/dx due to atmospheric 'drag'
         //function that could be written to compute atmospheric drag, but I've incorporated this
         //pretty well into the probabilities
      }

 };





#endif
