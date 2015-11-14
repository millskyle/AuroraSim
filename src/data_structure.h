#ifndef DATASTRUCTURE
#define DATASTRUCTURE

#include <vector> 
#include <math.h>
#include "utility_functions.h" 


using namespace std;

typedef vector<float> threevector;

vector<float> randoms;

class Electron ; 

class Simulation {
   public: 
      float box_sizex = 100.0; //.push_back(100.00);
      float box_sizey = 100.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);
      int N = 100;
      int tmax = 10000;
      float E_low_threshold = 10; 
      float dt = 0.001;
      float velocity_scale = 10.0;
      float m_e = 1.0 ;
      int timescale_red_emission = 3000;
      int timescale_green_emission = 200;
      int timescale_blue_emission = 2;
      const int wavelength_red = 500.0;
      const int wavelength_green = 600.0;
      const int wavelength_blue = 800.0;
      
      const float hplanck=1.0;
      const float clight=1.0;
      vector<Electron> electrons;


};

class Electron {
public:
   float x;
   Simulation* sim;
   float y;
   float z;
   float vx;
   float vy;
   float vz;
   float E;
   bool emitting;
   float emitting_time_left;
   float emitting_wavelength;
   float Fx;
   float Fy;
   float Fz;
   float t;
   
   
   float get_p_interaction() {
      return 0.02;
   }




   float get_p_emit_red() {
      return 0.3333;
   }
   float get_p_emit_green() {
      return 0.3333;
   }
   float get_p_emit_blue() {
      return 0.3333;
   }

   //times that each emission should be active for  
   float get_t_emit_red() {
      vector<float> rnd = gen_random(2);
      return rnd[0] * sim->timescale_red_emission;
   }
   float get_t_emit_green() {
      vector<float> rnd = gen_random(2);
      return rnd[0] * sim->timescale_green_emission;
   }
   float get_t_emit_blue() {
      vector<float> rnd = gen_random(2);
      return rnd[0] * sim->timescale_blue_emission;
   }

   
   int reset(Simulation *sim) {
      t = 0;
      x = sim->box_sizex * (float)rand() / RAND_MAX;
      y = sim->box_sizey * (float)rand() / RAND_MAX;
      z = sim->box_sizez * 2;
      randoms = gen_random(3);
      vx = 0.001 * sim->velocity_scale * randoms[0];
      vy = 0.001 * sim->velocity_scale * randoms[1];
      vz = sim->velocity_scale * randoms[2];
      E = 0.5 * sim->m_e * (pow(vx,2) + pow(vy,2) + pow(vz, 2));
      emitting = 0;
      emitting_time_left = 0;
      emitting_wavelength = 0;
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
   }


};



/*
int tt[60][60];

struct Image {

   const int dim = 500;

   float R = new float[100][100][100];
   float G = new float[100][100][100];
   float B = new float[100][100][100];


   int bin_emission(Emission e, int boxx, int boxy, int boxz) {
      float x = e.position[0];
      float y = e.position[1];
      float z = e.position[2];

      if ( e.energy==1 ) {
         G[ int(x/boxx * dim) ][ int(y/boxy*dim) ][int(z/boxx * dim)] +=e.intensity ;
      } else if (e.energy == 2) {
         R[ int(x/boxx*dim) ][ int(y/boxy*dim) ][int(z/boxx * dim)]+= e.intensity ;
      } else if (e.energy == 3) {
         B[ int(x/boxx*dim) ][ int(y/boxy*dim) ][int(z/boxx * dim)]+=e.intensity ;
      }

   }





};
*/










#endif
