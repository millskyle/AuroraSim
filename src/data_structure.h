#ifndef DATASTRUCTURE
#define DATASTRUCTURE

#include <vector> 
#include <math.h>


using namespace std;

typedef vector<float> threevector;




struct GlobalInfo {


};


struct Emission {
   float timestep;
   float origin;
   vector<float> position;
   float intensity;
   int energy; //aurora only cares about three energies, so we'll use 0-green 1-red 2-blue
   float altitude;
};


struct Emissions {
   vector<Emission> Emissions;
   int add_emission(Emission e) {
      Emissions.push_back(e);
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
