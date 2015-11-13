#ifndef DATASTRUCTURE
#define DATASTRUCTURE

#include <vector> 
using namespace std;

typedef vector<float> threevector;


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


struct Image {

   const int dim = 500;

   float R[500][500] = {{0}};
   float G[500][500] = {{0}};
   float B[500][500] = {{0}};


   int bin_emission(Emission e, int boxx, int boxy, int boxz) {
      float x = e.position[0];
      float y = e.position[1];
      float z = e.position[2];

      if ( e.energy==1 ) {
         G[ int(x/boxx * dim) ][ int(y/boxy*dim) ]++;
      } else if (e.energy == 2) {
         R[ int(x/boxx*dim) ][ int(y/boxy*dim) ]++;
      } else if (e.energy == 3) {
         B[ int(x/boxx*dim) ][ int(y/boxy*dim) ]++;
      }

   }






};










#endif
