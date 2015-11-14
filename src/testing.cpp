#ifndef TESTING
#define TESTING

#include <random>
#include "data_structure.h" 
using namespace std;



/*int random_emission_distribution() {
   ofstream ff;
   srand (time(NULL));
   ff.open("out.xyz");
   int N = 1000;
   float spread = 5.0;
   float rnd1,rnd2;
   vector<float> origin;
   origin.push_back(10.0);
   origin.push_back(10.0);
   origin.push_back(10.0);
   ff << N << "\n";
   ff << "title \n";
   float x,y,z,r,theta,phi;
   float radius = 10.0;

   Emissions all_emissions;
  // Image outimage;
   
   for (int i=0; i<N; i++) {
      x = 10 * radius * (  (float)rand() / RAND_MAX -0.5) ;
      y = 10 * radius * (  (float)rand() / RAND_MAX -0.5) ;
      z = 10 * radius * ((float)rand() / RAND_MAX -0.5);

      if ( ( x*x + y*y + z*z )  < radius * radius ) {
         ff << "A " << x << "\t  " << y << "\t  " << z << "  \n";
         cout << "A \t" << x << "\t  " << y << "\t  " << z << "  \n";
      } else {
         i = i-1;
      }


      Emission em;
      em.timestep = 0;
      em.position.push_back(x);
      em.position.push_back(y);
      em.position.push_back(z);
      em.intensity = 1.0;
      em.energy = (float)rand() / RAND_MAX;
      em.altitude = 100.0;
      all_emissions.add_emission(em);

//      outimage.bin_emission(em, 100,100,100);

      

      


   }

ff.close();


}

*/

#endif


