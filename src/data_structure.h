#ifndef DATASTRUCTURE
#define DATASTRUCTURE

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
      float box_sizex = 100.0; //.push_back(100.00);
      float box_sizey = 100.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);
      int N = 10000;
      int tmax = 10000;
      float E_low_threshold = 10; 
      const float dt = 0.001;
      float velocity_scale = 100.0;
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
   int interaction_count;
   int iteration_count;
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
   float ID;
   
   
   float get_p_interaction() {
      return 0.02;
   }




   float get_p_emit_red() {
      return 0.1;
   }
   float get_p_emit_green() {
      return 0.2;
   }
   float get_p_emit_blue() {
      return 0.7;
   }

   //times that each emission should be active for  
   float get_t_emit_red() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]) * sim->timescale_red_emission;
   }
   float get_t_emit_green() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]) * sim->timescale_green_emission;
   }
   float get_t_emit_blue() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]) * sim->timescale_blue_emission;
   }

   
   int reset() {
      t = 0;
      x = sim->box_sizex * (float)rand() / RAND_MAX;
      y = sim->box_sizey * (float)rand() / RAND_MAX;
      z = sim->box_sizez -0.0001;
      randoms = gen_random(3);
      vx = 0.1 * sim->velocity_scale * randoms[0];
      vy = 0.1 * sim->velocity_scale * randoms[1];
      vz = - (sim->velocity_scale + ( sim->velocity_scale * 0.1 * randoms[2]) );
      E = 0.5 * sim->m_e * (pow(vx,2) + pow(vy,2) + pow(vz, 2));
      emitting = 0;
      emitting_time_left = 0;
      emitting_wavelength = 0;
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
   }


};


class PhotonDensity {


   public:
      unsigned int resolution_x = 250; //divisions along the box dimension
      unsigned int resolution_y = 250;
      unsigned int resolution_z = 500;

      float *R  = new float[resolution_x*resolution_y*resolution_z];
      float *G  = new float[resolution_x*resolution_y*resolution_z];
      float *B  = new float[resolution_x*resolution_y*resolution_z];

   
      int incr_element(float *A ,int i, int j, int k) {
         A[i + resolution_x *( j + resolution_y * k) ]++;
         return 0;
      }

      float get_element(float *A, int i, int j, int k) {
         return A[i + resolution_x *( j + resolution_y * k) ];
      }
       
      
     float get_max_voxel(float *A) {
        float maxx = 0;
        for (int i=0; i< (resolution_x*resolution_y*resolution_z)    ; i++) {
           maxx = max(maxx,A[i]);
        }
        return maxx;
     }

     float scalar_divide(float *A, float s) {
        for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++)  {
           A[i] = A[i] / s;
        }
        return *(A);
     }

     int destroy(float *A) {
        delete [] A;
        return 1;
     }
    
     float normalize() {
        float maxR,maxG,maxB;


        maxR = get_max_voxel(R);
        maxG = get_max_voxel(G);
        maxB = get_max_voxel(B);
        float gmax; //global max brightness
        gmax = max(maxR,maxG);
        gmax = max(gmax,maxB);

        cout << "\"Brightest\" voxel had un-normed brightness of " << gmax << endl;

        *R = scalar_divide(R,gmax); 
        *G = scalar_divide(G,gmax); 
        *B = scalar_divide(B,gmax); 


     }

     
     
     int write_image(int t) {
        cout << "Normalizing for max brightness" << endl;
        normalize();

        ofstream img;
        img.open("output/" + to_string(t) + ".dat");

        img << resolution_y << " " << resolution_z << "\n";

        float *Rflat  = new float[resolution_y*resolution_z];
        float *Gflat  = new float[resolution_y*resolution_z];
        float *Bflat  = new float[resolution_y*resolution_z];


        float pixelsumR =0;
        float pixelsumG =0;
        float pixelsumB =0;

        float maxpixelsumR =1;
        float maxpixelsumG =1;
        float maxpixelsumB =1;

        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              pixelsumR=0;
              pixelsumG=0;
              pixelsumB=0;
              for (int i=0; i<resolution_x; i++)  {
                 pixelsumR += get_element(R, i, j, k);
                 pixelsumG += get_element(G, i, j, k);
                 pixelsumB += get_element(B, i, j, k);
              }
              Rflat[k*resolution_y + j ] = pixelsumR;
              Gflat[k*resolution_y + j ] = pixelsumG;
              Bflat[k*resolution_y + j ] = pixelsumB;
              maxpixelsumR = max(maxpixelsumR,pixelsumR);
              maxpixelsumG = max(maxpixelsumG,pixelsumG);
              maxpixelsumB = max(maxpixelsumB,pixelsumB);
/*                if (j==100) {
                   Rflat[j*resolution_y + k] = 100;
                   Gflat[j*resolution_y + k] = 100;
                   Bflat[j*resolution_y + k] = 100;
                } else {
                   
                   Rflat[j*resolution_y + k] = 0;
                   Gflat[j*resolution_y + k] = 0;
                   Bflat[j*resolution_y + k] = 0;
                }
*/

           }
        }
        
        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              img << Rflat[k*resolution_y + j]/maxpixelsumR << " "
                  << Gflat[k*resolution_y + j]/maxpixelsumG << " " 
                  << Bflat[k*resolution_y + j]/maxpixelsumB << "\n";
           }
        }
        
        img.close();

          

          return 0;


     }
     
        


     




};




#endif
