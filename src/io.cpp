#ifndef IO
#define IO

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>

// This is IO as in INPUT/OUTPUT.  Named this prior to simulating Jupiter's moon Io, so this is a misleading
// filename.  Sorry.


using namespace std;


class EnergyDensity {
   //class to keep track of the energy deposited into the system. It's a 1D array.
   public:
      int resolution_z ;
      float *energy_density; 
      ofstream energy_density_out;
      
      int init(int resz) {
         resolution_z = resz;
         energy_density = new float[resz];
         for (int i=0; i< resolution_z; i++) {
            energy_density[i] = 0;
         }
      }

      //zero it out
      int reset() {
         for (int i=0; i< resolution_z; i++) {
            energy_density[i] = 0;
         }
         return 0;
      }

      int increment(int i, float n) {
         // increment if it's a valid bin, otherwise don't 
         // (prevents seq faults if electrons move outside box)
         if (i>=0 && i<resolution_z) {
            energy_density[i]+=n;
         } 
         return 0;
      }

      int write_out(int t) {
         energy_density_out.open("energy/E_density_" + to_string(t) + ".dat");
         for (int i=0; i<resolution_z; i++) {
            energy_density_out << energy_density[i] << "\t" << i << "\n";
         }
         energy_density_out.close();
      }

    void cleanup() {
    //   delete[] energy_density; //this was giving me problems.  TO-DO: FIX THIS
    }

};




class PhotonDensity {
   //class that keeps track of the photon density (i.e. the emissions).
   public:
      int resolution_x = 400; //divisions along the box dimension
      int resolution_y = 400; //increasing these values significantly increases 
      int resolution_z = 400; //the time required to output the image.  Too large and images are huge.
        
      ofstream img;

      float optical_decay_power = 2.0; //light intensity drops off as distance to this power (vacuum would be 2.0) 
      Simulation* sim; //pointer to Simulation object

      int init(Simulation *s) {
         sim = s; //set the simulation pointer
         return 0;
      }
      
      //we need a 3d array for each colour.
      float *R  = new float[resolution_x*resolution_y*resolution_z];
      float *G  = new float[resolution_x*resolution_y*resolution_z];
      float *B  = new float[resolution_x*resolution_y*resolution_z];

      int reset() {
         //function to zero out photon density
         cout << "\r Zeroing out photon density                                     " ;
         cout.flush();
         for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++ ) {                 
                 R[i] = 0;
                 G[i] = 0;
                 B[i] = 0;
         }

      }

      int incr_element(float *A ,int i, int j, int k,float n) {      
         //function to increment voxel. Handles 3D array indexing
         A[i + resolution_x *( j + resolution_y * k) ] += n;
         return 0;
      }

      float get_element(float *A, int i, int j, int k) {
         //function to access voxel, handles 3D array indexing
         return A[i + resolution_x *( j + resolution_y * k) ];
      }
       
      
     float get_max_voxel(float *A) {
        //function to return the value of the maximum voxel in an array
        float maxx = 0;
        for (int i=0; i< (resolution_x*resolution_y*resolution_z)    ; i++) {
           maxx = max(maxx,A[i]);
        }
        return maxx;
     }

     float scalar_divide(float *A, float s) {
        //divide an array by a scalar
        for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++)  {
           A[i] = A[i] / s;
        }
        return *(A);
     }

     int destroy(float *A) {
        //destroy the photon density array
        delete [] A;
        return 1;
     }
    
    
     int write_image(int t) {
        //function to write the image out to an output file. Writes a flat file of RGB values which I 
        //then use Python & matplotlib to render an image.
        img.open("output/" + to_string(t) + ".dat");
        
        img << resolution_y << " " << resolution_z << "\n";

        float *Rflat  = new float[resolution_y*resolution_z];
        float *Gflat  = new float[resolution_y*resolution_z];
        float *Bflat  = new float[resolution_y*resolution_z];
        
        float opt_decay = 0;

        float pixelsumR =0;
        float pixelsumG =0;
        float pixelsumB =0;
        int element = -1 ;

        float maxpixelsumR =0;
        float maxpixelsumG =0;
        float maxpixelsumB =0;
        float gmax = 0;

        //collapse the image along the x dimension
        for (int k=0; k<resolution_z; k++) {
           for (int i=0; i<resolution_x; i++) {
              pixelsumR=0;
              pixelsumG=0;
              pixelsumB=0;
              for (int j=0; j<resolution_y; j++)  {
                 opt_decay = pow( 
                     (  1.1 - (float(j))/float(resolution_y)  ) , optical_decay_power);
                 pixelsumR += get_element(R, i, j, k) / (opt_decay);
                 pixelsumG += get_element(G, i, j, k) / (opt_decay);
                 pixelsumB += get_element(B, i, j, k) / (opt_decay);
              }
              Rflat[k*resolution_y + i ] = pixelsumR;
              Gflat[k*resolution_y + i ] = pixelsumG;
              Bflat[k*resolution_y + i ] = pixelsumB;
           }
        }

        element = -1;
        //get a global max, which we normalize by (dependent somewhat arbitrarily on electron density)
        gmax = (float)sim->N * pow(((float)resolution_x / (float)sim->box_sizex),3)  * 0.03 ;
        
        //write out the (now) 2D data to a file
        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              element++;
              
              img << Rflat[k*resolution_y + j]/gmax << " "
                  << Gflat[k*resolution_y + j]/gmax << " " 
                  << Bflat[k*resolution_y + j]/gmax << "\n";
           }
        }
        
       img.close();
       
       //delete the flattened arrays, don't need them anymore
       delete[] Rflat;
       delete[] Gflat;
       delete[] Bflat;
       return 0;

     }
    

    //function to clean up pixel arrays
    void cleanup() {
       delete[] R;
       delete[] G;
       delete[] B;
    }
       


};

#endif
