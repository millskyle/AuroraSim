#ifndef IO
#define IO

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>

using namespace std;



class PhotonDensity {


   public:
      int resolution_x = 600; //divisions along the box dimension
      int resolution_y = 600;
      int resolution_z = 600;

      float optical_decay_power = 1.0; //light intensity drops off as distance to this power (vacuum would be 2.0) 

      float *R  = new float[resolution_x*resolution_y*resolution_z];
      float *G  = new float[resolution_x*resolution_y*resolution_z];
      float *B  = new float[resolution_x*resolution_y*resolution_z];

      int reset() {
         cout << "\r Zeroing out photon density                                     " ;
         cout.flush();
         for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++ ) {                 
                 R[i] = 0;
                 G[i] = 0;
                 B[i] = 0;
         }

      }

   
      int incr_element(float *A ,int i, int j, int k,float n) {      
         A[i + resolution_x *( j + resolution_y * k) ] += n;
         /*if (i>1 && i<resolution_x-2) {
            A[(i+1) + resolution_x *( j + resolution_y * k) ] += n/4.0;
            A[(i-1) + resolution_x *( j + resolution_y * k) ] += n/4.0;
         }
         if (j>1 && j < resolution_y-2) {
            A[i + resolution_x *( (j+1) + resolution_y * k) ] += n/4.0;
            A[i + resolution_x *( (j-1) + resolution_y * k) ] += n/4.0;
         }
         if (k>1 && k < resolution_z-2) {
            A[i + resolution_x *( j + resolution_y * (k-1)) ] += n/4.0;
            A[i + resolution_x *( j + resolution_y * (k+1)) ] += n/4.0;
         }

      */




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
       // *R = scalar_divide(R,maxR);
        maxG = get_max_voxel(G);
//        *G = scalar_divide(R,maxR);
        maxB = get_max_voxel(B);
//        *B = scalar_divide(R,maxR);
        float gmax; //global max brightness
        gmax = max(maxR,maxG);
        gmax = max(gmax,maxB);

        cout << "\rBrightest voxel had un-normed brightness of " << gmax << "                             ";
        cout.flush();

        for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++)  {
           R[i] = R[i] / gmax;
           G[i] = G[i] / gmax;
           B[i] = B[i] / gmax;
        }
        
     }
     
     int write_image(int t) {
        cout << "\rNormalizing for max brightness                                 " ;
        cout.flush();
        normalize();

        ofstream img;
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

        float maxpixelsumR =1;
        float maxpixelsumG =1;
        float maxpixelsumB =1;
        float gmax = 0;

        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              pixelsumR=0;
              pixelsumG=0;
              pixelsumB=0;
              for (int i=0; i<resolution_x; i++)  {
                 opt_decay =  pow( 
                     (  1 - (float(i))/float(resolution_x)  ) , optical_decay_power);
                 pixelsumR += get_element(R, i, j, k) / (opt_decay);
                 pixelsumG += get_element(G, i, j, k) / (opt_decay);
                 pixelsumB += get_element(B, i, j, k) / (opt_decay);
              }
              Rflat[k*resolution_y + j ] = pixelsumR;
              Gflat[k*resolution_y + j ] = pixelsumG;
              Bflat[k*resolution_y + j ] = pixelsumB;
              maxpixelsumR = max(maxpixelsumR,pixelsumR);
              maxpixelsumG = max(maxpixelsumG,pixelsumG);
              maxpixelsumB = max(maxpixelsumB,pixelsumB);
           }
        }

        gmax = maxpixelsumB + maxpixelsumG + maxpixelsumR;
        element = -1;
        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              element++;
/*              if ( Rflat[element] > Gflat[element] ) {
                 Gflat[element] = 0;
              } else {
                 Rflat[element] =0; 
              }
              if (Rflat[element] > Bflat[element]) {
                 Bflat[element] = 0;
              } else {
                 Rflat[element] = 0;
              }
*/

              
//              img << Rflat[k*resolution_y + j]/maxpixelsumR << " "
//                  << Gflat[k*resolution_y + j]/maxpixelsumG << " " 
//                  << Bflat[k*resolution_y + j]/maxpixelsumB << "\n";

              img << Rflat[k*resolution_y + j]/gmax << " "
                  << Gflat[k*resolution_y + j]/gmax << " " 
                  << Bflat[k*resolution_y + j]/gmax << "\n";
           }
        }
        
        img.close();

          

          return 0;


     }
     

};




class OutputReport {
  public:
    vector<float> EvsT;
    int interactions;
    int respawns;

    ofstream initial_positions_file;
    int init() {
       initial_positions_file.open("output/initial_positions.txt");
    }


    int write(string filename) {
        ofstream report;
        report.open(filename);

       
       report << "Particle 3:\nTime-dependent energy" << endl;
       for (int t=0; t < EvsT.size()-1; t++ ) {
          report << t << "\t" << EvsT[t] << endl;
       }

       report << "\n" ;
       report << "Interactions per life:" << (float)interactions/(float)respawns << endl;
       report << "Interactions: " << (float)interactions << "       Respawns: "<<(float)respawns << endl;

       report.close();     
       return 0;

    }

};


#endif
