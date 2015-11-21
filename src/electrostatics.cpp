#ifndef ELECTRO
#define ELECTRO

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>

using namespace std;

class MagneticField {
   public:
      float strength = 1000.0;

      vector<float> at(float t,float x,float y,float z){
         vector<float> r;
//         r.push_back((strength*sin(x-50)/(x-50)));
//         r.push_back((strength*sin(x-50)/(x-50)));
         r.push_back(0);
         r.push_back(strength);
         r.push_back(0);
         return r;
      }
};


class ChargeDensity {
   public:
      int resolution_x = 40;
      int resolution_y = 40;
      int resolution_z = 80;

      float *p = new float[resolution_x*resolution_y*resolution_z];
    
      int incr_element(int i, int j, int k) {
         p[i + resolution_x *( j + resolution_y * k) ] ++ ;
         return 0;
      }

      float get_element(int i, int j, int k) {
         return p[i + resolution_x *( j + resolution_y * k) ];
      }

      int reset() {
         for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++) {
            p[i]=0;            
         }
      }
 


};


class ElectricField {

   public:
      bool planned;
      int Nx;
      int Ny;
      int Nz;

      float Lx;
      float Ly;
      float Lz;
      Simulation *sim;
      ChargeDensity *rho;
      float *Ex ;//= new float[Nx*Ny*Nz];
      float *Ey ;//= new float[Nx*Ny*Nz];
      float *Ez ;//= new float[Nx*Ny*Nz];

      int init(Simulation *simm, ChargeDensity *rhoo) {
         sim = simm;
         rho = rhoo;
         Nx = rho->resolution_x ;
         Ny = rho->resolution_y ;
         Nz = rho->resolution_z ;
         Lx = sim->box_sizex ; 
         Ly = sim->box_sizey ;
         Lz = sim->box_sizez ;

          Ex = new float[Nx*Ny*Nz];
          Ey = new float[Nx*Ny*Nz];
          Ez = new float[Nx*Ny*Nz];
         planned = 0;

      }


      int update_element(float* A, int i, int j, int k, float val) {
         A[i + Nx *( j + Ny * k) ] = val ;
         return 0;
      }
      float get_element(float* A, int i, int j, int k) {
         if ( i >= 0 && i < Nx && j >= 0 && j < Ny && k>=0 && k < Nz ) {
            return A[i + Nx *( j + Ny * k) ] ;
         } else {
            return 0;
         }
      }
   
      int compute() {

         float dx = Lx/float(Nx-1) ;
         float dy = Ly/float(Ny-1) ;
         float dz = Lz/float(Nz-1) ;

         float invNx = 1.0/(Nx*Nx*Nx);
         float invNy = 1.0/(Ny*Ny*Ny);
         float invNz = 1.0/(Nz*Nz*Nz);

         fftw_complex *in = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *out = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *in2 = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *out2 = new fftw_complex[Nx*Ny*Nz];

         fftw_complex *inPhiX = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *outEX = new fftw_complex[Nx*Ny*Nz];

         fftw_complex *inPhiY = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *outEY = new fftw_complex[Nx*Ny*Nz];

         fftw_complex *inPhiZ = new fftw_complex[Nx*Ny*Nz];
         fftw_complex *outEZ = new fftw_complex[Nx*Ny*Nz];


         float *phi = new float[Nx*Ny*Nz];
         double elem[2] = {0};

         fftw_plan forward;
         fftw_plan inverse;
         fftw_plan inverseX;
         fftw_plan inverseY;
         fftw_plan inverseZ;

         int element = -1;

         float wavenumber_x = 0;
         float wavenumber_y = 0;
         float wavenumber_z = 0;

         if (!planned) {
            forward = fftw_plan_dft_3d(Nx,Ny,Nz,in, out, FFTW_FORWARD, FFTW_MEASURE );
            inverse = fftw_plan_dft_3d(Nx,Ny,Nz,in2, out2, FFTW_BACKWARD, FFTW_MEASURE );

            inverseX = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiX, outEX, FFTW_BACKWARD, FFTW_MEASURE );
            inverseY = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiY, outEY, FFTW_BACKWARD, FFTW_MEASURE );
            inverseZ = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiZ, outEZ, FFTW_BACKWARD, FFTW_MEASURE );
         }
        

         element = -1 ;
         for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
               for (int k=0; k< Nz; k++) {
                  element++;
                  // poisson's equation: del^2 phi = - rho / epsilon_naught
                  elem[1]=0;
                  elem[0] = - (double)(rho->p[element]);
//                  if (elem[0]!=0){ cout << elem[0] << " ";}
                  in[element][0] = elem[0];
                  in[element][1] = elem[1];
               }
            }
         }

         fftw_execute(forward);

/*         for (int i=0; i<out.size(); i++) {
            if (out[i]==0) {
            } else {
            cout << out[i] << endl;
            }
         }
*/
         double tmpReal;
         double tmpImag;

         element = -1;
         for (int i=0; i<Nx; i++) {
                   if (2*i < Nx) {
                     wavenumber_x = (  (double)i*i) * invNx;
                  } else { 
                     wavenumber_x = (  (double)(Nx - i)*(Nx - i) ) * invNx;
                  }

           for (int j=0; j<Ny; j++) {
                   if (2*j < Ny) {
                     wavenumber_y = (  (double)j*j) * invNy;
                  } else { 
                     wavenumber_y = (  (double)(Ny - j)*(Ny - j) ) * invNy;
                  }
               for (int k=0; k<Nz; k++) {
                  element++;
                 
                  if (2*k < Nz) {
                     wavenumber_z = (  (double)k*k) * invNz;
                  } else { 
                     wavenumber_z = (  (double)(Nz - k)*(Nz - k) ) * invNz;
                  }


                  //flip the real and imaginary parts
//                  in2[i + Nx *( j + Ny * k)][0] = -sqrt(i*i + j*j + k*k)*out[i + Nx *( j + Ny * k)][0] / ((double)i*i+j*j+k*k+1e-16);
//                  in2[i + Nx *( j + Ny * k)][1] = -sqrt(i*i + j*j + k*k) * out[i + Nx *( j + Ny * k)][1] / (i*i+j*j+k*k+1e-16);
                  //(  (double)i*i*dx*dx + (double)j*j*dy*dy + (double)k*k*dz*dz + 1e-20) / (2*M_PI*2*M_PI) ;
//                  cout << denom << "\t" ;
                  
                 
          tmpReal = out[element][0] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
          tmpImag = out[element][1] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);

          tmpReal /= -sim->epsilon_naught;
          tmpImag /= -sim->epsilon_naught;

          inPhiX[element][1] =  tmpReal * wavenumber_x ;
          inPhiX[element][0] =  -tmpImag * wavenumber_x ;
          
          inPhiY[element][1] = tmpReal * wavenumber_y ;
          inPhiY[element][0] =  -tmpImag * wavenumber_y ;
          
          inPhiZ[element][1] = tmpReal * wavenumber_z ;
          inPhiZ[element][0] = - tmpImag * wavenumber_z ;

          // in2[element][0] = - out[element][0] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
          // in2[element][1] = - out[element][1] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20) ;
               }
            }
         }

         //fftw_execute(inverse);
         
         fftw_execute(inverseX);
         fftw_execute(inverseY);
         fftw_execute(inverseZ);

         element=-1;
         for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
               for (int k=0; k< Nz; k++) {
                  element++;
                  // poisson's equation: del^2 phi = - rho / epsilon_naught
//                  cout << i << "|" << j << "|" << k << endl;
//                  cout <<"!"<< (float)outEX[element][0] << "!" <<endl;
                  update_element(Ex, i, j, k, (float)outEX[element][0]) ;
                  update_element(Ey, i, j, k, (float)outEY[element][0]) ;
                  update_element(Ez, i, j, k, (float)outEZ[element][0]) ;
                  //Ex[element] = (float)outEX[element][0] ;


//                  if (phi[element]==0) {
//                  }else{
//                     cout << phi[i + Nx *( j + Ny * k) ] << "  " ;
//                  }
               }
            }
         }

////////////////////////////////////////
/*  NO LONGER NEEDED IF I'M USING SPECTRAL METHODS TO FIND E DIRECTLY.



         float tmp0;
         float tmp1;

//cout << phi[0] << endl;
//cout << out[0] << endl;

         for (int i=1; i<Nx-1; i++) {
            for (int j=1; j<Ny-1; j++) {
               for (int k=1; k<Nz-1; k++) {

//         cout << "HETERERERAKSEJRHASHDKLA"<< endl<<endl;
  //                cout << "tmp0 = get_element(phi,"<< i-1 << "," << j << "," << k <<");" << endl;
                  tmp0 = get_element(phi, i-1, j, k);
                  tmp1 = get_element(phi, i+1, j, k); 
                  update_element(Ex, i, j, k,  (tmp0-tmp1)/(2*dx)      );
                  
                  tmp0 = get_element(phi, i, j-1, k);
                  tmp1 = get_element(phi, i, j+1, k); 
                  update_element(Ey, i, j, k,  (tmp0-tmp1)/(2*dy)      );
                  
                  tmp0 = get_element(phi, i, j, k-1);
                  tmp1 = get_element(phi, i, j, k+1); 
                  update_element(Ez, i, j, k,  (tmp0-tmp1)/(2*dz)      );
//                  if (i*i + j*j + k*k < 3000) {
//                     update_element(Ez, i, j, k,  10000.0    );
//                     update_element(Ey, i, j, k,  10000.0    );
//                     update_element(Ex, i, j, k,  1000.0    );
//                  }
               }
            }
         }

*/
////////////////////////////////////////

fftw_destroy_plan(forward);
fftw_destroy_plan(inverse);
fftw_destroy_plan(inverseX);
fftw_destroy_plan(inverseY);
fftw_destroy_plan(inverseZ);


}

         









};




#endif

