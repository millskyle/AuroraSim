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
      float strength = 1.0;

      vector<float> at(float t,float x,float y,float z){
         vector<float> r;
//         r.push_back((strength*sin(x-50)/(x-50)));
//         r.push_back((strength*sin(x-50)/(x-50)));
         r.push_back(0);
         r.push_back(-strength);
         r.push_back(0);
         return r;
      }
};


class ChargeDensity {
   public:
      int resolution_x =128 ;
      int resolution_y = 128 ;
      int resolution_z = 128 ;

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
 
      void cleanup() {
         delete[] p;
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
      float *Ex = NULL ;//= new float[Nx*Ny*Nz];
      float *Ey  = NULL ;//= new float[Nx*Ny*Nz];
      float *Ez = NULL  ;//= new float[Nx*Ny*Nz];
      fftw_complex *in = NULL  ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *out = NULL   ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *in2 = NULL   ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *out2 = NULL   ;// = new fftw_complex[Nx*Ny*Nz];

      fftw_complex *inPhiX = NULL   ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *outEX = NULL   ;// = new fftw_complex[Nx*Ny*Nz];

      fftw_complex *inPhiY = NULL   ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *outEY  = NULL  ;// = new fftw_complex[Nx*Ny*Nz];

      fftw_complex *inPhiZ = NULL   ;// = new fftw_complex[Nx*Ny*Nz];
      fftw_complex *outEZ  = NULL  ;// = new fftw_complex[Nx*Ny*Nz];
      float dx;
      float dy ;
      float dz ;

      float invNx ;
      float invNy ;
      float invNz ;

      fftw_plan forward;
      fftw_plan inverseX;
      fftw_plan inverseY;
      fftw_plan inverseZ;

      int init(Simulation *simm, ChargeDensity *rhoo) {
         sim = simm;
         rho = rhoo;
         Nx = rho->resolution_x ;
         Ny = rho->resolution_y ;
         Nz = rho->resolution_z ;
         Lx = sim->box_sizex ; 
         Ly = sim->box_sizey ;
         Lz = sim->box_sizez ;        
     
         dx = Lx/float(Nx-1) ;
         dy = Ly/float(Ny-1) ;
         dz = Lz/float(Nz-1) ;

         invNx = 1.0/(Nx*Nx*Nx);
         invNy = 1.0/(Ny*Ny*Ny);
         invNz = 1.0/(Nz*Nz*Nz);
         
         Ex = new float[Nx*Ny*Nz];
         Ey = new float[Nx*Ny*Nz];
         Ez = new float[Nx*Ny*Nz];
         in = new fftw_complex[Nx*Ny*Nz];
         out = new fftw_complex[Nx*Ny*Nz];
         in2 = new fftw_complex[Nx*Ny*Nz];
         out2 = new fftw_complex[Nx*Ny*Nz];
                          
         inPhiX = new fftw_complex[Nx*Ny*Nz];
         outEX = new fftw_complex[Nx*Ny*Nz];

         inPhiY = new fftw_complex[Nx*Ny*Nz];
         outEY = new fftw_complex[Nx*Ny*Nz];

         inPhiZ = new fftw_complex[Nx*Ny*Nz];
         outEZ = new fftw_complex[Nx*Ny*Nz];

         forward = fftw_plan_dft_3d(Nx,Ny,Nz,in, out, FFTW_FORWARD, FFTW_MEASURE );
         inverseX = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiX, outEX, FFTW_BACKWARD, FFTW_MEASURE );
         inverseY = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiY, outEY, FFTW_BACKWARD, FFTW_MEASURE );
         inverseZ = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiZ, outEZ, FFTW_BACKWARD, FFTW_MEASURE );
        


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

         ofstream debugfile;
         debugfile.open("fftw.plt");


         double elem[2] = {0};

                  int element = -1;

         float wavenumber_x = 0;
         float wavenumber_y = 0;
         float wavenumber_z = 0;

         
         ofstream rho_out;
         rho_out.open("rho/rho_" + to_string(sim->t) + ".dat");

         element = -1 ;
         for (int k=0; k<Nz; k++) {
            for (int j=0; j<Ny; j++) {
               for (int i=0; i< Nx; i++) {
                  element++;
                  // poisson's equation: del^2 phi = - rho / epsilon_naught
                  elem[1]=0;
                  elem[0] = (double)rho->p[element] / (sim->N*(Lx*Ly*Lz)/(Nx*Ny*Nz)     );
                  if (k==Nz/2.0 ) {
                     rho_out << i << "\t" << j << "\t" << elem[0] << endl;
                  }
                  in[element][0] = elem[0];
                  in[element][1] = elem[1];
               }
            }
         }

         fftw_execute(forward);

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
               if (2*k < Nz) {
                  wavenumber_z = (  (double)k*k) * invNz;
               } else { 
                  wavenumber_z = (  (double)(Nz - k)*(Nz - k) ) * invNz;
               }
       
               element++;    

               tmpReal = out[element][0] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
               tmpImag = out[element][1] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
               if (i==Nx/2 && j==Ny/2) {
                  debugfile << tmpReal*tmpReal << "\t" << tmpImag*tmpImag << endl;
               }
                
               tmpReal /= -sim->epsilon_naught;
               tmpImag /= -sim->epsilon_naught;

               inPhiX[element][1] =  -tmpReal * wavenumber_x ;
               inPhiX[element][0] =  tmpImag * wavenumber_x ;
               
               inPhiY[element][1] = -tmpReal * wavenumber_y ;
               inPhiY[element][0] = tmpImag * wavenumber_y ;
               
               inPhiZ[element][1] = -tmpReal * wavenumber_z ;
               inPhiZ[element][0] = tmpImag * wavenumber_z ;

            }
         }
      }

      fftw_execute(inverseX);
      fftw_execute(inverseY);
      fftw_execute(inverseZ);

      element=-1;
            for(int k=0; k< Nz; k++) {
         for (int j=0; j<Ny; j++) {
      for (int i=0; i<Nx; i++) {
                  element++;
                  // poisson's equation: del^2 phi = - rho / epsilon_naught
                  update_element(Ex, i, j, k, (float)outEX[element][0]/Nx) ;
                  update_element(Ey, i, j, k, (float)outEY[element][0]/Ny) ;
                  update_element(Ez, i, j, k, (float)outEZ[element][0]/Nz) ;

             }
          }
       }

//      fftw_destroy_plan(forward);
//      fftw_destroy_plan(inverseX);
//      fftw_destroy_plan(inverseY);
//      fftw_destroy_plan(inverseZ);
rho_out.close();
//if (sim->t==25) { exit(1); }
}

         
void cleanup() {
 
         delete[] Ex;
         delete[] Ey;
         delete[] Ez;

         delete[] in;
         delete[] out;
         delete[] in2;
         delete[] out2;

         delete[] inPhiX;
         delete[] outEX;
         
         delete[] inPhiZ;
         delete[] outEZ;
         
         delete[] inPhiY;
         delete[] outEY;
         

}







};




#endif

