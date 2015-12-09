#ifndef ELECTRO
#define ELECTRO

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>
#include <thread>

using namespace std;

class MagneticField {
   //Here we could define complicated functions to compute the magnetic field
   //but I'll just keep it simple for now and use a constant.  B = -B0 yhat 
   public:
      float strength = 1.0;
      vector<float> at(int t,float x,float y,float z){
         vector<float> r;
         r.push_back(0);
         r.push_back(-strength);
         r.push_back(0);
         return r;
      }
};


class ChargeDensity {
   public:
      //resolution of the mesh used for Fourier transforms.  Should be (but doesn't have to be) a power of 2
      int resolution_x = 128 ;
      int resolution_y = 128 ;
      int resolution_z = 128 ;

      //pointer to the array holding the charge density
      float *p = new float[resolution_x*resolution_y*resolution_z]; ;

      int init(Simulation *simm) {
         int k = simm->k; //2^k particles in the simulation
         int kk = floor(k/3) + 1; 
         resolution_x = int(pow(2,kk));
         resolution_y = int(pow(2,kk));
         resolution_z = int(pow(2,kk));
         p = new float[resolution_x*resolution_y*resolution_z];
      }
    
      
      int incr_element(int i, int j, int k) {
         //function to increment a bin.  Takes care of the 1D/3D array indexing
         p[i + resolution_x *( j + resolution_y * k) ] ++ ;
         return 0;
      }

      float get_element(int i, int j, int k) {
         //function to get the value of the charge density from a bin
         //takes care of the 3D indexing of the 1D bin
         return p[i + resolution_x *( j + resolution_y * k) ];
      }

      int reset() {
         //zero out the charge density
         for (int i=0; i<(resolution_x*resolution_y*resolution_z); i++) {
            p[i]=0;            
         }
      }
 
      void cleanup() {
         //deallocate memory
         delete[] p;
      }

};


class ElectricField {

   public:
      bool planned;
      int fftw_threads_init() ;
      int Nx; //real space mesh size
      int Ny;
      int Nz;
      float Lx; //real space dimensions
      float Ly;
      float Lz;
      Simulation *sim; //pointer to Simulation object
      ChargeDensity *rho; //pointer to charge density object
      float *Ex = NULL ; //Electric field x component
      float *Ey  = NULL ;// Electric field y component
      float *Ez = NULL  ; //electric field z component
      fftw_complex *in = NULL; //input 1D array (charge densities)
      fftw_complex *out = NULL; //fourier transform of rho
      fftw_complex *inPhiX = NULL;  //input to the reverse fourier transform
      fftw_complex *outEX = NULL; //output from the inverse fourier transform
      fftw_complex *inPhiY = NULL; // ""
      fftw_complex *outEY  = NULL; // ""
      fftw_complex *inPhiZ = NULL; // ""
      fftw_complex *outEZ  = NULL; // ""
      
      float dx; //fourier space dimensions
      float dy;
      float dz;

      float invNx ; //fourier space mesh size
      float invNy ;
      float invNz ;

      // FFT plans // 
      fftw_plan forward;
      fftw_plan inverseX;
      fftw_plan inverseY;
      fftw_plan inverseZ;
      
      ofstream rho_out;
      ofstream debugfile;

      int init(Simulation *simm, ChargeDensity *rhoo) {
         // init() takes a pointer to the Simulation and ChargeDensity objects
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
                          
         inPhiX = new fftw_complex[Nx*Ny*Nz];
         outEX = new fftw_complex[Nx*Ny*Nz];

         inPhiY = new fftw_complex[Nx*Ny*Nz];
         outEY = new fftw_complex[Nx*Ny*Nz];

         inPhiZ = new fftw_complex[Nx*Ny*Nz];
         outEZ = new fftw_complex[Nx*Ny*Nz];
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
         //computes the electric field by using fourier transforms on the 
         //charge density. Updates the Ex, Ey, Ez arrays  
         debugfile.open("fftw.plt");
         double elem[2] = {0};
         int element = -1;
         float wavenumber_x = 0;
         float wavenumber_y = 0;
         float wavenumber_z = 0;
         
         rho_out.open("rho/rho_" + to_string(sim->t) + ".dat");
      
         //plan the FFTs
         forward = fftw_plan_dft_3d(Nx,Ny,Nz,in, out, FFTW_FORWARD, FFTW_MEASURE );
         inverseX = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiX, outEX, FFTW_BACKWARD, FFTW_MEASURE );
         inverseY = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiY, outEY, FFTW_BACKWARD, FFTW_MEASURE );
         inverseZ = fftw_plan_dft_3d(Nx,Ny,Nz,inPhiZ, outEZ, FFTW_BACKWARD, FFTW_MEASURE );
       
         //build the fft input arrays
         element = -1 ;
         for (int k=0; k<Nz; k++) {
            for (int j=0; j<Ny; j++) {
               for (int i=0; i< Nx; i++) {
                  element++;
                  elem[1]=0; //no imaginary component of charge 
                  elem[0] = (double)rho->p[element] / (sim->N*(Lx*Ly*Lz)/(Nx*Ny*Nz));
                  if (k==Nz/2.0 ) {  //write the charge density to a file for plotting.
                     rho_out << i << "\t" << j << "\t" << elem[0] << endl;
                  }
                  in[element][0] = elem[0];
                  in[element][1] = elem[1];
               }
            }
         }

         //execute the forward transform
         fftw_execute(forward);

         double tmpReal;
         double tmpImag;

         // build up the reverse Fourier data. This is complicated...
         element = -1;
         for (int i=0; i<Nx; i++) {
            if (2*i < Nx) {  //FFTW returns transforms in a weird order, from N/2 to N, then 0 to N/2
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

                  //divide the real part of the fourier transform by k^2 = (kx**2 + ky**2 + kz**2)
                  tmpReal = out[element][0] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
                  //divide the imag part of the fourier transform by k^2 = (kx**2 + ky**2 + kz**2)
                  tmpImag = out[element][1] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
                  
                  //divide by epsilon_naught defined in Simulation object
                  tmpReal /= -sim->epsilon_naught;
                  tmpImag /= -sim->epsilon_naught;


                  // FLIP THE REAL AND IMAGINARY parts, then multiply by wavevectors
                  inPhiX[element][1] =  -tmpReal * wavenumber_x ;
                  inPhiX[element][0] =  tmpImag * wavenumber_x ;
                  
                  inPhiY[element][1] = -tmpReal * wavenumber_y ;
                  inPhiY[element][0] = tmpImag * wavenumber_y ;
                  
                  inPhiZ[element][1] = -tmpReal * wavenumber_z ;
                  inPhiZ[element][0] = tmpImag * wavenumber_z ;

               }
            }
         }
        
      //do the three reverse transforms in different threads
      thread first(  fftw_execute, inverseX );
      thread second(  fftw_execute, inverseY );
      thread third(  fftw_execute, inverseZ );
      //wait for the threads
      first.join();
      second.join();
      third.join();

      element=-1;
      for(int k=0; k< Nz; k++) {
         for (int j=0; j<Ny; j++) {
            for (int i=0; i<Nx; i++) {
                  element++;
                  //update the electric field, normalizing by N since FFTW does f = Nx * F(f(x))
                  update_element(Ex, i, j, k, (float)outEX[element][0]/Nx) ;
                  update_element(Ey, i, j, k, (float)outEY[element][0]/Ny) ;
                  update_element(Ez, i, j, k, (float)outEZ[element][0]/Nz) ;

             }
          }
       }

      rho_out.close();
   }

         
   void cleanup() {
      fftw_destroy_plan(forward);
      fftw_destroy_plan( inverseX);
      fftw_destroy_plan(inverseY);
      fftw_destroy_plan(inverseZ);

      delete[] Ex;
      delete[] Ey;
      delete[] Ez;

      delete[] in;
      delete[] out;

      delete[] inPhiX;
      delete[] outEX;
      
      delete[] inPhiZ;
      delete[] outEZ;
      
      delete[] inPhiY;
      delete[] outEY;

   }



};




#endif

