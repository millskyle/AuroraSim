#ifndef DATASTRUCTURE
#define DATASTRUCTURE

#include <vector> 
#include <math.h>
#include "utility_functions.h" 
#include <fstream>
#include <string.h>
#include <fftw3.h>

using namespace std;

typedef vector<float> threevector;

vector<float> randoms;

class Electron ; 

class Simulation {
   public: 
      int N = 1000;
      int tmax = 1000;

      /* box size, in km */
      float box_sizex = 100.0; //.push_back(100.00);
      float box_sizey = 100.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);
      
      float dt=0.001;
   
      /* Mean energy of electrons when they begin, in keV */
      float E_mean = 10.0;
      
      float E_low_threshold = 10; 
      
      int timescale_red_emission = 50.; //110/dt;
      int timescale_green_emission = 50.;// 0.7/dt;
      int timescale_blue_emission = 50.;// 0.001/dt;


      float wavelength_red = 500.0;
      float wavelength_green = 600.0;
      float wavelength_blue = 800.0;
      float e_chg = 1e3; // elementary charge
      float epsilon_naught = 1e-9;
      
      vector<Electron> electrons;

      const float m_e = 1.0;
      float hc = 1.23984197 *10000; //*plank's constant in keV*nm 

      float number_density(float h) {
         return 7.86778379e6*exp(-h*0.176299767);
      }
      float bethe_energy_loss(float h,float v) {
         // dE/dx due to atmospheric 'drag'
         float effective_Z = 7.2; // 
         float Ionization = 92.0;
         float coeff = 4*M_PI*number_density(h)*pow(e_chg,4)*effective_Z / (m_e * pow(v,2));
         float parenthesis= log(  2*pow(v,2) / (m_e*Ionization  )) - 1.2329;
         return coeff*parenthesis;
      }

      bool in_sheet1(float x, float y) {
         int N = 4;
         float shift = 0.1*(float)rand()/RAND_MAX + 0.95;
         if (x < box_sizex*sin(shift * N * M_PI * (y - box_sizey/2) / box_sizey) 
             && x > box_sizex*sin(shift * N * M_PI * (y - box_sizey/2) / box_sizey) ) {
          return 1  ;
         } else {
            return 0;
         }

      }

 };


class OutputReport {
  public:
    vector<float> EvsT;
    int interactions;
    int respawns;


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

class Electron {
public:
   int interaction_count=0;
   int respawn_count=0;
   float x;
   Simulation* sim;
   float y;
   float z;
   float vx;
   float vy;
   float vz;
   float E=0;
   bool emitting;
   float emitting_time_left;
   float emitting_wavelength;
   float Fx;
   float Fy;
   float Fz;
   float t;
   float ID;
   float tmp;
   vector<float>  rands;
   

   float p_emit=0;
   float p_emit_r=0;
   float p_emit_g=0;
   float p_emit_b=0;

    
   float get_p_emit_red(float x) {
      return 
      1.97692801419e-17 * pow(x, 8 ) + 
     -1.24675920355e-14 * pow(x, 7 ) + 
     2.61694799395e-12 * pow(x, 6 ) + 
     -1.20009942187e-10 * pow(x, 5 ) + 
     -3.03037545164e-08 * pow(x, 4 ) + 
     4.58332800758e-06 * pow(x, 3 ) + 
     -0.000224277684039 * pow(x, 2 ) + 
     0.00397540318554 * pow(x, 1 ) + 
     -0.00815824073091 * pow(x, 0 ) ;
   }

   
   float get_p_emit_green(float x) {
     return 1.22309846233e-16 * pow(x, 8 ) + 
     -1.04247753054e-13 * pow(x, 7 ) + 
     3.62402159771e-11 * pow(x, 6 ) + 
     -6.54860670058e-09 * pow(x, 5 ) + 
     6.39054458804e-07 * pow(x, 4 ) + 
     -2.94617910697e-05 * pow(x, 3 ) + 
     0.000107184265124 * pow(x, 2 ) + 
     0.0318371583695 * pow(x, 1 ) + 
     0.0028840983908 * pow(x, 0 ) ;
   }

   
   float get_p_emit_blue(float x) {
      return 
         7.71746456929e-17 * pow(x, 8 ) + 
        -6.74045102078e-14 * pow(x, 7 ) + 
        2.40765868528e-11 * pow(x, 6 ) + 
        -4.49364430269e-09 * pow(x, 5 ) + 
        4.59098779432e-07 * pow(x, 4 ) + 
        -2.34482985981e-05 * pow(x, 3 ) + 
        0.000317362998942 * pow(x, 2 ) + 
        0.0130756074812 * pow(x, 1 ) + 
        0.00402470459253 * pow(x, 0 );
  
   }


   float calculate_probabilities(float h) {
        p_emit_r = get_p_emit_red(h);
        p_emit_g = get_p_emit_green(h);
        p_emit_b = get_p_emit_blue(h);
        p_emit = p_emit_r + p_emit_g + p_emit_b ;
        p_emit_r = p_emit_r / p_emit;
        p_emit_g = p_emit_g / p_emit;
        p_emit_b = p_emit_b / p_emit;
        //cout << h<<"   "   << p_emit_r << "  " << p_emit_g << "  " << p_emit_b << "  "    << p_emit << endl;

   }
   float get_p_interaction(float h) {
      return 0.4*(get_p_emit_red(h) + get_p_emit_green(h) +  get_p_emit_blue(h)) ;
   }

   //times that each emission should be active for  
   float get_t_emit_red() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_red_emission;
   }
   float get_t_emit_green() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_green_emission;
   }
   float get_t_emit_blue() {
      vector<float> rnd = gen_random(2);
      return (rnd[0]+1) * sim->timescale_blue_emission;
   }

   int random_collision() {
      vector<float> rands = gen_random(4);
      float theta = 2*M_PI*abs(rands[0]);
      vx = vx * cos(theta);
      vy = vy * sin(theta);
      vz = -sqrt( 2*E / sim->m_e - (vx*vx + vy*vy) );
    }

   int rescale_velocities() {
      //rescales the velocities to be consistent with the current kinetic energy
//      cout << "BEFORE: " << vx << "  " << vy << "  " << vz << endl;
      float R = 0.5*sim->m_e*(vx*vx + vy*vy + vz*vz) / E;
      vx = vx/sqrt(R);
      vy = vy/sqrt(R);
      vz = vz/sqrt(R);
//      cout << "BEFORE: " << vx << "  " << vy << "  " << vz << endl;
//      cout << " " << endl;
      

   }

   
   int reset() {
      if (ID==3) {
         cout << "Respawning" << endl;
         cout << "  at respawn: Fx=" << Fx << "   Fy=" << Fy << "   Fz="<<Fz << endl;
         cout << "  at respawn:  x=" <<  x << "    y=" <<  y << "    z="<< z << endl;
      }
      t = 0;
      x = sim->box_sizex *(float)rand() / RAND_MAX;
      y = sim->box_sizey* (float)rand() / RAND_MAX;

       
     
      randoms = gen_random(6);
      
      x = (randoms[3])*75.0 + sim->box_sizex/2.0;
      y = abs(randoms[4])*75.0;

/*      if ((float)rand()/RAND_MAX > 0.9) {
         y = 75.0;
      } else {
         y = 30.0;
      }
      x = 50.0;

      */
      

      //
      //z = 0.1 * sim->box_sizez* (float)rand() / RAND_MAX + 0.9*sim->box_sizez;
      z = 1.1*sim->box_sizez -0.0001;

  //    z = 180.0 ; //sim->box_sizez*(float)rand()/RAND_MAX ;  

      E = 100000*sim->E_mean * (abs(randoms[0])+1.0);
      vx = (0.0001 * sqrt(2*E/sim->m_e) * randoms[1]);
      vy = (0.0001 * sqrt(2*E/sim->m_e) * randoms[2]);
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -sqrt( 2*E / sim->m_e - tmp ) ;
      emitting = 0;
      emitting_time_left = 0;
      emitting_wavelength = 0;
      Fx = 0; 
      Fy = 0; 
      Fz = 0;
      respawn_count +=1;

   }


};


class PhotonDensity {


   public:
      int resolution_x = 100; //divisions along the box dimension
      int resolution_y = 400;
      int resolution_z = 800;

      float optical_decay_power = 2.0; //light intensity drops off as distance to this power (vacuum would be 2.0) 

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

        cout << "\rBrightest voxel had un-normed brightness of " << gmax << "                             ";
        cout.flush();

        *R = scalar_divide(R,gmax); 
        *G = scalar_divide(G,gmax); 
        *B = scalar_divide(B,gmax); 

            

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

        float maxpixelsumR =1;
        float maxpixelsumG =1;
        float maxpixelsumB =1;

        for (int k=0; k<resolution_z; k++) {
           for (int j=0; j<resolution_y; j++) {
              pixelsumR=0;
              pixelsumG=0;
              pixelsumB=0;
              for (int i=0; i<resolution_x; i++)  {
                 opt_decay =  pow( 
                     (  1 - float(i)/float(resolution_x)  ) , optical_decay_power);
                 pixelsumR += get_element(R, i, j, k) / (0.8*opt_decay);
                 pixelsumG += get_element(G, i, j, k) / (0.9*opt_decay);
                 pixelsumB += get_element(B, i, j, k) / (opt_decay);
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


class MagneticField {
   public:
      float strength = 1000.0;

      vector<float> at(float t,float x,float y,float z){
         vector<float> r;
//         r.push_back((strength*sin(x-50)/(x-50)));
//         r.push_back((strength*sin(x-50)/(x-50)));
         r.push_back(strength);
         r.push_back(strength);
         r.push_back(0);
         return r;
      }
};


class ChargeDensity {
   public:
      int resolution_x = 20;
      int resolution_y = 20;
      int resolution_z = 60;

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
      float *Ex = new float[Nx*Ny*Nz];
      float *Ey = new float[Nx*Ny*Nz];
      float *Ez = new float[Nx*Ny*Nz];

      int init(Simulation *simm, ChargeDensity *rhoo) {
         sim = simm;
         rho = rhoo;
         Nx = rho->resolution_x ;
         Ny = rho->resolution_y ;
         Nz = rho->resolution_z ;
         Lx = sim->box_sizex ; 
         Ly = sim->box_sizey ;
         Lz = sim->box_sizez ;

         float *Ex = new float[Nx*Ny*Nz];
         float *Ey = new float[Nx*Ny*Nz];
         float *Ez = new float[Nx*Ny*Nz];
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
                  
                 
          tmpReal = -out[element][0] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);
          tmpImag = -out[element][1] / (pow(wavenumber_x,2) + pow(wavenumber_y,2) + pow(wavenumber_z,2) + 1e-20);

          tmpReal /= sim->epsilon_naught;
          tmpImag /= sim->epsilon_naught;

          inPhiX[element][1] = -tmpReal * wavenumber_x ;
          inPhiX[element][0] =  tmpImag * wavenumber_x ;
          
          inPhiY[element][1] = -tmpReal * wavenumber_y ;
          inPhiY[element][0] =  tmpImag * wavenumber_y ;
          
          inPhiZ[element][1] = -tmpReal * wavenumber_z ;
          inPhiZ[element][0] =  tmpImag * wavenumber_z ;

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


}

         









};




#endif

