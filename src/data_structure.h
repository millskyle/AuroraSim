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
      int N = 100000;
      int tmax = 1000;

      /* box size, in km */
      float box_sizex = 100.0; //.push_back(100.00);
      float box_sizey = 100.0; //.push_back(100.00);
      float box_sizez = 200.0; //push_back(200.00);
      
      float dt=0.001;
   
      /* Mean energy of electrons when they begin, in keV */
      float E_mean = 10.0;
      
      float E_low_threshold = 10; 
      
      int timescale_red_emission = 110/dt;
      int timescale_green_emission = 0.7/dt;
      int timescale_blue_emission = 0.001/dt;
      float wavelength_red = 500.0;
      float wavelength_green = 600.0;
      float wavelength_blue = 800.0;
      float e_chg = 1.0; // elementary charge
      
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
      t = 0;
      x = sim->box_sizex * (float)rand() / RAND_MAX;
      y = sim->box_sizey * (float)rand() / RAND_MAX;
      z = 0.9*sim->box_sizez -0.0001;
      randoms = gen_random(3);
      E = 1000*sim->E_mean * (abs(randoms[0])+1.0);
      vx = (0.0001 * sqrt(2*E/sim->m_e) * randoms[1]);
      vy = (0.0001 * sqrt(2*E/sim->m_e) * randoms[2]);
      tmp = vx*vx + vy*vy; //calculate how much energy is taken by x,y velocity
      vz = -sqrt( 2*E / sim->m_e - tmp   ) ;
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
      unsigned int resolution_x = 250; //divisions along the box dimension
      unsigned int resolution_y = 250;
      unsigned int resolution_z = 500;

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
      float strength = 1.0;

      vector<float> at(float t,float x,float y,float z){
         vector<float> r;
         r.push_back((strength*sin(x-50)/(x-50)));
         r.push_back((strength*sin(x-50)/(x-50)));
         r.push_back(0);
         return r;
      }
};


class ChargeDensity {
   public:
      int resolution_x = 100;
      int resolution_y = 100;
      int resolution_z = 200;

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
   
      int compute(Simulation *sim, ChargeDensity *rho) {

         int Nx = rho->resolution_x ;
         int Ny = rho->resolution_y ;
         int Nz = rho->resolution_z ;
         
         float Lx = sim->box_sizex ; 
         float Ly = sim->box_sizey ;
         float Lz = sim->box_sizez ;
         
         float dx = Lx/float(Nx-1) ;
         float dy = Ly/float(Ny-1) ;
         float dz = Lz/float(Nz-1) ;
   
         vector<double> in(Nx*Ny*Nz, 0);
         vector<double> out(Nx*Ny*Nz, 0);

         fftw_plan q;

         q = fftw_plan_r2r_3d(Nx,Ny,Nz,in.data(), out.data(), FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00, 0 );
         
         for (int i=0; i<Nx; i++) {
            for (int j=0; j<Ny; j++) {
               for (int k=0; k< Nz; k++) {
                  // poisson's equation: del^2 phi = - rho / epsilon_naught
                  in[i + Nx *( j + Ny * k) ] = - (double)(rho->p[i + int(Nx) *( j + int(Ny) * k) ]);
               }
            }
         }

         cout << "\rExecuting FFTW           " ;
         cout.flush();
         fftw_execute(q);
         cout << "\rFFTW done                " ;
         cout.flush();

         for (int i=0; i<out.size(); i++) {
            if (out[i]==0) {
            } else {
            cout << out[i] << endl;
            }
         }




}

         























































};




#endif

