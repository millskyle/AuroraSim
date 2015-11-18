#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <fftw3.h>

#include "utility_functions.h"
#include "data_structure.h"

using namespace std;



int main() {
   //srand ( time(NULL));
   vector<float> randoms;
   Simulation sim;
   PhotonDensity photon_density;
   OutputReport report;
   MagneticField mag_field;
   ChargeDensity rho;
   ElectricField E_field;
   E_field.init(&sim, &rho);
   
   unsigned int voxelx,voxely,voxelz;

   vector<float> B;


   for (int i=0; i<sim.N; i++) {
      Electron ee;
      ee.sim = &sim;
      ee.reset();
      ee.ID = i;
      sim.electrons.push_back(ee);

   }
   cout << "Electron objects created" << endl;
 
   ofstream debug_xyz;
  // debug_xyz.open("out.xyz");


   float rnd;
   float p_interaction = 0.1;
   int t=0;
   cout << endl;
   cout << "Beginning to integrate: " << endl;
   for (t=0; t< sim.tmax; t++ ) {
      if (t%100==0) {
         photon_density.write_image(t);
         photon_density.reset();

      }


      if (t%5==0) { 
         cout << "\r   t = " << t << "                                                    ";
         cout.flush(); 
      }
 
      for (int i=0; i<sim.N ; i++) {
         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
           //if electron has left the cell, or has not enough energy remaining, we're done tracking it.
         //Reset it.
         while (  //e->x <= 0 || e->x >= sim.box_sizex
          //  || e->y <= 0 || e->y >= sim.box_sizey
           // || 
            e->z <= 0 //|| e->z >= sim.box_sizez
            || e->E <= (sim.hc / e->sim->wavelength_red ) ) {
            e->reset();
            e->t = t;
         }
        
        voxelx = int((float)photon_density.resolution_x / (float)sim.box_sizex * e->x);
        voxely = int((float)photon_density.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)photon_density.resolution_z / (float)sim.box_sizez * e->z);

        if ( e->z < sim.box_sizez 
          && e->x >= 0 
          && e->x <= sim.box_sizex
          && e->y >= 0
          && e->y <= sim.box_sizey    ) {
         
         e->calculate_probabilities(e->z);

         rnd = (float)rand()/RAND_MAX;
         if (rnd < e->p_emit) {
            e->interaction_count++;
            rnd = (float)rand()/RAND_MAX;
            if (rnd < e->p_emit_r) {
               //emit red
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_red();
               e->emitting_wavelength = sim.wavelength_red;
            } else if (rnd < (e->p_emit_r + e->p_emit_g)) {
               //emit green
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_green();
               e->emitting_wavelength = sim.wavelength_green;
            } else { 
               //emit blue
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_blue();
               e->emitting_wavelength = sim.wavelength_blue;
            }

         }

         //Check to see if electron will emit energy this timestep:
         if (e->emitting==1) {
            e->E -= sim.hc / e->emitting_wavelength;
            e->random_collision();
//            cout << sim.hplanck * sim.clight / e->emitting_wavelength<< endl;
            if(e->ID==3) {
               report.EvsT.push_back(e->E);
            }
            e->emitting_time_left -= 1;
            

            if (e->emitting_wavelength == sim.wavelength_red ) {
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.2126*1.0) ;
                photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.7152*79.0/255.0) ;
            } else if (e->emitting_wavelength == sim.wavelength_green ) {
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.2126*189.0/255.0) ;
                photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.7152) ;
            } else if (e->emitting_wavelength == sim.wavelength_blue ) {
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.2126*70./255.) ;
                photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.0722) ;
            }

            
            if (e->emitting_time_left < 1) { //if it's done emitting
               e->emitting = 0;
               e->emitting_wavelength = 0;
            }
         }

        }

        voxelx = int((float)rho.resolution_x / (float)sim.box_sizex * e->x);
        voxely = int((float)rho.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)rho.resolution_z / (float)sim.box_sizez * e->z);
         
        rho.incr_element(voxelx,voxely,voxelz);   



         
         e->rescale_velocities();


         //add any applicable forces to e->Fx, e->Fy, e->Fz :

         //rescale the velocities according to the new energy

         randoms = gen_random(3);

         B = mag_field.at(e->x,e->y,e->z,t);

//         cout << B[0] << " " << B[1] << "  " << B[2]  << endl;

         e->Fx += (e->vy * B[2] - e->vz * B[1]) / 1000. ;
         e->Fy += (e->vz * B[0] - e->vx * B[2]) /1000. ;
         e->Fz += (e->vx * B[1] - e->vy * B[0]) /1000.;

         e->Fx += (sim.e_chg * E_field.get_element(E_field.Ex,voxelx,voxely,voxelz));
         e->Fy += (sim.e_chg * E_field.get_element(E_field.Ey,voxelx,voxely,voxelz));
         e->Fz += (sim.e_chg * E_field.get_element(E_field.Ez,voxelx,voxely,voxelz));

         
         //Perform equation of motion integration:

         e->vx += 0.5 * (e->Fx / sim.m_e) * sim.dt ;
         e->vy += 0.5 * (e->Fy / sim.m_e) * sim.dt ;
         e->vz += 0.5 * (e->Fz / sim.m_e) * sim.dt ;

         //e->E-= sim.bethe_energy_loss(e->z, sqrt(e->vx*e->vx + e->vy*e->vy + e->vz*e->vx))*e->vz * sim.dt ;

         e->x += e->vx * sim.dt;
         e->y += e->vy * sim.dt;
         e->z += e->vz * sim.dt;
         
//if (e->ID==3) cout << endl << e->z << "   "  << e->p_emit_r << "  " << e->p_emit_g << "  " << e->p_emit_b << "  "<< e->p_emit ;
        
//         if (e->ID==3) cout << e->x << "\t" << e->y << "\t" << e->z << "\tE: " << e->E << "\n" ;
//         if (e->ID==3) cout << e->Fx << "\t" << e->Fy << "\t" << e->Fz << "\n" ;

//         cout << e->x << "\t" << e->y << "\t" << e->z << "\n" ; 
//         cout << e->vx << "\t" << e->vy << "\t" << e->vz << "\n" ; 
         
         
      }

      cout << "Computing electric field                  " << endl ;
      cout.flush();
      E_field.compute();
      rho.reset();

   }

   cout << endl;

photon_density.write_image(t);


debug_xyz.close();


report.interactions = sim.electrons[2].interaction_count;
report.respawns = sim.electrons[2].respawn_count;
report.write("output/report.log");






cout << "\n\n" << endl;

return 0;
}
