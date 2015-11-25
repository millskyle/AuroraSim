#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */

#include "utility_functions.h"
#include "simulation.h"
#include "electrostatics.cpp"
#include "io.cpp"
#include "electron.h"

using namespace std;



int main() {
   srand ( time(NULL));
   vector<float> randoms;
   Simulation sim;
   sim.seed(); //generate random numbers that will be the same throughout program execution
   PhotonDensity photon_density;
   photon_density.init(&sim);
   EnergyDensity energy_density;
   energy_density.init(photon_density.resolution_z);
//   OutputReport report;
//   report.init();
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
     sim.t = t;
      if (t%5==0 && t>19) {
         photon_density.write_image(t);
         photon_density.reset();
         energy_density.write_out(t);

      }


   /*   if (t%5==0) { 
          cout << "\r                                             " ;
         cout << "t = " << t << "     ";
          cout << endl;
         cout.flush(); 
      }
      */
 
      for (int i=0; i<sim.N ; i++) {
         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
           //if electron has left the cell, or has not enough energy remaining, we're done tracking it.
         //Reset it.
         while (  
            e->z!=e->z //(that will evaluate true if e->z == nan) 
            ||e->z <= 0 
            || e->E <= (sim.hc / e->sim->wavelength_red ) ) 
          {
            e->reset();
            e->t = t;
         }
        
//        if (t==300) {
//           report.initial_positions_file << e->x << " " << e->y <<  endl;
//        }

        voxelx = int((float)photon_density.resolution_x / (float)sim.box_sizex * e->x);
        voxely = int((float)photon_density.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)photon_density.resolution_z / (float)sim.box_sizez * e->z);

        
         e->calculate_probabilities(e->z);

         rnd = (float)rand()/RAND_MAX;
         if (rnd/3.0 < e->p_emit) {
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
            e->E -= 100*sim.hc / e->emitting_wavelength;
            energy_density.increment(voxelz,100*sim.hc / e->emitting_wavelength);
//            e->random_collision();
//            cout << sim.hplanck * sim.clight / e->emitting_wavelength<< endl;
//            if(e->ID==0) {
//               report.EvsT.push_back(e->E);
//              }
            e->emitting_time_left -= 1;
            
        if ( e->z < sim.box_sizez
          && e->z >=0
          && e->x >= 0 
          && e->x <= sim.box_sizex
          && e->y >= 0
          && e->y <= sim.box_sizey    ) {
 

            if (e->emitting_wavelength == sim.wavelength_red ) {
//                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.2126) ;
//                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.2126) ;
//                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.2126) ;
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.9*1.0) ;
                photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.4*1.0) ;
                photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.9*79.0/255.0) ;
            } else if (e->emitting_wavelength == sim.wavelength_green ) {
//                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.7152) ;
//                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.7152) ;
//                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.7152) ;
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.9*189.0/255.0) ;
                photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.9) ;
            } else if (e->emitting_wavelength == sim.wavelength_blue ) {
//                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.0722) ;
//                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.0722) ;
//                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.0722) ;
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.9*70./255.) ;
                photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,1.8) ;
            }

        }

            
            if (e->emitting_time_left < 1) { //if it's done emitting
               e->emitting = 0;
               e->emitting_wavelength = 0;
            }
         }

        


        e->rescale_velocities();
       
        voxelx = int((float)rho.resolution_x / (float)sim.box_sizex * e->x); 
        voxely = int((float)rho.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)rho.resolution_z / (float)sim.box_sizez * e->z);
       
        if (
             voxelx >= 0 
          && voxelx<rho.resolution_x 
          && voxely >= 0 
          && voxely<rho.resolution_y 
          && voxelz >= 0 
          && voxelz<rho.resolution_z) {
        
        rho.incr_element(voxelx,voxely,voxelz);  
        
        }
         


         //add any applicable forces to e->Fx, e->Fy, e->Fz :

         e->Fx = 0;
         e->Fy = 0;
         e->Fz = 0;
         //rescale the velocities according to the new energy

//         randoms = gen_random(3);

         B = mag_field.at(e->x,e->y,e->z,t);


//         e->Fx += (e->vy * B[2] - e->vz * B[1]) / 1e-1 ;
//         e->Fy += (e->vz * B[0] - e->vx * B[2]) /1e-1;
//         e->Fz += (e->vx * B[1] - e->vy * B[0]) /1e-1 ;
 
         

         e->Fx += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ex,voxelx,voxely,voxelz));
         e->Fy += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ey,voxelx,voxely,voxelz));
         e->Fz += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ez,voxelx,voxely,voxelz));

 //        e->Fx = 0;
 //        e-> Fy = 0;
 //        e->Fz = 0;


/*         if ( e->ID == 4 ) {
            cout << "x=" << e->x << "    ey="<<e->y<<"      ez="<<e->z << endl;
            cout << "F_Ex=" <<E_field.get_element(E_field.Ex,voxelx,voxely,voxelz) << endl;
            cout << "F_Bx=" << (e->vy * B[2] - e->vz * B[1]) / 100. << endl;
         }
*/
         //cout << "F_Ex=" << (sim.e_chg * E_field.get_element(E_field.Ex,voxelx,voxely,voxelz)) << endl;
         
         //Perform equation of motion integration:

         

         e->vx += 0.5 * (e->Fx / sim.m_e) * sim.dt ;
         e->vy += 0.5 * (e->Fy / sim.m_e) * sim.dt ;
         e->vz += 0.5 * (e->Fz / sim.m_e) * sim.dt ;

         //e->E-= sim.bethe_energy_loss(e->z, sqrt(e->vx*e->vx + e->vy*e->vy + e->vz*e->vx))*e->vz * sim.dt ;

         e->x += e->vx * sim.dt;
         e->y += e->vy * sim.dt;
         e->z += e->vz * sim.dt;
         
//if (e->ID==3) cout << endl << e->z << "   "  << e->p_emit_r << "  " << e->p_emit_g << "  " << e->p_emit_b << "  "<< e->p_emit ;

//         periodic_boundary_conditions:
           while (e->x < 0) {e->x += sim.box_sizex;}
           while (e->y < 0) {e->y += sim.box_sizey;}
//           while (e->z < 0) {e->z += sim.box_sizez;}

//           while (e->x > sim.box_sizex) {e->x -= sim.box_sizex;}
//           while (e->y > sim.box_sizey) {e->y -= sim.box_sizey;}
//           while (e->z > sim.box_sizez) {e->z -= sim.box_sizez;}
           

        
         if (e->ID==3) cout << t << "\t" << e->x << "\t" << e->y << "\t" << e->z 
                            << "\tE: " << e->E <<  "\t"  << e->Fx << "\t" << e->Fy << "\t" << e->Fz 
                            << "\t" << e->vx << "\t" << e->vy << "\t" << e->vz <<  "\n" ;

//         cout << e->x << "\t" << e->y << "\t" << e->z << "\n" ; 
//         cout << e->vx << "\t" << e->vy << "\t" << e->vz << "\n" ; 
         
         
      }

      if (t%1==0 ) { E_field.compute(); }
      rho.reset();

   }

   cout << endl;

photon_density.write_image(t);


debug_xyz.close();


//report.interactions = sim.electrons[2].interaction_count;
//report.respawns = sim.electrons[2].respawn_count;
//report.write("output/report.log");

energy_density.cleanup();

photon_density.cleanup();
E_field.cleanup();
rho.cleanup();


cout << "\n\n" << endl;

return 0;
}
