#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <thread>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */

//my own stuff
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
   sim.init();
   PhotonDensity photon_density;
   photon_density.init(&sim);
   EnergyDensity energy_density;
   energy_density.init(photon_density.resolution_z);
   MagneticField mag_field;
   ChargeDensity rho;
   ElectricField E_field;
   E_field.init(&sim, &rho);
   
   int voxelx,voxely,voxelz;
   int Evoxelx,Evoxely,Evoxelz;
   vector<float> B;
   float By=-1.0;
   float Bx=0.0;
   float Bz=0.0;
   float tmpfloat;


   for (int i=0; i<sim.N; i++) {
      Electron ee;
      ee.sim = &sim;
      ee.reset();
      ee.ID = i;
      sim.electrons.push_back(ee);

   }
   cout << "Electron objects created" << endl;
 
   ofstream debug_xyz;

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
         energy_density.reset();

      }

cout << sim.N << endl;

      for (int i=0; i<sim.N ; i++) {

         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
         
         //if electron has left the cell, or has not enough energy remaining, we're done tracking it.
         //Reset it.
         while (  
            e->z!=e->z //(that will evaluate true if e->z == nan) 
            || e->z < 0 
            || e->z >= 1.05*(sim.box_sizez)
            || e->E <= (sim.hc / e->sim->wavelength_red )  
          ) {
            e->reset();
            e->t = t;
         }

          
         voxelx = int((float)photon_density.resolution_x / (float)sim.box_sizex * e->x);
         voxely = int((float)photon_density.resolution_y / (float)sim.box_sizey * e->y);
         voxelz = int((float)photon_density.resolution_z / (float)sim.box_sizez * e->z);
         Evoxelx = voxelx;
         Evoxely = voxely;
         Evoxelz = voxelz;
        
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
         if (e->emitting==1 && e->dead_counter == 0) {
            e->E -= sim.E_loss_factor*sim.hc / e->emitting_wavelength;
            energy_density.increment(voxelz,sim.E_loss_factor*sim.hc / e->emitting_wavelength);
            e->emitting_time_left -= 1;
            
            if (  e->z < sim.box_sizez
                  && e->z >=0
                  && e->x >= 0 
                  && e->x <= sim.box_sizex
                  && e->y >= 0
                  && e->y <= sim.box_sizey    ) {

               if (e->emitting_wavelength == sim.wavelength_red ) {
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.8) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.16) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.4) ;
               } else if (e->emitting_wavelength == sim.wavelength_green ) {
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.18) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.69) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.4) ;
               } else if (e->emitting_wavelength == sim.wavelength_blue ) {
                  photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz,0.4) ;
                  photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz,0.4) ;
                  photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz,0.2) ;
               }

            }

            if (e->emitting_time_left < 1) { //if it's done emitting
               e->emitting = 0;
               e->emitting_wavelength = 0;
            }
         }


//CHARGE DENSITY CALCULATION///////////////////////////////
        voxelx = int((float)rho.resolution_x / (float)sim.box_sizex * e->x); 
        voxely = int((float)rho.resolution_y / (float)sim.box_sizey * e->y);
        voxelz = int((float)rho.resolution_z / (float)sim.box_sizez * e->z);
        if (
             voxelx >= 0 
          && voxelx<rho.resolution_x 
          && voxely >= 0 
          && voxely<rho.resolution_y 
          && voxelz >= 0 
          && voxelz<rho.resolution_z) 
        {
            rho.incr_element(voxelx,voxely,voxelz);  
        }
///////////////////////////////////////

         
         
//FORCES////////////////////////
         e->Fx = 0;
         e->Fy = 0;
         e->Fz = 0;
         //magnetic force F = v x B
         e->Fx += (e->vy * Bz - e->vz * By) / 1e-1 ;
         e->Fy += (e->vz * Bx - e->vx * Bz) / 1e-1;
         e->Fz += (e->vx * By - e->vy * Bx) / 1e-1 ;
         //electric force F = q E
         e->Fx += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ex,voxelx,voxely,voxelz));
         e->Fy += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ey,voxelx,voxely,voxelz));
         e->Fz += 1.0e12/sim.N*(sim.e_chg * E_field.get_element(E_field.Ez,voxelx,voxely,voxelz));
////////////////////////////////
         

//INTEGRATION////////////////////////////////////
         //Perform equation of motion integration:
         e->vx +=  (e->Fx / sim.m_e) * sim.dt ;
         e->vy +=  (e->Fy / sim.m_e) * sim.dt ;
         e->vz +=  (e->Fz / sim.m_e) * sim.dt ;

         e->x += e->vx * sim.dt;
         e->y += e->vy * sim.dt;
         e->z += e->vz * sim.dt;
//////////////////////////////////////////////////
         


///////////PERIODIC BOUNDARY CONDITIONS/////
           while (e->x < 0) {e->x += sim.box_sizex;}
           while (e->y < 0) {e->y += sim.box_sizey;}
           while (e->x > sim.box_sizex) {e->x -= sim.box_sizex;}
           while (e->y > sim.box_sizey) {e->y -= sim.box_sizey;}
///////////////////////////////////////////           


//print out a bunch of info for one of the particles so we can see how simulation is progressing 
         if (e->ID==3) { 
            cout << t 
                 << "\t" 
                 << e->x 
                 << "\t" 
                 << e->y 
                 << "\t" 
                 << e->z 
                 << "\tE: " 
                 << e->E 
                 <<  "\t"  
                 << e->Fx 
                 << "\t" 
                 << e->Fy 
                 << "\t" 
                 << e->Fz 
                 << "\t" 
                 << e->vx 
                 << "\t" 
                 << e->vy 
                 << "\t" 
                 << e->vz 
                 << "\n" ;
         }
         
      } // end of loop over electrons

      //recompute the electric field
      if (t%sim.E_field_recalc==0 ) {
         E_field.compute();
      }
      rho.reset(); //zero out the charge density for next time


   } //end of loop over time
   
   cout << endl;

   cout << "\nCLEANING UP" << endl;

   energy_density.cleanup();
   photon_density.cleanup();
   E_field.cleanup();
   rho.cleanup();

   cout << "\n\n" << endl;

   return 0;

}








