#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <cmath>
#include <time.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */

#include "testing.cpp"
#include "utility_functions.h"
#include "data_structure.h"

using namespace std;



int main() {
   srand ( time(NULL));
   vector<float> randoms;
   Simulation sim;
   PhotonDensity photon_density;
   
   unsigned int voxelx,voxely,voxelz;


   for (int i=0; i<sim.N; i++) {
      Electron ee;
      ee.sim = &sim;
      ee.reset();
      ee.ID = i;
      sim.electrons.push_back(ee);

   }
   cout << "Electron objects created" << endl;
 
   ofstream debug_xyz;
   debug_xyz.open("out.xyz");


   float rnd;
   float p_interaction = 0.1;


   for (int t=0; t< sim.tmax; t++ ) {
      debug_xyz << sim.N << "\n t=" << t << endl;
      

      for (int i=0; i<sim.N ; i++) {
         Electron* e = &sim.electrons[i]; //make a pointer to the electron we're dealing with
         debug_xyz << " A " << e->x << " " << e->y << " " << e->z << "\n" ;
         rnd = (float)rand()/RAND_MAX;
         if (rnd < e->get_p_interaction() ) {
            rnd = (float)rand()/RAND_MAX;
            if (rnd < e->get_p_emit_red()) {
               //emit red
               e->emitting = 1;
               e->emitting_time_left = e->get_t_emit_red();
               e->emitting_wavelength = sim.wavelength_red;
            } else if (rnd < (e->get_p_emit_red() + e->get_p_emit_green()) ) {
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
            e->E -= sim.hplanck * sim.clight / e->emitting_wavelength;
            e->emitting_time_left -= 1;

            voxelx = floor((float)photon_density.resolution_x / (float)sim.box_sizex * e->x);
            voxely = floor((float)photon_density.resolution_y / (float)sim.box_sizey * e->y);
            voxelz = floor((float)photon_density.resolution_z / (float)sim.box_sizez * e->z);



            if (e->emitting_wavelength == sim.wavelength_red ) {
                photon_density.incr_element(photon_density.R, voxelx,voxely,voxelz) ;
            } else if (e->emitting_wavelength == sim.wavelength_green ) {
                photon_density.incr_element(photon_density.G, voxelx,voxely,voxelz) ;
            } else if (e->emitting_wavelength == sim.wavelength_blue ) {
                photon_density.incr_element(photon_density.B, voxelx,voxely,voxelz) ;
            }

            
            if (e->emitting_time_left < 1) { //if it's done emitting
               e->emitting = 0;
               e->emitting_wavelength = 0;
            }
         }
         //if electron has exited the bottom of the cell, or has not enough energy remaining.
         if (e->z < 0 || e->E < (sim.hplanck * sim.clight / e->sim->wavelength_red ) ) {
            e->reset();
            e->t = t;
//            cout << "resetting electron" << endl;
         }


         //add any applicable forces to e->Fx, e->Fy, e->Fz :

         
         
         //Perform equation of motion integration:

         e->vx += 0.5 * (e->Fx / sim.m_e) * sim.dt*sim.dt ;
         e->vy += 0.5 * (e->Fy / sim.m_e) * sim.dt*sim.dt ;
         e->vz += 0.5 * (e->Fz / sim.m_e) * sim.dt*sim.dt ;

         e->x += e->vx * sim.dt;
         e->y += e->vy * sim.dt;
         e->z += e->vz * sim.dt;

         
         
         
         
      }

   }




debug_xyz.close();











return 0;
}
