/*
    This file is part of Wavy.

    Wavy is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Wavy is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Wavy.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2013, Cedric Beny
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <fftw3.h>

typedef struct {
  int x, y;
} icoo;

typedef struct {
  double x, y;
} coo;

typedef struct {
  int      nx, ny;         // dimensions of the 2d grid
  double   lx, ly;         // lx = nx * dx etc...
  int      N;              // N = nx*ny
  fftw_complex  *buf;      // complex wavefunction 
  double  *prob;           // probability density (derived from buf)
  double   dx, dy, dt;     // lattice spacing and time step
  double   mass;

  double   norm, max;

  double  *pot;            // energetic potential
 
  fftw_complex  *prob_fft; // stores fft of prob

  fftw_plan forward, backward;
  fftw_plan prob_fw, prob_bw;

  double time_unit, max_dt;
  
  coo pos, vel;            // expected position and velocity
} particle;


void particle_init( particle *f, double size, int precx, int precy );
void particle_fourierize( particle *f );
void particle_unfourierize( particle *f );
void prob_fourierize( particle *f );
void prob_unfourierize( particle *f );
void particle_add_coherent( particle *f, double x0, double y0, double vx, double vy, double omega );
void particle_make_probability( particle *f );
void particle_normalize( particle *f );
void particle_represent( particle *f, unsigned char *pixmap, int disp_pot );
void move_potential( particle *f );
void move_kinetic( particle *f );
double variance( particle *f );
icoo measure_position( particle *f, double sigma, fftw_complex * meas_fft );

#endif
