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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <fftw3.h>

#include "particle.h"

inline double frand( double a, double b ) {
  return (b-a)*(random()/(RAND_MAX+1.0)) + a;
}


/* Initialize a particle structure, with zero wavefunction */

void particle_init( particle *f, double size, int precx, int precy ) {
  f->nx = precx; 
  f->ny = precy;
  f->dx = size / f->nx;

  f->dy = f->dx;
  f->lx = f->nx * f->dx;
  f->ly = f->ny * f->dy;
  f->N = f->nx * f->ny;

  f->buf = fftw_malloc( sizeof(fftw_complex) * f->N );
  f->prob = fftw_malloc( sizeof( double ) * f->N );
  f->pot = calloc( f->N, sizeof( double ));
  f->prob_fft = fftw_malloc( sizeof(fftw_complex) * f->ny * ( f->nx / 2 + 1 ));
  if( !f->buf || !f->prob || !f->pot || !f->prob_fft ) {
    printf( "Not enough memory! Aborting.\n" );
    exit( 1 );
  }


  f->forward = fftw_plan_dft_2d( f->ny, f->nx, f->buf, f->buf, FFTW_FORWARD, FFTW_MEASURE );
  f->backward = fftw_plan_dft_2d( f->ny, f->nx, f->buf, f->buf, FFTW_BACKWARD, FFTW_MEASURE );

  f->prob_fw = fftw_plan_dft_r2c_2d( f->ny, f->nx, f->prob, f->prob_fft, FFTW_MEASURE );
  f->prob_bw = fftw_plan_dft_c2r_2d( f->ny, f->nx, f->prob_fft, f->prob, FFTW_MEASURE );

  int i;
  for( i = 0; i < f->N; i++ ) {
    f->buf[i][0] = f->buf[i][1] = 0.0;
  }

  // some reasonable default settings

  f->mass = 10.0;
  f->time_unit = 2.0; 
  f->max_dt = 0.2;  
  f->dt = 0.05;
}

/* Go from position to momentum representation... */
 
void particle_fourierize( particle *f ) {
  fftw_execute( f->forward );
}

/* ... and back */

void particle_unfourierize( particle *f ) {
  fftw_execute( f->backward );

  double Z = f->nx * f->ny;
  int i;
  for( i = 0; i < f->N; i++ ) {
    f->buf[i][0] /= Z;
    f->buf[i][1] /= Z;
  }
}

/* Compute the Fourier transform of the probability density */ 

void prob_fourierize( particle *f ) {
  fftw_execute( f->prob_fw );
}

void prob_unfourierize( particle *f ) {
  fftw_execute( f->prob_bw );

  double Z = f->nx * f->ny;
  int i;
  for( i = 0; i < f->N; i++ ) {
    f->prob[i] /= Z;
  }
}


/* Superpose a coherent state to the particle's present state */

void particle_add_coherent( particle *f, double x0, double y0, 
			    double vx, double vy, double omega ) {
  int i, j;
 
  double Z = sqrt(( f->mass * omega ) / M_PI );

  for( j = 0; j < f->ny; j++) {
    for( i = 0; i < f->nx; i++) {

      double x = i * f->dx;
      double y = j * f->dy;
      double a = x - x0;
      double b = y - y0;

      double amp = Z * exp( -( a*a + b*b ) * ( f->mass * omega / 2 ));
      double phase = f->mass * ( vx * x + vy * y );

      f->buf[ i + j*f->nx ][0] += cos( phase ) * amp;
      f->buf[ i + j*f->nx ][1] += sin( phase ) * amp;      
    }
  }
}


/* Compute the position probability density of the particle and compute its norm 
   and max value  */

void particle_make_probability( particle *f ) {
  f->max = 0.0;
  f->norm = 0.0;

  int n;
  for( n = 0; n < f->N; n++ ) {
    double r = f->buf[n][0];
    double i = f->buf[n][1];
    
    f->prob[n] = r*r + i*i;
    if( f->prob[n] > f->max ) f->max = f->prob[n];
    f->norm += f->prob[n];
  }
  f->norm *= f->dx * f->dy;
}


/* Divides the particle's wavefunction by f->norm */

void particle_normalize( particle *f ) {
  int n;
  double sn = sqrt( f->norm );
  for( n = 0; n < f->N; n++ ) {
    f->buf[n][0] /= sn;
    f->buf[n][1] /= sn;
  }
}


/* create the 2D image in given 24 bits RGB array `pixmap' */

#define bound(x,a,b)  ( (x) < (a) ? (a) : (x) > (b) ? (b) : (x) )

void particle_represent( particle *f, unsigned char *pixmap, int disp_pot ) {
  
  particle_make_probability( f );

  float max = f->max * 0.8;
  float smax = sqrtf( max );

  // we should do these calculations in a shader...
  int i, j;
  for( j = 0; j < f->ny; j++ ) {
    for( i = 0; i < f->nx; i++ ) {
      int n = i + j * f->nx;

      float real = f->buf[n][0] / smax;
      float imag = f->buf[n][1] / smax;
      float val = f->prob[n] / max;
      

      int pot = 0;
      if( disp_pot ) {
        pot = ( sinf( f->pot[n] / 3.0 ) + 1.0 ) * 7;
      }

      //float p = ( real + 1 ) / 2;   // faster
      float p = fabsf( atan2f( real, imag )) / (M_PI);  
      int r = sqrtf(val * p * 1.0) * 255.0;
      int g = sqrtf(val * ( 1 - p ) * 0.3) * 255.0;
      int b = sqrtf(val * ( 1 - p ) * 1.0) * 255.0;
      //int r = val * sqrtf(cosf(p)) * ( 255.0 + 100 );
      //int g = val * sqrtf(1-cosf(p)) * ( 110.0 + 100 );
      //int b = val * sqrtf(1-cosf(p)) * ( 255.0 + 0 );
      
      pixmap[ 3*n + 0 ] = bound( r, pot, 255 );
      pixmap[ 3*n + 1 ] = bound( g, pot, 255 );
      pixmap[ 3*n + 2 ] = bound( b, pot, 255 );
    }
  }
}

/* Evolves the wavefunction in position representation for a short time 
   with only the potential term of the Hamiltonian. */

void move_potential( particle *f ) {

  int i, j;
  int n = 0;

  for( j = 0; j < f->ny; j++ ) {
    for( i = 0; i < f->nx; i++ ) {

      float phase = f->dt * f->pot[n];
      double c = cosf( phase );
      double s = sinf( phase );

      double a = f->buf[n][0];
      double b = f->buf[n][1];

      f->buf[n][0] = a * c - b * s;
      f->buf[n][1] = a * s + b * c;

      n++;
    }
  }
}


/* Evolves the wavefunction in momentum representation for a short time 
   with only kinetic term of the Hamiltonian. */

void move_kinetic( particle *f ) {

  int i, j;
  int n = 0;
  double fac = sqrt( f->dt / ( 2 * f->mass ));
  double dpx =  2 * M_PI / f->lx * fac;
  double dpy =  2 * M_PI / f->lx * fac;

  // FIXME: we could precaculate some stuff here

  for( j = 0; j < f->ny; j++ ) {
    double py = ( j > f->ny/2 ? j - f->ny : j ) * dpy;
    double py2 = py * py;

    for( i = 0; i < f->nx; i++ ) {
      double px = ( i > f->nx/2 ? i - f->nx : i ) * dpx;
         
      float phase = px * px + py2;
      double c = cosf( phase );
      double s = sinf( phase );

      double a = f->buf[n][0];
      double b = f->buf[n][1];

      f->buf[n][0] = a * c - b * s;
      f->buf[n][1] = a * s + b * c;

      n++;
    }
  }
}

/* computes the variance of the f->prob array */

double variance( particle *f ) {

  if( !f->norm ) {
    f->norm = 0.0;
    int n;
    for( n = 0; n < f->N; n++ ) f->norm += f->prob[n];
    f->norm *= f->dx * f->dy;
  }

  double xave = 0.0, yave = 0.0;
  int i, j;
  for( j = 0; j < f->ny; j++ ) {
    for( i = 0; i < f->nx; i++ ) {
      int n = i + j*f->nx;
      xave += i * f->dx * f->prob[n];
      yave += j * f->dy * f->prob[n];
    }
  }
  xave *= f->dx * f->dy / f->norm;
  yave *= f->dx * f->dy / f->norm;
  double r2ave = 0.0;

  for( j = 0; j < f->ny; j++ ) {
    for( i = 0; i < f->nx; i++ ) {
      int n = i + j*f->nx;
      double x = i * f->dx - xave;
      double y = j * f->dy - yave;
      r2ave += ( x*x + y*y ) * f->prob[n];
    }
  } 
  r2ave *= f->dx * f->dy / f->norm;

  return r2ave;
}

/* Implements a measurement of position with precision sigma */

icoo measure_position( particle *f, double sigma, 
    fftw_complex * meas_fft ) {

  int i, j, n;
  double s2 = sigma * sigma;

  // compute the exact position probabilities

  for( n = 0; n < f->N; n++ ) {
    double r = f->buf[n][0];
    double i = f->buf[n][1];
    f->prob[n] = r*r + i*i;
  }

  // get the probs for the approx position POVM by convolution

  prob_fourierize( f );
  for( n = 0; n < f->ny * ( f->nx / 2 + 1 ); n++ ) {
    double a = f->prob_fft[n][0];
    double b = f->prob_fft[n][1];
    double x = meas_fft[n][0];
    double y = meas_fft[n][1];
    f->prob_fft[n][0] = a * x - b * y;
    f->prob_fft[n][1] = a * y + b * x;
  }
  prob_unfourierize( f );

  // then select a position randomly according to the distribution

  double max = 0.0;
  for( n = 0; n < f->N; n++ ) {
    if( f->prob[n] > max ) max = f->prob[n];
  }

  icoo outcome;
  for( ;; ) {
    outcome.x = f->nx * ( random() / ( RAND_MAX + 1.0 ));
    outcome.y = f->ny * ( random() / ( RAND_MAX + 1.0 ));
    if( frand( 0, max ) < f->prob[ outcome.x + outcome.y*f->nx ]) break;
  }

  // and collapse the wavefunction

  for( j = -f->ny/2; j < f->ny/2; j++ ) {
    for( i = -f->nx/2; i < f->nx/2; i++ ) {

      double x = i * f->dx;
      double y = j * f->dy;
      double amp = exp( -( x*x + y*y ) / ( 2 * s2 ));  // FIXME shouldn't it be 4 in instead of 2??

      int i1 = ( outcome.x + i + f->nx ) % f->nx;
      int j1 = ( outcome.y + j + f->ny ) % f->ny;

      int n = i1 + j1*f->nx;
      f->buf[n][0] *= amp;
      f->buf[n][1] *= amp;      
    }
  }

  particle_make_probability( f );
  particle_normalize( f );

  return outcome;
}
