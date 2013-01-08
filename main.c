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
#include <time.h>
#include <fftw3.h>
#include <GL/glut.h>
#include <sys/time.h>

double frand( double a, double b ) {
  return (b-a)*(random()/(RAND_MAX+1.0)) + a;
}

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
} particle;


particle wave;
GLuint texobj;
unsigned char *pixmap;
double omega;
int observe;
struct timeval realtime;
fftw_complex * meas_fft;  // fourier transform of the POVM element, ready for convolution
double meas_sigma;

int disp_pot;
struct timeval disp_pot_time;

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

  fftw_plan_with_nthreads( 2 ); 

  f->forward = fftw_plan_dft_2d( f->ny, f->nx, f->buf, f->buf, FFTW_FORWARD, FFTW_MEASURE );
  f->backward = fftw_plan_dft_2d( f->ny, f->nx, f->buf, f->buf, FFTW_BACKWARD, FFTW_MEASURE );

  f->prob_fw = fftw_plan_dft_r2c_2d( f->ny, f->nx, f->prob, f->prob_fft, FFTW_MEASURE );
  f->prob_bw = fftw_plan_dft_c2r_2d( f->ny, f->nx, f->prob_fft, f->prob, FFTW_MEASURE );

  int i;
  for( i = 0; i < f->N; i++ ) {
    f->buf[i][0] = f->buf[i][1] = 0.0;
  }

  // some reasonable default settings

  wave.mass = 10.0;
  wave.time_unit = 2.0; 
  wave.max_dt = 0.2;  
  wave.dt = 0.05;
}

inline void particle_fourierize( particle *f ) {
  fftw_execute( f->forward );
}

inline void particle_unfourierize( particle *f ) {
  fftw_execute( f->backward );

  double Z = f->nx * f->ny;
  int i;
  for( i = 0; i < f->N; i++ ) {
    f->buf[i][0] /= Z;
    f->buf[i][1] /= Z;
  }
}

inline void prob_fourierize( particle *f ) {
  fftw_execute( f->prob_fw );
}

inline void prob_unfourierize( particle *f ) {
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

/* Displays the wavefunction in a window */

#define bound(x,a,b)  ( (x) < (a) ? (a) : (x) > (b) ? (b) : (x) )

void display_hook() { 

  // create the 2D image in 24 bits RGB array `pixmap'

  particle *f = &wave;

  particle_make_probability( f );

  double max = f->max * 0.8;
  double smax = sqrt( max );

  int i, j;
  for( j = 0; j < f->ny; j++ ) {
    for( i = 0; i < f->nx; i++ ) {
      int n = i + j * f->nx;

      double real = f->buf[n][0] / smax;
      //double imag = f->buf[n][1] / smax;
      double val = f->prob[n] / max;
      
      double p = ( real + 1 ) / 2;

      int pot = 0;
      if( disp_pot ) {
	pot = ( sinf( f->pot[n] / 3.0 ) + 1.0 ) * 15;
      }

      //double p = fabs( atan2( real, imag )) / M_PI;  // makes more sense but is slow.. should precalculate
      int r = val * p * 255.0;
      int g = val * ( 1 - p ) * 110.0;
      int b = val * ( 1 - p ) * 255.0;
      
      pixmap[ 3*n + 0 ] = bound( r, pot, 255 );
      pixmap[ 3*n + 1 ] = bound( g, pot, 255 );
      pixmap[ 3*n + 2 ] = bound( b, pot, 255 );
    }
  }

  // display the particle's pixmap

  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, wave.nx, wave.ny, 0, GL_RGB, GL_UNSIGNED_BYTE, pixmap ); 

  glBindTexture( GL_TEXTURE_2D, texobj );
  glBegin( GL_QUADS );
  glTexCoord2f( 0.0, 0.0 );
  glVertex2f( -1.0, -1.0 );
  glTexCoord2f( 1.0, 0.0 );
  glVertex2f( 1.0, -1.0 );
  glTexCoord2f( 1.0, 1.0 );
  glVertex2f( 1.0, 1.0 );
  glTexCoord2f( 0.0, 1.0 );
  glVertex2f( -1.0, 1.0 );
  glEnd();

  glFlush();
  glutSwapBuffers();
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

void measure_position( particle *f ) {
  double sigma = meas_sigma;

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

  int i0, j0;
  for( ;; ) {
    i0 = f->nx * ( random() / ( RAND_MAX + 1.0 ));
    j0 = f->ny * ( random() / ( RAND_MAX + 1.0 ));
    if( frand( 0, max ) < f->prob[ i0 + j0*f->nx ]) break;
  }

  // and collapse the wavefunction

  for( j = -f->ny/2; j < f->ny/2; j++ ) {
    for( i = -f->nx/2; i < f->nx/2; i++ ) {

      double x = i * f->dx;
      double y = j * f->dy;
      double amp = exp( -( x*x + y*y ) / ( 2 * s2 ));  // FIXME shouldn't it be 4 in instead of 2??

      int i1 = ( i0 + i + f->nx ) % f->nx;
      int j1 = ( j0 + j + f->ny ) % f->ny;

      int n = i1 + j1*f->nx;
      f->buf[n][0] *= amp;
      f->buf[n][1] *= amp;      
    }
  }

  particle_make_probability( f );
  particle_normalize( f );
}

void mouse_hook( int button, int state, int x, int y ) 
{
  if( state == GLUT_DOWN ) observe = 1;
  if( state == GLUT_UP ) observe = 0; 
}

void update_time_step() {
  struct timeval now;
  gettimeofday( &now, 0 );
  wave.dt = (( now.tv_sec - realtime.tv_sec ) + ( now.tv_usec - realtime.tv_usec ) * 1e-6 ) * wave.time_unit;
  realtime = now;
  wave.dt = bound( wave.dt, 0.0, wave.max_dt );

}

void idle_hook() { 

  // stop displaying the potential after half a second
  if( disp_pot ) {
    struct timeval now;
    gettimeofday( &now, 0 );
    if(( now.tv_sec - disp_pot_time.tv_sec ) + ( now.tv_usec - disp_pot_time.tv_usec ) * 1e-6 > 1.02 ) {
      disp_pot = 0;
    }  
  }

  // match simulation time step to real time elapsed since last update
  update_time_step();

  // observe the particle 
  if( observe ) {
    measure_position( &wave );
  }

  // implement potential step
  move_potential( &wave );

  // go to fourier space
  particle_fourierize( &wave );
  
  // implement kinetic step
  move_kinetic( &wave );
  
  // go back to real space
  particle_unfourierize( &wave );

  // tell GLUT to refresh the display
  glutPostRedisplay();
}

void keyboard_hook( unsigned char key, int x, int y ) {
  
  // digit keys select a potential

  if( (int)key >= 48 && (int)key <= 57 ) {

    int i, j;
    int n = 0;

    // briefly display the potential
    gettimeofday( &disp_pot_time, 0 );
    disp_pot = 1;
    
    for( j = 0; j < wave.ny; j++ ) {
      for( i = 0; i < wave.nx; i++ ) {

	double x = ( i - wave.nx/2 ) * wave.dx;
	double y = ( j - wave.ny/2 ) * wave.dy;
	double r2 = x*x + y*y;

	wave.pot[n] = 0.0;

	switch( key ) {
	case '0': // free particle
	  break;

	case '1': // harmonic potential
	  wave.pot[n] = r2 / 2 * wave.mass * omega * omega;
	  break;

	case '2': // quartic potential
	  wave.pot[n] = r2 * r2 * 0.1;
	  break;

	case '3': // logarithmic potential
	  wave.pot[n] = r2 == 0.0 ? 0.0 : 60 * log(sqrt(r2));
	  break;

	case '4': // 1/r potential
	  wave.pot[n] = r2 == 0.0 ? 0.0 : - 150.0 / sqrt(r2);
	  break;

	case '5': // 1/r impurity
	  wave.pot[n] = r2 == 0.0 ? 0.0 : 1.0 / sqrt(r2);
	  break;

	case '7': // box with impurity
	  wave.pot[n] = r2 == 0.0 ? 0.0 : 1.0 / sqrt(r2);

	case '6': // smooth box
	  {
	    double a, e = 4.0, b = 4.0;
	    a = wave.ly * 0.4;
	    wave.pot[ i + j*wave.nx ] += (y > a ? pow(b*(y-a),e) : (y < -a ? pow(b*(-a-y),e) : 0.0)) * 10.0;
	    a = wave.lx * 0.4;
	    wave.pot[ i + j*wave.nx ] += (x > a ? pow(b*(x-a),e) : (x < -a ? pow(b*(-a-x),e) : 0.0)) * 10.0;
	  }
	  break;
	}

	n++;
      }
    }
  } else {
    switch( key ) {
    case '=':
    case '+':
      wave.time_unit *= 1.5; break;
    case '-':
    case '_':
      wave.time_unit /= 1.5; break;
    }
  }
}

int main( int argcp, char **argv ) {

  // initialize multithreading for the fast fourier transform
 
  if( !fftw_init_threads()) {
    printf( "Problem initializing FFTW threads!\n" );
  }

  // seed the random generator with time in second
	
  srandom( time( 0 ));

  // create the particle 

  particle_init( &wave, 
	      10.0,     // size of the toroidal box
	      256,      // number of vertices across the box vertically
	      256 );    // number of vertices across the box horizontally

  wave.mass = 10.0;
  wave.time_unit = 3.0; // set the simulation time spent per real second
  wave.max_dt = 0.2;    // time step beyond which we give up real time

  // the wavefunction is not being observed initially

  observe = 0;

  // create the initial potential

  omega = 1.0;
  keyboard_hook( '2', 0, 0 );
  disp_pot = 0;

  // Create the initial wave function

  particle_add_coherent( &wave,  wave.lx/2, wave.ly/2 + wave.ly/4, 2.0, 0.5, omega );

  // Prepare the POVM measurement
  {
    meas_sigma = 2.0; // precision of the measurement
    
    int i, j;
    double *aux = fftw_malloc( sizeof(double) * wave.N );
    meas_fft = fftw_malloc( sizeof(fftw_complex) * wave.ny * ( wave.nx/2 + 1 ));
    for( j = -wave.ny/2; j < wave.ny/2; j++ ) {
      for( i = -wave.nx/2; i < wave.nx/2; i++ ) {
	double x = i * wave.dx;
	double y = j * wave.dy;
	int n = (( i + wave.nx ) % wave.nx) + (( j + wave.ny ) % wave.ny) * wave.nx;
	double amp = exp( -( x*x + y*y ) / ( meas_sigma * meas_sigma )); // FIXME: what of the factor 2?
	aux[n] = amp;
      }
    }

    fftw_plan plan = fftw_plan_dft_r2c_2d( wave.ny, wave.nx, aux, meas_fft, FFTW_ESTIMATE );
    fftw_execute( plan );
    // FIXME: free the data used buy this plan
  }

  // Initialize GLUT for the display and event handling

  double resol = 800;

  glutInit( &argcp, argv ); 
  glutInitWindowSize( resol * wave.nx/wave.ny, resol ); 
  glutInitWindowPosition( 0, 0 );
  
  glutCreateWindow( "Quantum wavefunction in 2 dimensions" );
  
  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
 
  glutDisplayFunc( display_hook );
  glutIdleFunc( idle_hook );
  glutMouseFunc( mouse_hook );
  glutKeyboardFunc( keyboard_hook );

  // Allocate the 2D image

  pixmap = malloc( sizeof( unsigned char ) * 3 * wave.N );
  if( !pixmap ) {
    printf( "Could not allocate memory for the display, aborting.\n" );
    exit( 1 );
  }

  // Which will be used as texture so we can scale it
  
  glEnable( GL_TEXTURE_2D );
  glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_FASTEST );
  glGenTextures( 1, &texobj );
  glBindTexture( GL_TEXTURE_2D, texobj );  
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

  // get real time

  gettimeofday( &realtime, 0 );

  // Start the GLUT event handling

  glutMainLoop(); 

  return 0;
};

