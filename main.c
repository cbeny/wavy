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
#ifdef MAC
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <sys/time.h>

#include "particle.h"

#define bound(x,a,b)  ( (x) < (a) ? (a) : (x) > (b) ? (b) : (x) )


/* some global variables for the glut hook functions */

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


/* Displays the wavefunction in the current glut window */

void display_hook() { 

  // create the 2D image in 24 bits RGB array `pixmap'
  particle_represent( &wave, pixmap, disp_pot );

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


/* What happens when the user uses his mouse in our window */

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


/* When there is no user input, carry on with the simulation */

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
    measure_position( &wave, meas_sigma, meas_fft );
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


/* Some keyboard control */

void keyboard_hook( unsigned char key, int x, int y ) {
  
  switch( key ) {
	case 'q': // quit
	  printf( "exiting.\n");
	  exit( 0 );
	  break;
  }

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
  fftw_plan_with_nthreads( 2 ); 

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
  keyboard_hook( '1', 0, 0 );
  disp_pot = 0;

  printf("Welcome to Wavy. Press keys 0-4 for something to happen.\n");

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

