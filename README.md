Wavy
====

Real-time simulation of a two-dimensional quantum particle using the split-step method. 

<img src="http://www.qimlr.org/assets/waves.jpg">

The color is brighter where the particle is more likely to be found. The color hue encodes the phase information. 


Usage
-----

If any mouse button is pressed, the particle's position is continuously monitored. This forces its localization, but also randomizes its position and speed to some extent. 

The '+' and '-' keys accelerate or slow down the simulation.

The Number keys from 0 to 3 select different potentials.

Press q to quit.


Required libraries
------------------

This software requires the multithreaded fftw3 library and gl/glut libraries.


How to build on Mac/Linux
-------------------------

Type 'make' in the source directory

This creates the stand-alone binary called 'wavy'

