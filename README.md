Wavy
====

Real-time simulation of a two-dimensional quantum particle using the split-step method. 

<img src="http://dl.dropbox.com/u/10480705/wavy_screenshot1.jpg">

The color is brighter where the particle is more likely to be found. The color hue encodes the phase information. 


Usage
-----

If any mouse button is pressed, the particle's position is continuously monitored. This forces its localization, but also randomizes its position and speed to some extent. 

The '+' and '-' keys accelarate or slow down the simulation.

The Number keys from 0 to 3 select different potentials.


Required libraries
------------------

This software requires the fftw and gl/glut libraries.


How to build on Linux
---------------------

Type 'make' in the source directory

This creates the stand-alone binary called 'wavy'

