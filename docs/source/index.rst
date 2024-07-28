spacevlbi
=========

**spacevlbi** is a Python package for simulating and optimising a space-based 
VLBI mission. This package enables multiple space telescopes to be modelled by 
propagating their orbital and attitude state. Elements of the spacecraft 
design that impact when observations can be performed can also be included in 
the simulation to assess and mitigate their impact on the science return of 
the mission (e.g. source visibility including Sun/Earth/Moon avoidance, star 
trackers, radiators, solar panels, ground station access times, etc.).

A ground-based array of radio antenna can also be modelled, enabling the (u,v) 
coverage that the full interferometer can achieve of a given source(s) to be 
calculated.

Although the package has been developed specifically for space VLBI 
applications, it can also be used more generally for modelling other types of 
astronomy mission and assessing the impact of the spacecraft design on the 
science return.

See the :doc:`Guide` section for further information.

.. note::

   This project is under active development. Please email benhudson@tudelft.nl 
   if you have any questions or submit a pull request on the git repository.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Guide
   Station
   TimeLoop
   Orbit
   Attitude
   Constraints
   Observation
   Figures
   Optimisation




