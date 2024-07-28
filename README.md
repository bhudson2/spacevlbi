spacevlbi
=========
Python package for simulating and optimising a space-based VLBI mission. This package enables multiple space telescopes to be modelled by propagating their orbital and attitude state. Elements of the spacecraft design that impact when observations can be performed can also be included in the simulation to assess and mitigate their impact on the science return of the mission (e.g. source visibility including Sun/Earth/Moon avoidance, star trackers, radiators, solar panels, ground station access times, etc.).

A ground-based array of radio antenna can also be modelled, enabling the (u,v) coverage that the full interferometer can achieve of a given source(s) to be calculated.

Although the package has been developed specifically for space VLBI applications, it can also be used more generally for modelling other types of astronomy mission and assessing the impact of the spacecraft design on the science return.

Please email benhudson@tudelft.nl if you have any questions.

Installation
------------
The latest version is available on [PyPi](https://pypi.org/project/spacevlbi/). Ensure that pip is installed and run the following command:

`pip install spacevlbi`

Installing with pip will install/update all of the required libraries automatically ([numpy](http://www.numpy.org/), [poliastro](https://www.poliastro.space/), [matplotlib](http://www.matplotlib.org/), [astropy](http://www.astropy.org/)).

Documentation
-------------
More detailed documentation for the package is available [here](https://spacevlbi.readthedocs.io/en/latest/).

Provided in the Examples folder is the script ExampleSetup.py. This script shows how the package can be used to model a VLBI array with a single space element. This example is based upon the preliminary concept for the [Black Hole Explorer (BHEX)](https://www.blackholeexplorer.org/) mission. The script ExampleSpaceTelescope.py shows how an object of the SpaceTelescope class can be defined.

License
-------
spacevlbi is licensed under GPLv3. See LICENSE.txt for more details.

Citation
--------
If you use spacevlbi in your publication, please cite: 