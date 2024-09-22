spacevlbi
=========
Python package for simulating and optimising a space-based VLBI mission. This package enables multiple space telescopes to be modelled by propagating their orbital and attitude state. Elements of the spacecraft design that impact when observations can be performed can also be included in the simulation to assess and mitigate their impact on the science return of the mission (e.g. source visibility including Sun/Earth/Moon avoidance, star trackers, radiators, solar panels, ground station access times, etc.).

A ground-based array of radio antenna can also be modelled, enabling the (u,v) coverage that the full interferometer can achieve of a given source(s) to be calculated.

Although the package has been developed specifically for space VLBI applications, it can also be used more generally for modelling other types of astronomy mission and assessing the impact of the spacecraft design on the science return.

Installation
------------
The latest version is available on [PyPi](https://pypi.org/project/spacevlbi/). Ensure that pip is installed and run the following command:

`pip install spacevlbi`

Installing with pip will install/update all of the required libraries automatically ([numpy](http://www.numpy.org/), [poliastro](https://www.poliastro.space/), [matplotlib](http://www.matplotlib.org/), [astropy](http://www.astropy.org/)).

Structure
---------
In this repository you can find:
- `spacevlbi/` directory containing the core Python libraries
- `Examples/` directory containing example scripts to demonstrate the use of the package
- `docs/` directory containing readthedocs configuration files

Documentation
-------------
More detailed documentation for the package is available [here](https://spacevlbi.readthedocs.io/en/latest/).

Provided in the Examples folder is the script ExampleSetup.py. This script shows how the package can be used to model a VLBI array with a single space element. This example is based upon the preliminary concept for the [Black Hole Explorer (BHEX)](https://www.blackholeexplorer.org/) mission. The script ExampleSpaceTelescope.py shows how an object of the SpaceTelescope class can be defined.

Author(s)
---------
This software has been developed by
**Ben Hudson** ![ORCID logo](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png) [0000-0002-3368-1864](https://orcid.org/0000-0002-3368-1864), Technische Universiteit Delft

License
-------
spacevlbi is licensed under GPLv3. See LICENSE.txt for more details.

Technische Universiteit Delft hereby disclaims all copyright interest in the program "spacevlbi". spacevlbi is a python package to simulate space-based VLBI missions written by the Author(s).  
Ben Hudson, Faculty of Aerospace Engineering, Technische Universiteit Delft.

&copy; 2024, B. Hudson

Citation
--------
If you use spacevlbi in your publication, please cite: 

Hudson, B., 2024. spacevlbi. 4TU.ResearchData. Software. https://doi.org/10.4121/21711227 

Would you like to contribute?
-----------------------------
If you have any comments, feedback, or recommendations, feel free to **reach out** by sending an email to benhudson@tudelft.nl

If you would like to contribute directly, you are welcome to **fork** this repository.
