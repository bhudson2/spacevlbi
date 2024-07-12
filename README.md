spacevlbi
=========

Python package for simulating and optimising a space-based VLBI mission. This package enables multiple space telescopes to be modelled by propagating their orbital and attitude state. Elements of the spacecraft design that impact when observations can be performed can also be included in the simulation to assess and mitigate their impact on the science return of the mission (e.g. science antenna pointing, star trackers, radiators, solar panels, ground station access times, etc.).

A ground-based array of radio antenna can also be modelled, enabling the (u,v) coverage that the full interferometer can achieve of a given source(s) to be calculated.

Please email benhudson@tudelft.nl if you have any questions.

Installation
------------

The latest version is available on `PyPi <https://pypi.org/project/spacevlbi/>`_. Ensure that pip is installed and run the following command:

.. code-block:: bash

    pip install spacevlbi

Installing with pip will install/update all of the required libraries automatically (`numpy <http://www.numpy.org/>`_, `poliastro <https://www.poliastro.space/>`_, `matplotlib <http://www.matplotlib.org/>`_, `astropy <http://www.astropy.org/>`).

Documentation
-------------
More detailed documentation for the package is available `here <https://spacevlbi.readthedocs.io/en/latest/>`_.

License
-------
ehtim is licensed under GPLv3. See LICENSE.txt for more details.