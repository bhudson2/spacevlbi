"""
Library for simulating space-based VLBI missions (spacevlbi)

Copyright 2024 Ben Hudson

Licensed under the GNU General Public License, Version 3.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.gnu.org/licenses/gpl-3.0.en.html

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# ExampleSetup.py
#
# An example implementation of the spacevlbi package to demonstrate how the
# functions can be executed to perform simulation of a space-based VLBI mission.
# The example given is a simplified representation of the Black Hole Explorer
# (BHEX) mission - https://www.blackholeexplorer.org/

# @author: BenHudson - 22/09/2024

from spacevlbi import Station
from spacevlbi.TimeLoop import TimeLoop
from spacevlbi import Figures
from spacevlbi.Optimisation import Optimisation
import numpy as np
from poliastro.util import Time
from ExampleSpaceTelescope import BaselineBHEX
import matplotlib

# Plot configurations
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams.update({'font.size': 20})

###############################################################################
# Simulation Time Definition
###############################################################################

initTime = Time("2025-01-01 00:00", scale="utc")  # start time of simulation
timeStep = 100  # simulation time step, sec
simLength = 86400  # length of simulation, sec

###############################################################################
# Observation Parameters
###############################################################################

obsFreq = 320e9  # frequency observations will be conducted at, Hz
# Calculate (u,v) coverage for full celestial sphere? NOTE. Functional
# constraints cannot be modelled in all-sky mode.
allsky = 0

# M87*
sourceRa = 187.705930  # target source right ascension, deg
sourceDec = 12.391123  # target source declination, deg

# Sgr A*
#sourceRa = 266.25  # target source right ascension, deg
#sourceDec = -29.0078 # target source declination, deg

intTime = 0  # integration time, sec
dutyCycle = 0  # time between the start of one scan and the next, sec
bandwidth = 32e9  # bandwidth of observations, sec

###############################################################################
# Ground Telescope Definition
###############################################################################

# See GroundTelescope class for definition of input parameters
elev = 15
ALMA = Station.GroundTelescope("ALMA", 73, 0.68, 76, \
                np.array([-23.02920517, -67.75475214, 5.07441706]), elev, initTime)
IRAM = Station.GroundTelescope("IRAM", 30, 0.47, 226, \
                np.array([37.06613845, -3.39260427, 2.92014542]), elev, initTime)
APEX = Station.GroundTelescope("APEX", 12, 0.61, 118, \
                np.array([-23.00577926, -67.75914002, 5.10474176]), elev, initTime)
JCMT = Station.GroundTelescope("JCMT", 15, 0.52, 345, \
                np.array([  19.82283802, -155.47702794, 4.12030997]), elev, initTime)
LMT = Station.GroundTelescope("LMT", 32.5, 0.28, 371, \
                np.array([ 18.98577444, -97.31477955, 4.59354009]), elev, initTime)
SMA = Station.GroundTelescope("SMA", 14.7, 0.75, 285, \
                np.array([  19.82422848, -155.47754761, 4.11529924]), elev, initTime)
SMT = Station.GroundTelescope("SMT", 10, 0.60, 291, \
                np.array([  32.70161115, -109.89124471, 3.15926115]), elev, initTime)
SPT = Station.GroundTelescope("SPT", 6, 0.60, 118, \
                np.array([-90,0,0]), elev, initTime)
NOEMA = Station.GroundTelescope("NOEMA", 52,   0.50, 270,\
                np.array([44.63495268,  5.90666814, 2.7590566 ]), elev, initTime)
HAY = Station.GroundTelescope("HAY", 52, 0.50, 270, \
                np.array([ 42.62394744, -71.48877143, 0.11466984]), elev, initTime)
PV = Station.GroundTelescope("PV", 15, 0.47, 226, \
                np.array([37.06613701, -3.39260427, 2.92002487]), elev, initTime)

###############################################################################
# Spacecraft Definition
###############################################################################

# See ExampleSpaceTelescope.py script for spacecraft definition and the 
# SpaceTelescope class for definition of input parameters
sc1 = BaselineBHEX(initTime)

###############################################################################
# Ground Station Definition
###############################################################################

# See GroundStation class for definition of input parameters
Haleakala = Station.GroundStation("Haleakala", \
                np.array([20.42,-156.15,3.052]), 20, initTime);
Lasilla = Station.GroundStation("La Silla", \
                np.array([-29.15,-70.44,2.4]), 20, initTime);
Achaea = Station.GroundStation("Achaea", \
                np.array([37.48, 22.42, 0]), 20, initTime);
Perth = Station.GroundStation("Perth", \
                np.array([-31.57, 115.51, 0]), 20, initTime);
    
###############################################################################
# Initialise Station Arrays
###############################################################################

spaceTelescopes = [sc1]
groundTelescopes = [ALMA, LMT, SMA, APEX, JCMT, SPT, SMT, NOEMA, HAY, IRAM]
groundStations = [Haleakala, Lasilla, Achaea, Perth]

###############################################################################
# Simulation
###############################################################################

# Execute time loop (i.e. run simulation)
spaceTelescopes, groundTelescopes, groundStations, simTime = TimeLoop(initTime, \
                simLength, timeStep, spaceTelescopes, groundTelescopes, \
                groundStations, obsFreq, sourceRa, sourceDec, dutyCycle, \
                intTime, allsky) 
    
###############################################################################
# Plot Figures
###############################################################################

# Plot space telescope orbits
Figures.OrbitPlot(spaceTelescopes)

# Plot (u,v) coverage
Figures.UvPlot(spaceTelescopes, groundTelescopes, allsky, 1)

# Plot attitude sphere (user can control which elements (e.g. Earth, Sun, Moon,
# antenna, etc.) are included in plot with additional arguments)
Figures.AttitudeSphere(spaceTelescopes, 0, 120, 30)

# Plot incidence angle of Sun on solar panels
Figures.SolarPanelIncidence(spaceTelescopes, simTime, 0)

# Plot elevation of space telescope from ground station
Figures.GroundStationElevation(spaceTelescopes, groundStations, simTime, 0)

# Plot Earth map with space telescope ground track and ground station visibility
Figures.ElevationEarthMap(spaceTelescopes, groundStations, 0)

###############################################################################
# Optimisation
###############################################################################

# Run Optimisation function to determine optimal position for a star tracker
# with a Sun and Earth minimum exclusion angle of 30 degrees.
results = Optimisation(spaceTelescopes, 0, 30, 30, 0)
