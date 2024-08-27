# ExampleSetup.py
#
# An example implementation of the spacevlbi package to demonstrate how the
# functions can be executed to perform simulation of a space-based VLBI mission.
# The example given is a simplified representation of the Black Hole Explorer
# (BHEX) mission - https://www.blackholeexplorer.org/

# @author: BenHudson - 27/08/2024

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
dutyCycle = 0  # time between the start of one integration and the next, sec
bandwidth = 32e9  # bandwidth of observations, sec

###############################################################################
# Ground Telescope Definition
###############################################################################

# See GroundTelescope class for definition of input parameters
elev = 15
ALMA = Station.GroundTelescope("ALMA", 73, 0.68, 76, \
                np.array([2225.0613,-5440.0617,-2481.6812]), elev, initTime)
IRAM = Station.GroundTelescope("IRAM", 30, 0.47, 226, \
                np.array([5088.9678,-301.6812,3825.0122]), elev, initTime)
APEX = Station.GroundTelescope("APEX", 12, 0.61, 118, \
                np.array([2225.0395, -5441.1976, -2479.3034]), elev, initTime)
JCMT = Station.GroundTelescope("JCMT", 15, 0.52, 345, \
                np.array([-5464.5847, -2493.0012, 2150.6540]), elev, initTime)
LMT = Station.GroundTelescope("LMT", 32.5, 0.28, 371, \
                np.array([-768.7156, -5988.5071, 2063.3549]), elev, initTime)
SMA = Station.GroundTelescope("SMA", 14.7, 0.75, 285, \
                np.array([-5464.5555, -2492.9280, 2150.7972]), elev, initTime)
SMT = Station.GroundTelescope("SMT", 10, 0.60, 291, \
                np.array([-1828.7962, -5054.4068, 3427.8652]), elev, initTime)
SPT = Station.GroundTelescope("SPT", 6, 0.60, 118, \
                np.array([0.8098, -0.8169, -6359.5687]), elev, initTime)
NOEMA = Station.GroundTelescope("NOEMA", 52,   0.50, 270,\
                np.array([4524.0004, 468.0421, 4460.5098]), elev, initTime)
HAY = Station.GroundTelescope("HAY", 52, 0.50, 270, \
                np.array([1492.341, -4457.234,  4296.933]), elev, initTime)
PV = Station.GroundTelescope("PV", 15, 0.47, 226, \
                np.array([5088.9678, -301.6812,  3825.012]), elev, initTime)

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
Svalbard = Station.GroundStation("Svalbard", \
                np.array([1258.4, 346.3, 6222.2]), 5, initTime);
Haleakala = Station.GroundStation("Haleakala", \
                np.array([-5466.003,-2404.290, 2242.294]), 5, initTime);
Lasilla = Station.GroundStation("La Silla", \
                np.array([1838.689,5259.299,3099.28]), 5, initTime);
Nemea = Station.GroundStation("Nemea", \
                np.array([4654.281,1947.909,3888.707]), 5, initTime);
Perth = Station.GroundStation("Perth", \
                np.array([-2384.691,4860.073,-3361.166]), 5, initTime);
    
###############################################################################
# Initialise Station Arrays
###############################################################################

spaceTelescopes = [sc1]
groundTelescopes = [ALMA, LMT, SMA, APEX, JCMT, SPT, SMT, NOEMA, HAY, IRAM]
groundStations = [Haleakala, Lasilla, Nemea, Perth]

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

###############################################################################
# Optimisation
###############################################################################

# Run Optimisation function to determine optimal position for a star tracker
# with a Sun and Earth minimum exclusion angle of 30 degrees.
results = Optimisation(spaceTelescopes, 0, 30, 30, 0)
