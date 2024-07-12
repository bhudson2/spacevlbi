# -*- coding: utf-8 -*-
#
# ExampleSetup.py
#
# An example implementation of the spacevlbi package to demonstrate how the
# functions can be executed to perform simulation of a space-based VLBI mission.
# The example given is a simplified representation of the Black Hole Explorer
# (BHEX) - https://www.blackholeexplorer.org/

# @author: BenHudson - 05/07/2024

from spacevlbi import Station
from spacevlbi.TimeLoop import TimeLoop
from spacevlbi import Figures
import numpy as np
from poliastro.util import Time
from ExampleSpaceTelescope import BaselineBHEX

###############################################################################
# Simulation Time Definition
###############################################################################

initTime = Time("2022-03-11 05:05", scale="utc")  # start time of simulation
timeStep = 1000  # simulation time step, sec
simLength = 86400  # length of simulation, sec

###############################################################################
# Observation Parameters
###############################################################################

obsFreq = 320e9  # frequency observations will be conducted at, Hz
# Calculate (u,v) coverage for full celestial sphere? NOTE. Functional
# constraints cannot be modelled in all-sky mode.
allsky = 0
sourceRa = 187.705930  # target source right ascension, deg
sourceDec = 12.391123  # target source declination, deg
intTime = 300  # integration time, sec
dutyCycle = 0  # time between the start of one integration and the next, sec
bandwidth = 32e9  # bandwidth of observations, sec

###############################################################################
# Ground Telescope Definition
###############################################################################

# See GroundTelescope class for definition of input parameters
ALMA = Station.GroundTelescope("ALMA", 73, 0.68, 1, 1, 76, \
                np.array([2225.0613,-5440.0617,-2481.6812]), 15, initTime)
IRAM = Station.GroundTelescope("IRAM", 30, 0.47, 1, 1, 226, \
                np.array([5088.9678,-301.6812,3825.0122]), 15, initTime)
APEX = Station.GroundTelescope("APEX", 12, 0.61, 1, 1, 118, \
                np.array([2225.0395, -5441.1976, -2479.3034]), 15, initTime)
JCMT = Station.GroundTelescope("JCMT", 15, 0.52, 1, 1, 345, \
                np.array([-5464.5847, -2493.0012, 2150.6540]), 15, initTime)
LMT = Station.GroundTelescope("LMT", 32.5, 0.28, 1, 1, 371, \
                np.array([-768.7156, -5988.5071, 2063.3549]), 15, initTime)
SMA = Station.GroundTelescope("SMA", 14.7, 0.75, 1, 1, 285, \
                np.array([-5464.5555, -2492.9280, 2150.7972]), 15, initTime)
SMT = Station.GroundTelescope("SMT", 10, 0.60, 1, 1, 291, \
                np.array([-1828.7962, -5054.4068, 3427.8652]), 15, initTime)
SPT = Station.GroundTelescope("SPT", 6, 0.60, 1, 1, 118, \
                np.array([0.8098, -0.8169, -6359.5687]), 15, initTime)
NOEMA = Station.GroundTelescope("NOEMA", 52,   0.50, 1, 1, 270,\
                np.array([4524.0004, 468.0421, 4460.5098]), 15, initTime)
HAY = Station.GroundTelescope("HAY", 52, 0.50, 1, 1, 270, \
                np.array([1492.341, -4457.234,  4296.933]), 15, initTime)

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
                np.array([1285.03,3443.22, 6216.94]), 5, initTime);
    
###############################################################################
# Initialise Station Arrays
###############################################################################

spaceTelescopes = [sc1]
groundTelescopes = [ALMA, IRAM, NOEMA, LMT, SMA, HAY, APEX, JCMT, SPT, SMT]
groundStations = [Svalbard]

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
Figures.AttitudeSphere(spaceTelescopes, 0, 45, 30)

# Plot incidence angle of Sun on solar panels
Figures.SolarPanelIncidence(spaceTelescopes, 0, simTime)

# Plot elevation of space telescope from ground station
Figures.GroundStationElevation(spaceTelescopes, 0, groundStations, simTime)