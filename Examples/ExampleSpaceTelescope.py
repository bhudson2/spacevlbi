# -*- coding: utf-8 -*-
#
# ExampleSpaceTelescope.py
#
# An example instance of the SpaceTelescope class. The modelled spacecraft is
# representative of the preliminary concept of operations of the Black Hole 
# Explorer (BHEX). More detailed elements of the design defined here are for
# example only and are not intended to accurately represent BHEX's final design
# - https://www.blackholeexplorer.org/
#
# @author: BenHudson - 05/07/2024

from spacevlbi import Station
import numpy as np

def BaselineBHEX(initTime):
    """Example SpaceTelescope object implementation for the Black Hole Explorer
    (BHEX).

    :param initTime: Simulation start datetime, defaults to None
    :type initTime: str    
    :return: sc: SpaceTelescope object representing BHEX
    :rtype sc: SpaceTelescope
    """
    
###############################################################################
#   General spacecraft properties
###############################################################################

    name = "BHEX";
    mass = 300;  # Spacecraft mass, kg
    areaDrag = 4;  # Spacecraft drag area, m^2
    areaSolar = 4;  # Spacecraft solar area, m^2
    cD = 2.2;  # Spacecraft drag coefficient
    cR = 2.2;  # Spacecraft solar reflectivity coefficient
    
###############################################################################
#   Orbit Configuration
###############################################################################

    sma = 26563.88  # semi-major axis, km
    ecc = 0.001  # Eccentricity
    inc = 90  # Inclination, deg
    ra = 247.7  # Right ascension of the ascending node, deg
    aop = 0  # Argument of perigee in degrees
    ta = 90  #  Starting true anomaly in degrees
    
###############################################################################
#   Attitude configuration
###############################################################################

    # Body-fixed axis to be pointed at target source. Should be the same as 
    # RadioPayload antennaBoresight to point antenna at target
    pointingVector = np.array([0,0,1]) 
    # Constraint vector in body frame. Must be perpendicular to pointingVector.
    constraintVector = np.array([0,-1,0])

###############################################################################
#   Payload configuration
###############################################################################

    payloadName = "VLBI"
    antennaDiameter = 3.5  # antenna diameter, metres
    antennaBoresight = np.array([0,0,1])  # antenna boresight in body frame
    antennaEff = 1  # Antena efficiency
    antSunExcl = 90   # Antenna - Sun exclusion angle, deg
    antEarthExcl = 5  # Antenna - Earth limb exclusion angle, deg
    antMoonExcl = 0  # Antenna - Moon exclusion angle, deg
    sysTemp = 4  # System noise temperature, Kelvin
    corrEff = 1  # Correlator efficiency
    clockEff = 1  # Clock efficiency
    
    # Initialise RadioPayload object
    payload1 = Station.RadioPayload(payloadName, antennaDiameter, antennaEff, \
                 corrEff, clockEff, sysTemp, antennaBoresight, antSunExcl, \
                 antEarthExcl, antMoonExcl)
    radioPayloads = [payload1]

###############################################################################
#   Equipment Model configuration
###############################################################################
    
    strModel = 0  # Model star trackers?
    radModel = 0  # Model radiators?
    panelModel = 1  # Model solar panels?
    commsModel = 0  # Model comms systems?
    
###############################################################################
#   Star tracker configuration
###############################################################################

    # Note. During the simulation, if a star tracker is blinded at the current
    # time step, an observation cannot be performed.
    strName1 = "STR1"
    strBoresight1 = np.array([0,0,-1])  # Star tracker boresight in body frame
    strSunExcl1 = 30   # Star tracker - Sun exclusion angle, deg
    strEarthExcl1 = 30  # Star tracker - Earth limb exclusion angle, deg
    strMoonExcl1 = 0  # Star tracker - Moon exclusion angle, deg
    
    # Initialise StarTracker object
    starTracker1 = Station.StarTracker(strName1, strBoresight1, strSunExcl1, \
                                      strEarthExcl1, strMoonExcl1)
    
    strName2 = "STR2"
    strBoresight2 = np.array([0,0,-1])  # Star tracker boresight in body frame
    strSunExcl2 = 30   # Star tracker - Sun exclusion angle, deg
    strEarthExcl2 = 30  # Star tracker - Earth limb exclusion angle, deg
    strMoonExcl2 = 0  # Star tracker - Moon exclusion angle, deg
    
    # Initialise StarTracker object
    starTracker2 = Station.StarTracker(strName2, strBoresight2, strSunExcl2, \
                                       strEarthExcl2, strMoonExcl2)
    
    starTrackers = [starTracker1, starTracker2]
    
###############################################################################
#   Solar panel configuration
###############################################################################
    
    panelName1 = "Solar Panel 1"
    panelNorm1 = np.array([0,0,-1])  # Normal vector in body frame
    # Initialise SolarPanel object
    panel1 = Station.SolarPanel(panelName1, panelNorm1)
                 
    panelName2 = "Solar Panel 2"
    panelNorm2 = np.array([0,1,0])  # Normal vector in body frame
    # Initialise SolarPanel object
    panel2 = Station.SolarPanel(panelName2, panelNorm2)
    
    solarPanels = [panel1, panel2]
    
###############################################################################
#   Radiator surface configuration
###############################################################################

    # Note. During the simulation, if a radiator is blinded at the current
    # time step, an observation cannot be performed.
    radName = "Rad1"
    radNorm = np.array([0.707,0,-0.707]);  # Normal vector in body frame
    radSunExcl = 30  # Radiator - Sun exclusion angle, deg
    radEarthExcl = 30  # Radiator - Earth exclusion angle, deg
    radMoonExcl = 30  # Radiator - Moon exclusion angle, deg
    # Initialise Radiator object
    rad1 = Station.Radiator(radName, radNorm, radSunExcl, radEarthExcl, \
                 radMoonExcl)
        
    radiators = [rad1]
    
###############################################################################
#   Comms configuration
###############################################################################
    
    commsName = "Optical Terminal"
    commsNorm = np.array([0,0,-1])  # Normal vector in body frame
    commsFov = 88  # Beamwidth/gimbal limit (half angle from normal vector), deg
    groundReqObs = 0  # Is a ground station required insight for observations?
    # Initialise CommsSystem object
    comms1 = Station.CommsSystem(commsName,  commsNorm, commsFov, groundReqObs)
    
    commsSystems = [comms1]
    
###############################################################################
#   Generate SpaceTelescope
###############################################################################

    sc = Station.SpaceTelescope(name, mass, areaDrag, areaSolar, cD, cR, \
                 initTime, sma, ecc, inc, ra, aop, ta,  pointingVector, \
                 constraintVector, radioPayloads, starTrackers, radiators, \
                 commsSystems, solarPanels, strModel, radModel, commsModel, \
                 panelModel)
    return sc