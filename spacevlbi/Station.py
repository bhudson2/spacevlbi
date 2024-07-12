# -*- coding: utf-8 -*-
#
# Station.py
#
# This file contains the definitions of three classes used throughout spacevlbi:
# SpaceTelescope, GroundTelescope and GroundStation. Also included are classes
# for defining spacecraft components that can be modelled within the tool.
#
# @author: BenHudson - 09/07/2024

from astropy.coordinates import SkyCoord,GCRS,ITRS
from astropy import constants as const
import numpy as np
from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.bodies import Earth

###############################################################################
#   SpaceTelescope
###############################################################################

class SpaceTelescope:
    
    def __init__(self, name, mass, areaDrag, areaSolar, cD, cR, initTime, a, \
                 ecc, inc, ra, aop, ta,  pointingVector, constraintVector, \
                 radioPayloads, starTrackers, radiators, commsSystems, \
                 solarPanels, strModel, radModel, commsModel, panelModel):
        """SpaceTelescope object initialisation function.

           Args:
                name (str): Space telescope name
                mass (float): Spacecraft mass, kg
                areaDrag (float): Spacecraft drag area, m**2
                areaSolar (float): Spacecraft solar area, m**2
                cD (float): Drag coefficient
                cR (float): Solar reflectivity
                initTime (Time): Simulation start time
                a (float): Semi-major axis, km
                ecc (float): Eccentricity
                inc (float): Inclination, deg
                aop (float): Argument of Perigee, deg
                ta (float): Initial true anomaly, deg
                pointingVector (float): Body-fixed vector to point at target 
                                        source (should be same as RadioPayload
                                        antennaBoresight)
                constraintVector (float): Body-fixed vector to point in
                                          direction perpendicular to target
                                          source
                radioPayloads (obj): Array of RadioPayload objects
                starTrackers (obj): Array of StarTracker objects
                radiators (obj): Array of Radiators objects
                commsSystems (obj): Array of CommsSystem objects
                solarPanels (obj): Array of SolarPanel objects
                strModel (BOOL): Flag indicating whether StarTrackers should be
                                 modelled
                radModel (BOOL): Flag indicating whether Radiators should be
                                 modelled
                commsModel (BOOL): Flag indicating whether CommsSystems should be
                                 modelled
                panelModel (BOOL): Flag indicating whether SolarPanels should be
                                 modelled

           Returns:
               SpaceTelescope (obj): Instance of SpaceTelescope
        """   
        self.name = name
        self.pointingVector = np.array([pointingVector])
        self.constraintVector = np.array([constraintVector])
        self.a = a << u.km
        self.ecc = ecc << u.one
        self.inc = inc << u.deg
        self.ra = ra << u.deg
        self.aop = aop << u.deg
        self.ta = ta << u.deg
        self.orbit = Orbit.from_classical(Earth, self.a, self.ecc, self.inc, \
                      self.ra, self.aop, self.ta, initTime)
        self.eciPosition = np.array([self.orbit.r*1000])<< u.m
        self.eciVelocity = np.array([self.orbit.v*1000])<< (u.m / u.s)
        eci =SkyCoord(x=self.eciPosition[0,0],y=self.eciPosition[0,1],\
                      z=self.eciPosition[0,2], unit='m',frame='gcrs',\
                      representation_type='cartesian', obstime = initTime)
        ecef = (eci.transform_to(ITRS)).cartesian
        self.ecefPosition = np.array([[ecef.x.value, ecef.y.value, \
                      ecef.z.value]]) << u.m
        self.earthLimbAngle = np.array([0])
        self.sunLimbAngle = np.array([0])
        self.moonLimbAngle = np.array([0])
        self.eclipse = np.array([0])
        self.mass = mass
        self.areaDrag = areaDrag
        self.areaSolar = areaSolar
        self.cD = cD
        self.cR = cR
        self.radioPayloads = radioPayloads
        self.strModel = strModel
        self.starTrackers = starTrackers
        self.panelModel = panelModel
        self.solarPanels = solarPanels
        self.radModel = radModel
        self.radiators = radiators
        self.commsModel = commsModel
        self.commsSystems = commsSystems
        self.sunInertial = np.array([[0,0,0]])
        self.moonInertial = np.array([[0,0,0]])
        self.sourceInertial = np.array([[0,0,0]])
        self.sunBody = np.array([[0,0,0]])
        self.moonBody = np.array([[0,0,0]])
        self.earthBody = np.array([[0,0,0]])
        self.baselines = []
        self.lostBaselines = []
        self.sourceVisibility = np.array([0])
        self.simTime = 0
        
###############################################################################
#   GroundTelescope
###############################################################################

class GroundTelescope:
    
    def __init__(self, name, diameter, antennaEff, corrEff, clockEff, sysTemp, \
               ecefPosition, minEl, initTime):
        """GroundTelescope object initialisation function.

           Args:
                name (str): Space telescope name
                diameter (float): Antenna diameter, metres
                antennaEff (float): Antenna efficiency
                corrEff (float): Correlator efficiency
                clockEff (float): Clock efficiency
                sysTemp (float): System temperature, Kelvin
                ecefPosition (float): Ground telescope position, metres
                minEl (float): Minimum elevation at which observation can be
                               observation can be performed, deg
                initTime (Time): Simulation start time

           Returns:
               GroundTelescope (obj): Instance of GroundTelescope
        """   
        self.name = name
        self.diameter = diameter
        self.sysTemp = sysTemp
        self.efficiency = antennaEff * corrEff * clockEff
        self.SEFD = (2*const.k_B*sysTemp)/(self.efficiency * np.pi *\
            (diameter/2)**2)/10**-26;  # Jy
        self.ecefPosition = np.array([ecefPosition*1000 << u.m])
        self.minEl = minEl
        # Convert ECEF to ECI positions
        ecef =SkyCoord(x=ecefPosition[0],y=ecefPosition[1],z=ecefPosition[2],\
                       unit='m',frame='itrs',representation_type='cartesian',\
                           obstime = initTime)
        eci = (ecef.transform_to(GCRS)).cartesian       
        self.eciPosition = np.array([[eci.x.value, eci.y.value, eci.z.value]])\
                       << u.m
        self.baselines = []
        self.lostBaselines = []
        self.elevation = np.array([0])
        self.elevationFlag = np.array([0])
        self.simTime = 0

###############################################################################
#   GroundStation
###############################################################################

class GroundStation:
    
    def __init__(self, name, ecefPosition, minEl, initTime):
        """GroundStation object initialisation function.

           Args:
                name (str): Space telescope name
                ecefPosition (float): Ground telescope position, metres
                minEl (float): Minimum elevation at which observation can be
                               observation can be performed, deg
                initTime (Time): Simulation start time

           Returns:
               GroundStation (obj): Instance of GroundStation
        """ 
        self.name = name
        self.ecefPosition = ecefPosition*1000
        self.minEl = minEl
        # Convert ECEF to ECI positions
        ecef =SkyCoord(x=ecefPosition[0],y=ecefPosition[1],z=ecefPosition[2],\
                       unit='m',frame='itrs',representation_type='cartesian',\
                           obstime = initTime)
        eci = (ecef.transform_to(GCRS)).cartesian       
        self.eciPosition = np.array([eci.x.value, eci.y.value, eci.z.value]) \
                        << u.m
        self.satRange = np.array([0]) << u.m
        self.satElev = np.array([0])
        self.simTime = 0
        
###############################################################################
#   Spacecraft Components
###############################################################################

class RadioPayload:
    
    def __init__(self, name, diameter, antennaEff, corrEff, clockEff, sysTemp,\
                 antBoresight, antSunExcl, antEarthExcl, \
                 antMoonExcl):
        """RadioPayload object initialisation function. Used to model simple
           elements of a radio telescope payload. If RadioPayload's antenna's
           Sun, Earth and/or Moon exclusion angles are violated during
           simulation, observation cannot take place.

           Args:
                name (str): Antenna name 
                diameter (float: Antenna diameter
                antennaEff (float): Antenna efficiency
                corrEff (float): Correlator efficiency
                clockEff (float): Clock efficiency
                sysTemp (float): System temperature, Kelvin
                antBoresight (float): Antenna boresight in body-frame
                antSunExcl (float): Exclusion angle between Antenna
                                    boresight and the Sun's limb, deg
                antEarthExcl (float): Exclusion angle between Antenna
                                    boresight and the Earth's limb, deg
                antMoonExcl (float): Exclusion angle between Antenna
                                    boresight and the Moon's limb, deg

           Returns:
               RadioPayload (obj): Instance of Antenna
        """   
        self.name = name
        self.diameter = diameter
        self.sysTemp = sysTemp
        self.efficiency = antennaEff * corrEff * clockEff;
        self.SEFD = (2*const.k_B*sysTemp)/(self.efficiency * np.pi *\
            (diameter/2)**2)/10**-26  # Jy
        self.antBoresight = antBoresight
        self.antSunExcl = antSunExcl
        self.antEarthExcl = antEarthExcl
        self.antMoonExcl = antMoonExcl
        # Antenna normal in inertial attitude system
        self.antennaInertial = np.array([0,0,0])
        # Havw the Sun/Earth/Moon exclusion angles been violated?
        self.sunFlag = np.array([0])
        self.earthFlag = np.array([0])
        self.moonFlag = np.array([0])
        

class StarTracker:
    
    def __init__(self, name, strBoresight, strSunExcl, strEarthExcl, \
                 strMoonExcl):
        """StarTracker object initialisation function. If StarTracker is 
           blinded during simulation, observations cannot take place.

           Args:
                name (str): Star tracker name 
                strBoresight (float): Star tracker boresight in body-frame
                strSunExcl (float): Exclusion angle between star tracker
                                    boresight and the Sun's limb, deg
                strEarthExcl (float): Exclusion angle between star tracker
                                    boresight and the Earth's limb, deg
                strMoonExcl (float): Exclusion angle between star tracker
                                    boresight and the Moon's limb, deg

           Returns:
               StarTracker (obj): Instance of StarTracker
        """   
        self.name = name
        self.strBoresight = strBoresight
        self.strSunExcl = strSunExcl
        self.strEarthExcl = strEarthExcl
        self.strMoonExcl = strMoonExcl
        # Flag indicating whether star tracker is blinded (1==blinded)
        self.strBlindFlag = np.array([0])
        # StarTracker normal in inertial attitude system
        self.strInertial = np.array([[0,0,0]])
        
    
class Radiator:
    
    def __init__(self, name, radNorm, radSunExcl, radEarthExcl, \
                 radMoonExcl):
        """Radiator object initialisation function. If Radiator is 
           blinded during simulation, observations cannot take place.

           Args:
                name (str): Radiator name
                radNorm (float): Radiator normal vector in body-frame
                radSunExcl (float): Exclusion angle between Radiator normal
                                    vector and the Sun's limb, deg
                radEarthExcl (float): Exclusion angle between Radiator normal 
                                    vector and the Earth's limb, deg
                radMoonExcl (float): Exclusion angle between Radiator normal
                                    vector and the Moon's limb, deg

           Returns:
               Radiator (obj): Instance of Radiator
        """   
        self.name = name
        self.radNorm = radNorm
        self.radSunExcl = radSunExcl
        self.radEarthExcl = radEarthExcl
        self.radMoonExcl = radMoonExcl
        # Radiator normal in inertial attitude system
        self.radInertial = np.array([[0,0,0]])
        # Flag indicating whether radiator is blinded (1==blinded)
        self.radBlindFlag = np.array([0])

        
class SolarPanel:
    
    def __init__(self, name, panelNorm):
        """SolarPanel object initialisation function. SolarPanels have no
           impact on the observations.

           Args:
                name (str): Radiator name
                panelNorm (float): Panel normal vector in body-frame

           Returns:
               SolarPanel (obj): Instance of SolarPanel
        """   
        self.name = name
        self.panelNorm = panelNorm
        # Panel normal in inertial attitude system
        self.solarPanelInertial = np.array([[0,0,0]])
        # Angle between panel normal and solar direction
        self.betaAngle = np.array([0])
        
        
class CommsSystem:
    
    def __init__(self, name, commsNorm, commsFov, groundReqObs):
        """CommsSystem object initialisation function. CommsSystem can be used
           to model systems such as optical terminal or flight radios and
           perform analysis on the access times between the spacecraft and the 
           ground.

           Args:
                name (str): CommsSystem name
                commsNorm (float): CommsSystem normal vector in body-frame
                commsFov (float): Defines the field of view of the comms system
                                  within which, if a ground station is insight,
                                  a link session can be maintained. Could be
                                  the half angle beamwidth of a radio or the
                                  gimbal limit of an optical terminal, etc.
                groundReqObs (BOOL): Ground station contact required for 
                                    observations?


           Returns:
                CommsSystem (obj): Instance of CommsSystem
        """   
        self.name = name
        self.commsNorm = commsNorm
        self.commsFov = commsFov
        self.groundReqObs = groundReqObs
        # Comms system normal in inertial attitude system
        self.commsInertial = np.array([[0,0,0]])
        # Angle between comms system normal and ground station direction
        self.commsGroundAngle = np.array([0])
        # Array of flags indicating which ground stations are insight
        self.groundStationInsight = np.array([0])
        

        

