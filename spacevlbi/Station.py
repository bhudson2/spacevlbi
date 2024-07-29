# Station.py
#
# This file contains the definitions of three classes used throughout spacevlbi:
# SpaceTelescope, GroundTelescope and GroundStation. Also included are classes
# for defining spacecraft components that can be modelled within the tool.
#
# @author: BenHudson - 28/07/2024

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
    """Class to model a space telescope within the spacevlbi package.

    :param name: Space telescope name, defaults to None
    :type name: str 
    :param mass: Spacecraft mass in kg, defaults to None
    :type mass: float
    :param areaDrag: Spacecraft drag area in m^2, defaults to None
    :type areaDrag: float
    :param areaSolar: Spacecraft solar area in m^2, defaults to None
    :type areaSOlar: float
    :param cD: Drag coefficient, defaults to None
    :type cD: float
    :param cR: Solar reflectivity, defaults to None
    :type cR: float
    :param initTime: Simulation start time, defaults to None
    :type initTime: str
    :param a: Orbit semi-major axis in km, defaults to None
    :type a: float
    :param ecc: Orbit eccentricity, defaults to None
    :type ecc: float
    :param inc: Orbit inclination in degrees, defaults to None
    :type inc: float
    :param aop: Orbit argument of perigee in degrees, defaults to None
    :type aop: float
    :param ta: Orbit true anomaly in degrees, defaults to None
    :type ta: float
    :param pointingVector: Body-fixed vector to point at target source 
        (should be same as RadioPayload.antennaBoresight), defaults to None
    :type pointingVector: numpy.ndarray
    :param constraintVector: Body-fixed vector to point in direction 
        perpendicular to target source. Must be perpendicular to pointingVector,
        defaults to None
    :type constraintVector: numpy.ndarray
    :param rollAngle: Roll angle of space telescope in degrees about the 
        pointingVector. Measured from the plane containing the celestial 
        north pole and the target source direction. Clockwise direction is 
        positive when viewing along the pointingVector direction. Roll 
        angle control is achieved by defining a list of the following 
        format: [transition 1 time in sec, Roll angle 1 in degrees, 
        transition 2 time in sec, roll angle 2, etc.], defaults to None
    :type rollAngle: list
    :param radioPayloads: Array of RadioPayload objects, defaults to None
    :type radioPayloads: list
    :param starTrackers: Array of StarTracker objects, defaults to None
    :type starTrackers: list
    :param reqStarTrackers: Number of star trackers that must not be
        blinded to conduct an observation, defaults to None
    :type reqStarTrackers: int
    :param radiators: Array of Radiator objects, defaults to None
    :type radiators: list
    :param commsSystems: Array of CommsSystem objects, defaults to None
    :type commsSystems: list
    :param solarPanels: Array of SolarPanel objects, defaults to None
    :type solarPanels: list
    :param strModel: Flag indicating whether StarTrackers should be
        modelled, defaults to None
    :type strModel: bool
    :param radModel: Flag indicating whether Radiators should be
        modelled, defaults to None
    :type radModel: bool
    :param commsModel: Flag indicating whether CommsSystems should be
        modelled, defaults to None
    :type commsModel: bool
    :param panelModel: Flag indicating whether SolarPanels should be
        modelled, defaults to None
    :type panelModel: bool
    :return: SpaceTelescope object
    :rtype: SpaceTelescope
    """
        
    def __init__(self, name, mass, areaDrag, areaSolar, cD, cR, initTime, a, \
                 ecc, inc, ra, aop, ta,  pointingVector, constraintVector, \
                 rollAngle, radioPayloads, starTrackers, reqStarTrackers, \
                 radiators, commsSystems, solarPanels, strModel, radModel, 
                 commsModel, panelModel):
               
        self.name = name
        self.pointingVector = np.array([pointingVector])
        self.constraintVector = np.array([constraintVector])
        self.rollAngle = rollAngle
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
        self.reqStarTrackers = reqStarTrackers
        self.panelModel = panelModel
        self.solarPanels = solarPanels
        self.radModel = radModel
        self.radiators = radiators
        self.commsModel = commsModel
        self.commsSystems = commsSystems
        self.sunInertial = np.array([[0,0,0]])  # Sun vector in inertial frame
        self.moonInertial = np.array([[0,0,0]])  # Moon vector in inertial frame
        self.sourceInertial = np.array([[0,0,0]])  # Earth vector in inertial frame
        self.rSun = np.array([[0,0,0]])  # Sat - Sun vector
        self.rMoon = np.array([[0,0,0]])  # Sat - Moon vector
        self.sunBody = np.array([[0,0,0]])  # Sun vector in body-fixed frame
        self.moonBody = np.array([[0,0,0]])  # Moon vector in body-fixed frame
        self.earthBody = np.array([[0,0,0]])  # Earth vector in body-fixed frame
        self.baselines = []  # VLBI baselines
        self.lostBaselines = []  # VLBI baselines lost due to functional constraints
        # Flag indicating whether the source is in view of the spacecraft
        self.sourceVisibility = np.array([0])
        
        
###############################################################################
#   GroundTelescope
###############################################################################

class GroundTelescope:
    """Class to model a ground telescope within the spacevlbi package.
    
    :param name: Ground telescope name, defaults to None
    :type name: str 
    :param diameter: Antenna diameter in metres, defaults to None
    :type diameter: float
    :param apertureEff: Aperture efficiency due to surface design and losses
        from surface deformations, defaults to None
    :type apertureEff: float
    :param sysTemp: System temperature in Kelvin, defaults to None
    :type sysTemp: float
    :param ecefPosition: Ground telescope position in Earth Centered, Earth
        Fixed (ECEF) frame in km, defaults to None
    :type ecefPosition: numpy.ndarray
    :param minEl: Minimum elevation at which observation can be performed
        in degrees, defaults to None
    :type minEl: float
    :param initTime: Simulation start time, defaults to None
    :type initTime: str
    :return: GroundTelescope object
    :rtype: GroundTelescope
    """
    
    def __init__(self, name, diameter, apertureEff, sysTemp, \
               ecefPosition, minEl, initTime):

        self.name = name
        self.diameter = diameter
        self.sysTemp = sysTemp
        self.apertureEff = apertureEff
        self.SEFD = (2*const.k_B*sysTemp)/(self.apertureEff * np.pi *\
            (diameter/2)**2)/(10**-26);  # Jy
        self.ecefPosition = np.array([ecefPosition*1000 << u.m])
        self.minEl = minEl
        # Convert ECEF to ECI positions
        ecef =SkyCoord(x=ecefPosition[0],y=ecefPosition[1],z=ecefPosition[2],\
                       unit='m',frame='itrs',representation_type='cartesian',\
                           obstime = initTime)
        eci = (ecef.transform_to(GCRS)).cartesian       
        self.eciPosition = np.array([[eci.x.value, eci.y.value, eci.z.value]])\
                       << u.m
        self.baselines = []  # VLBI baselines
        self.lostBaselines = []  # VLBI baselines lost due to functional constraints
        # Elevation of source from ground telescope in topocentric frame
        self.elevation = np.array([0])
        self.elevationFlag = np.array([0])


###############################################################################
#   GroundStation
###############################################################################

class GroundStation:
    """Class to model a ground station within the spacevlbi package.
    
    :param name: Ground telescope name, defaults to None
    :type name: str 
    :param ecefPosition: Ground station position in Earth Centered, Earth
        Fixed (ECEF) frame in km, defaults to None
    :type ecefPosition: numpy.ndarray
    :param minEl: Minimum elevation at which link with spacecraft can be 
        achieved in degrees, defaults to None
    :type minEl: float
    :param initTime: Simulation start time, defaults to None
    :type initTime: str
    :return: GroundStation object
    :rtype: GroundStation
    """
    
    def __init__(self, name, ecefPosition, minEl, initTime):
        
        self.name = name
        self.ecefPosition = ecefPosition*1000
        self.minEl = minEl
        # Convert ECEF to ECI positions
        ecef =SkyCoord(x=ecefPosition[0],y=ecefPosition[1],z=ecefPosition[2],\
                       unit='m',frame='itrs',representation_type='cartesian',\
                           obstime = initTime)
        eci = (ecef.transform_to(GCRS)).cartesian       
        self.eciPosition = np.array([eci.x.value, eci.y.value, eci.z.value])*1000 \
                        << u.m
        # Range of spacecraft from ground station
        self.satRange = np.array([0]) << u.m
        # Elevation of spacecraft from ground station
        self.satElev = np.array([0])
        
        
###############################################################################
#   Spacecraft Components
###############################################################################

class RadioPayload:
    """Class to model a radio astronomy payload onboard a SpaceTelescope.
    
    :param name: Radio payload name, defaults to None
    :type name: str 
    :param diameter: Antenna diameter in metres, defaults to None
    :type diameter: float
    :param apertureEff: Aperture efficiency due to surface design and losses
        from surface deformations, defaults to None
    :type apertureEff: float
    :param sysTemp: System temperature in Kelvin, defaults to None
    :type sysTemp: float        
    :param antBoresight: Antenna boresight in body-fixed frame, defaults to
        None
    :type antBoresight: numpy.ndarray
    :param antSunExcl: Exclusion angle between antenna boresight and the 
        Sun's limb in degrees, defaults to None
    :type antSunExcl: float
    :param antEarthExcl: Exclusion angle between antenna boresight and the 
        Earth's limb in degrees, defaults to None
    :type antEarthExcl: float
    :param antMoonExcl: Exclusion angle between antenna boresight and the 
        Moon's limb in degrees, defaults to None
    :type antMoonExcl: float
    :return: RadioPayload object
    :rtype: RadioPayload
    """
    
    def __init__(self, name, diameter, apertureEff, sysTemp,\
                 antBoresight, antSunExcl, antEarthExcl, \
                 antMoonExcl):
 
        self.name = name
        self.diameter = diameter
        self.sysTemp = sysTemp
        self.apertureEff = apertureEff
        self.SEFD = (2*const.k_B*sysTemp)/(self.apertureEff * np.pi *\
            (diameter/2)**2)/(10**-26)  # Jy
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
    """Class to model a star tracker onboard a SpaceTelescope.
    
    :param name: Star tracker name, defaults to None
    :type name: str       
    :param strBoresight: Star tracker boresight in body-fixed frame, 
        defaults to None
    :type strBoresight: numpy.ndarray
    :param strSunExcl: Exclusion angle between star tracker boresight and 
        the Sun's limb in degrees, defaults to None
    :type strSunExcl: float
    :param strEarthExcl: Exclusion angle between star tracker boresight and
        the Earth's limb in degrees, defaults to None
    :type strEarthExcl: float
    :param strMoonExcl: Exclusion angle between star tracker boresight and 
        the Moon's limb in degrees, defaults to None
    :type strMoonExcl: float
    :return: StarTracker object
    :rtype: StarTracker
    """
    
    def __init__(self, name, strBoresight, strSunExcl, strEarthExcl, \
                 strMoonExcl):
        
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
    """Class to model a radiator onboard a SpaceTelescope.
    
    :param name: Radiator name, defaults to None
    :type name: str       
    :param radNorm: Radiator normal vector in body-fixed frame, 
        defaults to None
    :type radNorm: numpy.ndarray
    :param radSunExcl: Exclusion angle between radiator normal and 
        the Sun's limb in degrees, defaults to None
    :type radSunExcl: float
    :param radEarthExcl: Exclusion angle between radiator normal and
        the Earth's limb in degrees, defaults to None
    :type radEarthExcl: float
    :param radMoonExcl: Exclusion angle between radiator normal and 
        the Moon's limb in degrees, defaults to None
    :type radMoonExcl: float
    :return: Radiator object
    :rtype: Radiator
    """
    
    def __init__(self, name, radNorm, radSunExcl, radEarthExcl, \
                 radMoonExcl):
        
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
    """Class to model a solar panel onboard a SpaceTelescope.
    
    :param name: Solar panel name, defaults to None
    :type name: str       
    :param panelNorm: Solar panel normal vector in body-fixed frame, 
        defaults to None
    :type panelNorm: numpy.ndarray
    :return: SolarPanel object
    :rtype: SolarPanel
    """
    
    def __init__(self, name, panelNorm):
        
        self.name = name
        self.panelNorm = panelNorm
        # Panel normal in inertial attitude system
        self.solarPanelInertial = np.array([[0,0,0]])
        # Angle between panel normal and solar direction
        self.betaAngle = np.array([0])
        
        
class CommsSystem:
    """Class to model a communications system onboard a SpaceTelescope.
    
    :param name: Comms system name, defaults to None
    :type name: str       
    :param commsNorm: Comms system normal vector in body-fixed frame, 
        defaults to None
    :param commsFov: Defines the field of view of the comms system within 
        which, if a ground station is insight, a link session can be
        maintained. Could be the half angle beamwidth of a radio or the gimbal
        limit of an optical terminal, etc., defaults to None
    :type commsFov: float
    :param groundReqObs: Ground station contact required for observations? 
        Defaults to None
    :type groundReqObs: bool
    :return: CommsSystem object
    :rtype: CommsSystem
    """
    
    def __init__(self, name, commsNorm, commsFov, groundReqObs):
        
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
        

        

