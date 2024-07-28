# Constraints.py
#
# Implementation of the functional constraints that apply to the space telescopes
# and the ground telescopes that will impact when observations can be performed
# in spacevlbi.
#
# @author: BenHudson - 28/07/2024

from astropy import constants as const
from astropy import units as u
import numpy as np
from numpy.linalg import norm
from numpy import vstack, dot, degrees, radians, cos, sin, arccos, arctan

###############################################################################
# ObsLimits
###############################################################################

def ObsLimits(spaceTelescopes, groundTelescopes, groundStations, sourceRa, \
              sourceDec, rSun, rMoon):
    """Calculate whether any of the functional constraints are stopping an
    observation from being performed at the current time step.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundTelescopes: Array of GroundTelescope objects, defaults to None
    :type groundTelescopes: list
    :param groundStations: Array of GroundStation objects, defaults to None
    :type groundStations: list
    :param sourceRa: Right ascension of target source in degrees, defaults to None
    :type sourceRa: float
    :param sourceDec: Declination of target source in degrees, defaults to None
    :type sourceDec: float
    :param rSun: Sun position vector in ECI frame in metres, defaults to None
    :type rSun: numpy.ndarray
    :param rMoon: Moon position vector in ECI frame in metres, defaults to None
    :type rMoon: numpy.ndarray
    :return: Array of SpaceTelescope objects
    :rtype: list
    :return: Array of GroundTelescope objects
    :rtype: list
    """
    
    # Iterate through space telescopes and calculate functional constraints
    if spaceTelescopes:
        for j in range(len(spaceTelescopes)):
            position = spaceTelescopes[j].eciPosition[-1,:] << u.m
            # Is the target source in view of the spacecraft?
            vis = SourceVisibility(sourceRa, sourceDec, position)
            spaceTelescopes[j].sourceVisibility = vstack((spaceTelescopes[j].sourceVisibility, \
                                vis))
                
            # Calculate angle of Earth limb from spacecraft position vector
            earthLimb = arctan(const.R_earth / norm(position))
            spaceTelescopes[j].earthLimbAngle = vstack((spaceTelescopes[j].earthLimbAngle, \
                                earthLimb))
                
            # Calculate angle of Sun limb from spacecraft position vector
            satSun = rSun - position
            sunLimb = arctan(const.R_sun / norm(satSun))
            spaceTelescopes[j].sunLimbAngle = vstack((spaceTelescopes[j].sunLimbAngle, \
                                sunLimb))
                
            # Calculate angle of Moon limb from spacecraft position vector
            moonRad = 1737.4e3 << u.m
            satMoon = rMoon - position
            moonLimb = arctan(moonRad / norm(satMoon))
            spaceTelescopes[j].moonLimbAngle = vstack((spaceTelescopes[j].moonLimbAngle, \
                                moonLimb))
                
            # Update Sun and Moon vectors in SpaceTelescope
            spaceTelescopes[j].rSun = vstack((spaceTelescopes[j].rSun, satSun))
            spaceTelescopes[j].rMoon = vstack((spaceTelescopes[j].rMoon, satMoon))
            
            sunInertial = spaceTelescopes[j].sunInertial[-1,:]
            moonInertial = spaceTelescopes[j].moonInertial[-1,:]
            earthInertial = (-spaceTelescopes[j].eciPosition[-1,:]/ \
                                norm(spaceTelescopes[j].eciPosition[-1,:]))
                
            # Is the antenna pointing within Sun, Moon or Earth exclusion zones?
            radioPayloads = spaceTelescopes[j].radioPayloads
            if radioPayloads != []:
                for k in range(len(radioPayloads)):
                    antennaInertial = radioPayloads[k].antennaInertial[-1,:]
                
                    sunAngle = abs(arccos(dot(antennaInertial, sunInertial)/ \
                        (norm(antennaInertial)*norm(sunInertial)))) - sunLimb
                    moonAngle = abs(arccos(dot(antennaInertial, moonInertial)/ \
                        (norm(antennaInertial)*norm(moonInertial)))) - moonLimb
                    earthAngle = abs(arccos(dot(antennaInertial, earthInertial)/ \
                        (norm(antennaInertial)*norm(earthInertial))))-earthLimb
            
                    if degrees(sunAngle.value) < radioPayloads[k].antSunExcl:
                        radioPayloads[k].sunFlag = vstack((radioPayloads[k].sunFlag,0))
                    else:
                        radioPayloads[k].sunFlag = vstack((radioPayloads[k].sunFlag,1))
                    if degrees(moonAngle.value) < radioPayloads[k].antMoonExcl:
                        radioPayloads[k].moonFlag = vstack((radioPayloads[k].moonFlag,0))
                    else:
                        radioPayloads[k].moonFlag= vstack((radioPayloads[k].moonFlag,1))
                    if degrees(earthAngle.value) < radioPayloads[k].antEarthExcl:
                        radioPayloads[k].earthFlag = vstack((radioPayloads[k].earthFlag,0))
                    else:
                        radioPayloads[k].earthFlag = vstack((radioPayloads[k].earthFlag,1))
                else:
                    radioPayloads[k].sunFlag = vstack((radioPayloads[k].sunFlag,1))
                    radioPayloads[k].moonFlag = vstack((radioPayloads[k].moonFlag,1))
                    radioPayloads[k].earthFlag = vstack((radioPayloads[k].earthFlag,1))
                spaceTelescopes[j].radioPayloads = radioPayloads
    
            # Are the star trackers pointing within Sun, Moon or Earth exclusion
            # angles?
            starTrackers = spaceTelescopes[j].starTrackers
            if starTrackers != [] and spaceTelescopes[j].strModel == 1:
                for k in range(len(starTrackers)):
                    strInertial = starTrackers[k].strInertial[-1,:]
                    strSunAngle = abs(arccos(dot(strInertial, sunInertial)/ \
                        (norm(strInertial)*norm(sunInertial))))
                    strMoonAngle = abs(arccos(dot(strInertial, moonInertial)/ \
                        (norm(strInertial)*norm(moonInertial))))
                    strEarthAngle = abs(arccos(dot(strInertial, earthInertial)/ \
                        (norm(strInertial)*norm(earthInertial)))) - earthLimb;
                    if (degrees(strSunAngle).value<starTrackers[k].strSunExcl) or \
                            (degrees(strMoonAngle).value<starTrackers[k].strMoonExcl) or \
                            (degrees(strEarthAngle.value)<starTrackers[k].strEarthExcl):
                       starTrackers[k].strBlindFlag = vstack((starTrackers[k].strBlindFlag,0))
                    else:
                       starTrackers[k].strBlindFlag = vstack((starTrackers[k].strBlindFlag,1))
                starTrackers[k].starTrackers = starTrackers
                    
            # Are the radiators pointing within Sun, Moon or Earth exclusion angles?
            radiators = spaceTelescopes[j].radiators
            if radiators != [] and spaceTelescopes[j].radModel == 1:
                for k in range(len(radiators)):   
                    radInertial = radiators[k].radInertial[-1,:]
                    radSunAngle = abs(arccos(dot(radInertial, sunInertial)/ \
                        (norm(radInertial)*norm(sunInertial))));
                    radMoonAngle = abs(arccos(dot(radInertial, moonInertial)/ \
                        (norm(radInertial)*norm(moonInertial))));
                    radEarthAngle = abs(arccos(dot(radInertial, earthInertial)/ \
                        (norm(radInertial)*norm(earthInertial)))) - earthLimb;
                    if (degrees(radSunAngle.value)<radiators[k].radSunExcl) or \
                            (degrees(radMoonAngle.value)<radiators[k].radMoonExcl) or \
                            (degrees(radEarthAngle.value)<radiators[k].radEarthExcl):
                       radiators[k].radBlindFlag = vstack((radiators[k].radBlindFlag,0))
                    else:
                       radiators[k].radBlindFlag = vstack((radiators[k].radBlindFlag,1))
                spaceTelescopes[j].radiators = radiators
    
            # Are there any ground station limitations on observations?
            commsSystems = spaceTelescopes[j].commsSystems
            if commsSystems != [] and spaceTelescopes[j].commsModel == 1:
                for k in range(len(commsSystems)):
                    # Is live downlink of the data being performed?
                    if commsSystems[k].groundReqObs == 1:
                        # Initialise angle parameters
                        angle = np.zeros((1, len(groundStations)))
                        insight = np.zeros((1, len(groundStations)))
                        # Iterate through ground stations
                        for l in range(len(groundStations)):
                            # Calculate elevation of ground station to spacecraft
                            groundECI = groundStations[l].eciPosition[-1,:]
                            el = groundStations[l].satElev[-1,j]
                            minEl = groundStations[l].minEl
                            
                            # Caculate angle between spacecraft comms system normal
                            # vector and the ground station
                            commsInertial = commsSystems[k].commsInertial[-1,:]
                            spaceGround = groundECI - position
                            spaceGround = spaceGround / norm(spaceGround)
                            commsGroundAngle = degrees(abs(arccos(dot(commsInertial, \
                                        spaceGround)/ (norm(commsInertial)* \
                                        norm(spaceGround))))).value
                            commsFov = commsSystems[k].commsFov
                            angle[0,l] = commsGroundAngle
        
                            # Sspacecraft must be above minEl and the comms terminal 
                            # beamwidth / gimbal limit must be pointed within commFOV
                            # of the ground station for a link to be maintained
                            if (el < 90 and el > minEl) and commsGroundAngle <= commsFov:
                                insight[0,l] = 1;
                            else:
                                insight[0,l] = 0;
                                
                        if len(commsSystems[k].commsGroundAngle) == 1:
                            commsSystems[k].commsGroundAngle = \
                                np.append(commsSystems[k].commsGroundAngle, \
                                        np.zeros((1,len(groundStations)-1)))
                        if len(commsSystems[k].groundStationInsight) == 1:
                            commsSystems[k].groundStationInsight = \
                                np.append(commsSystems[k].groundStationInsight, \
                                        np.zeros((1,len(groundStations)-1)))
                        
                        commsSystems[k].commsGroundAngle = \
                            vstack((commsSystems[k].commsGroundAngle, angle))
                        commsSystems[k].groundStationInsight = \
                            vstack((commsSystems[k].groundStationInsight, insight))
                spaceTelescopes[j].commsSystems = commsSystems
    
    if groundTelescopes:
        for j in range(len(groundTelescopes)):
            # Is the elevation between the ground station and target source
            # within limits?
            position = groundTelescopes[j].eciPosition[-1,:]
            el = Elevation(position, sourceRa, sourceDec);
            groundTelescopes[j].elevation = vstack((groundTelescopes[j].elevation,el));
            minEl = groundTelescopes[j].minEl
            if el < 90 and el > minEl:
                groundTelescopes[j].elevationFlag = vstack((groundTelescopes[j].\
                                                            elevationFlag,1));
            else:
                groundTelescopes[j].elevationFlag = vstack((groundTelescopes[j].\
                                                            elevationFlag,0));

    return spaceTelescopes, groundTelescopes


###############################################################################
# SourceVisibility
###############################################################################

def SourceVisibility(ra, dec, position):
    """Calculate whether source at given right ascension and declination is in
    view of a space telescope (i.e. not blocked by the Earth).

    :param ra: Right ascension of target source in degrees, defaults to None
    :type ra: float
    :param dec: Declination of target source in degrees, defaults to None
    :type dec: float
    :param position: Spacecraft ECI position vector in metres, defaults to None
    :type position: numpy.ndarray
    :return: Flag indicating whether source is in view of spacecraft.
        In View == 1, Obstructed == 0.
    :rtype: bool
    """
    
    R_Earth = const.R_earth
    # Calculate Earth-source unit vector
    sourceECI = np.array([sin(np.pi/2 - radians(dec)) * cos(radians(ra)), \
                          sin(np.pi/2 - radians(dec)) * sin(radians(ra)), \
                          cos(np.pi/2 - radians(dec))])
    earthSource = (sourceECI / np.linalg.norm(sourceECI)) << u.m
    
    if (np.dot(position.value, earthSource.value) < \
        -np.sqrt(np.linalg.norm(position.value)**2 - R_Earth.value**2)):
        visibility = 0;
    else:
        visibility = 1;
    
    return visibility

###############################################################################
# Elevation
###############################################################################

def Elevation(position, ra, dec):
    """Calculate elevation of a target source from the horizon at a ground 
    telescope location in the topocentric frame.

    :param position: GroundTelescope ECI position vector in metres, defaults to 
        None
    :param ra: Right ascension of target source in degrees, defaults to None
    :type ra: float
    :param dec: Declination of target source in degrees, defaults to None
    :type dec: float
    :type position: numpy.ndarray
    :return: Elevation of target source from horizon in degrees
    :rtype: float
    """
    
    # Calculate source vector in ECI
    sourceECI = np.array([cos(radians(dec)) * cos(radians(ra)), \
                          cos(radians(dec)) * sin(radians(ra)), \
                          sin(radians(dec))])
    earthSourceECI = (sourceECI / norm(sourceECI)) * const.kpc

    # Generate vector from ground telescope to source
    groundSource = earthSourceECI - position
    
    # Calculate elevation of source from ground telescope
    elRad = dot(groundSource, position)/(norm(groundSource)* norm(position));
    el = degrees(0.5*np.pi - arccos(elRad).value);
    
    return el

###############################################################################
# ObsMask
###############################################################################

def ObsMask(telescope1, telescope2):
    """For two specific antenna, calculate whether an observation can be
    performed at the current time step.

    :param telescope1: SpaceTelescope or GroundTelescope object, defaults to None
    :type telescope1: SpaceTelescope or GroundTelescope
    :param telescope2: SpaceTelescope or GroundTelescope object, defaults to None
    :type telescope2: SpaceTelescope or GroundTelescope
    :return: Flag indicating whether observation can be performed by 
        these two antenna at the current time step
    :rtype: bool
    """
    
    # Check telescope1 functional constraint flags
    if telescope1.__class__.__name__ == 'SpaceTelescope':
        
        radioPayloads = telescope1.radioPayloads
        starTrackers = telescope1.starTrackers
        radiators = telescope1.radiators
        commsSystems = telescope1.commsSystems
        
        radioFlag = 1
        if radioPayloads != []:    
            for j in range(len(radioPayloads)):
                sunFlag = radioPayloads[j].sunFlag[-1,:]
                earthFlag = radioPayloads[j].earthFlag[-1,:]
                moonFlag = radioPayloads[j].moonFlag[-1,:]
                radioFlag = radioFlag * sunFlag * earthFlag * moonFlag

        strFlag = 1
        strReq = telescope1.reqStarTrackers
        strCount = 0
        if starTrackers != [] and telescope1.strModel == 1:    
            for j in range(len(starTrackers)):
                blindFlag = starTrackers[j].strBlindFlag[-1,:]
                if blindFlag == 0:
                    strCount = strCount + 1
            
            if strCount < strReq:
                strFlag = 0
            else:
                strFlag = 1
        
        radFlag = 1
        if radiators != [] and telescope1.radModel == 1:    
            for j in range(len(radiators)):
                blindFlag = radiators[j].radBlindFlag[-1,:]
                radFlag = radFlag * blindFlag

        # Ground station visibility limits
        commsFlag = 1
        if commsSystems != [] and telescope1.commsModel == 1:    
            for j in range(len(commsSystems)):
                # Is live downlink of the data being performed?
                if commsSystems[j].groundReqObs == 1:
                    flag = any(commsSystems[j].groundStationInsight[-1,:])
                else:
                    flag = 1
                commsFlag = commsFlag * flag

        flag1 = radioFlag * strFlag * radFlag * commsFlag
    else:
        flag1 = 1;
    
    # Check telescope1 functional constraint flags
    if telescope2.__class__.__name__ == 'SpaceTelescope':
        
        radioPayloads = telescope2.radioPayloads
        starTrackers = telescope2.starTrackers
        radiators = telescope2.radiators
        commsSystems = telescope2.commsSystems
        
        radioFlag = 1
        if radioPayloads != []:    
            for j in range(len(radioPayloads)):
                sunFlag = radioPayloads[j].sunFlag[-1,:]
                earthFlag = radioPayloads[j].earthFlag[-1,:]
                moonFlag = radioPayloads[j].moonFlag[-1,:]
                radioFlag = radioFlag * sunFlag * earthFlag * moonFlag

        strFlag = 1
        if starTrackers != [] and telescope2.strModel == 1:    
            for j in range(len(starTrackers)):
                blindFlag = starTrackers[j].strBlindFlag[-1,:]
                strFlag = strFlag * blindFlag
        
        radFlag = 1
        if radiators != [] and telescope2.radModel == 1:    
            for j in range(len(radiators)):
                blindFlag = radiators[j].radBlindFlag[-1,:]
                radFlag = radFlag * blindFlag

        # Ground station visibility limits
        commsFlag = 1
        if commsSystems != [] and telescope2.commsModel == 1:    
            for j in range(len(commsSystems)):
                # Is live downlink of the data being performed?
                if commsSystems[j].groundReqObs == 1:
                    flag = any(commsSystems[j].groundStationInsight[-1,:])
                else:
                    flag = 1
                commsFlag = commsFlag * flag

        flag2 = radioFlag * strFlag * radFlag * commsFlag
    else:
        flag2 = 1;

    # If any flags are 0 then observation cannot take place
    obsFlag = flag1 * flag2

    return obsFlag
    