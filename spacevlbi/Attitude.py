# -*- coding: utf-8 -*-
"""
Attitude.py

Contains the function required to propagate a space telescope's attitude state
in spacevlbi.

@author: BenHudson - 05/07/2024
"""

import numpy as np
from numpy.linalg import norm, inv
from numpy import cross, matmul, vstack, dot, degrees, radians, cos, sin, arccos
from astropy import constants as const
from astropy import units as u

def AttitudePropagation(spaceTelescopes, rSun, rMoon, sourceRa, sourceDec):
    """Propagate space telescope's attitude state. This function calculates
    attitude matrix required to point the antenna at the target source. The
    constraint axis is pointed in a direction perpendicular to the
    antenna direction. This function could be updated to provide further
    control over the constraint axis by changing the definition of r2.

    :param spaceTelescopes: Array of spaceTelescope objects, defaults to None
    :type spaceTelescopes: SpaceTelescope
    :param rSun: Sun position vector in ECI frame in metres, defaults to None
    :type rSun: float
    :param rMoon: Moon position vector in ECI frame in metres, defaults to None
    :type rMoon: float
    :param sourceRa: Right ascension of target source in degrees, defaults to None
    :type sourceRa: float
    :param sourceDec: Declination of target source in degrees, defaults to None
    :type sourceRa: float
    :return: spaceTelescopes: Array of spaceTelescope objects
    :rtype spaceTelescopes: SpaceTelescope
    """
    
    # Iterate through space telescopes
    for j in range(len(spaceTelescopes)):
        # Calculate unit vector of source in spacecraft inertial frame
        sourceECI = np.array([cos(radians(sourceDec))*cos(radians(sourceRa)), \
                              cos(radians(sourceDec))*sin(radians(sourceRa)), \
                              sin(radians(sourceDec))])
        earthSourceECI = (sourceECI / norm(sourceECI)) * const.kpc
        # Extract spacecraft position vector
        rECI = spaceTelescopes[j].eciPosition[-1,:] << u.m
        rSource = (earthSourceECI - rECI)/np.linalg.norm(earthSourceECI - rECI)
        # Update space telescope object
        spaceTelescopes[j].sourceInertial = vstack((spaceTelescopes[j].sourceInertial, \
                            rSource))
            
        # Calculate Sun and Moon vectors with respect to spacecraft
        rSatSun = (rSun - rECI)/norm(rSun - rECI)
        rSatMoon = (rMoon - rECI)/norm(rMoon - rECI)
        spaceTelescopes[j].sunInertial = vstack((spaceTelescopes[j].sunInertial, \
                            rSatSun))
        spaceTelescopes[j].moonInertial = vstack((spaceTelescopes[j].moonInertial, \
                            rSatMoon))
            
        # Calculate attitude matrix to point antenna at target
        b1 = spaceTelescopes[j].pointingVector  # Body pointing axis
        b2 = spaceTelescopes[j].constraintVector  # Body constraint axis
        
        r1 = rSource;  # Inertial pointing axis
        r2 = cross(r1, b1) / norm(cross(r1, b1))  # Inertial constraint axis
        
        # Calculate attitude matrix
        attMat = TRIAD(r1,r2,b1,b2)
        
        # Propagate spacecraft attitude, rotate spacecraft component vectors
        # Rotate antenna
        radioPayloads = spaceTelescopes[j].radioPayloads
        if radioPayloads != []:
            for k in range(len(radioPayloads)):
                antennaInertial = matmul(inv(attMat), radioPayloads[k].antBoresight.reshape((3,1)))
                radioPayloads[k].antennaInertial = vstack((radioPayloads[k].antennaInertial, \
                                antennaInertial.reshape((1,3))))
            spaceTelescopes[j].radioPayloads = radioPayloads
        
        # Rotate star trackers
        starTrackers = spaceTelescopes[j].starTrackers
        if starTrackers != [] and spaceTelescopes[j].strModel == 1:
            for k in range(len(starTrackers)):
                strInertial = matmul(inv(attMat), starTrackers[k].strBoresight.reshape((3,1)))
                starTrackers[k].strInertial = vstack((starTrackers[k].strInertial, \
                                strInertial.reshape((1,3))))
            spaceTelescopes[j].starTrackers = starTrackers
                
        # Rotate solar panels and calculate beta angle
        solarPanels = spaceTelescopes[j].solarPanels
        if solarPanels != [] and spaceTelescopes[j].panelModel == 1:
            for k in range(len(solarPanels)):            
                panelInertial = matmul(inv(attMat), solarPanels[k].panelNorm.reshape((3,1)))
                solarPanels[k].solarPanelInertial = vstack((solarPanels[k].solarPanelInertial, \
                                panelInertial.reshape((1,3))))
                betaAngle = abs(arccos(dot(rSatSun, panelInertial) / (norm(rSatSun) \
                                * norm(panelInertial))))
                solarPanels[k].betaAngle = vstack((solarPanels[k].betaAngle,\
                                degrees(betaAngle)))
            spaceTelescopes[j].solarPanels = solarPanels
           
        # Rotate radiator surface
        radiators = spaceTelescopes[j].radiators
        if radiators != [] and spaceTelescopes[j].radModel == 1:
            for k in range(len(radiators)):      
                radInertial = matmul(inv(attMat), radiators[k].radNorm.reshape((3,1)))
                radiators[k].radInertial = vstack((radiators[k].radInertial, \
                                radInertial.reshape((1,3))))
            spaceTelescopes[j].radiators = radiators
            
        # Rotate comms terminal
        commsSystems = spaceTelescopes[j].commsSystems
        if commsSystems != [] and spaceTelescopes[j].commsModel == 1:
            for k in range(len(commsSystems)):
                commsInertial = matmul(inv(attMat), commsSystems[k].commsNorm.reshape((3,1)))
                commsSystems[k].commsInertial = vstack((commsSystems[k].commsInertial, \
                                commsInertial.reshape((1,3))))
            spaceTelescopes[j].commsSystems = commsSystems

        # Rotate inertial vectors into spacecraft body frame
        sunBody = matmul(attMat, rSatSun.reshape((3,1)))
        spaceTelescopes[j].sunBody = vstack((spaceTelescopes[j].sunBody, \
                            sunBody.reshape((1,3))))
        moonBody = matmul(attMat, rSatMoon.reshape((3,1))) 
        spaceTelescopes[j].moonBody = vstack((spaceTelescopes[j].moonBody, \
                            moonBody.reshape((1,3))))
        earthBody = matmul(attMat, (rECI / norm(rECI)).reshape((3,1)))
        spaceTelescopes[j].earthBody = vstack((spaceTelescopes[j].earthBody, \
                            earthBody.reshape((1,3))))
        
    return spaceTelescopes
    

def TRIAD(r1,r2,b1,b2):
    """TRIAD algorithm to calculate the attitude matrix from two pairs of body
    and inertial constraint vectors.

    :param r1: Inertial pointing vector, defaults to None
    :type r1: SpaceTelescope
    :param r2: Inertial constraint vector, defaults to None
    :type r2: float
    :param b1: Body pointing vector, defaults to None
    :type b1: float
    :param b2: Body constraint vector, defaults to None
    :type b2: float
    :return: attMat: Attitude Matrix
    :rtype attMat: np.array
    """
    
    rCross = cross(r1,r2) / norm(cross(r1,r2))
    bCross = cross(b1,b2) / norm(cross(b1,b2))

    a1 = matmul(b1.reshape((3,1)), r1.reshape((1,3)))
    a2 = cross(b1, bCross).reshape((3,1))
    a3 = cross(r1, rCross).reshape((1,3))
    a4 = matmul(bCross.reshape((3,1)), rCross.reshape((1,3)))
        
    attMat = a1 + matmul(a2, a3) + a4
        
    return attMat

