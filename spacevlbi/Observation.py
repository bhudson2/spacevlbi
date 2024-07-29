# Observation.py
#
# Functions required to simulate VLBI observations of target source and calculate
# baselines in spacevlbi.
#
# @author: BenHudson - 28/07/2024

from spacevlbi.Constraints import ObsMask
from astropy import constants as const
import numpy as np
from numpy.linalg import norm
from numpy import dot, cross, radians, cos, sin
from spacevlbi import Constraints

###############################################################################
# Baselines
###############################################################################

def Baselines(i, spaceTelescopes, groundTelescopes, sourceRa, sourceDec, \
              frequency, simLength, timeStep, allsky):
    """Calculate baselines formed by an array of a given source (or evaluated
    for a range of locations across the celestial sphere if allsky == 1).
    Baselines that couldn't be formed due to the impact of a functional
    constraint are also calculated. 

    :param i: Current time step, defaults to None
    :type i: int
    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundTelescopes: Array of GroundTelescope objects, defaults to None
    :type groundTelescopes: list
    :param sourceRa: Right ascension of target source in degrees, defaults to None
    :type sourceRa: float
    :param sourceDec: Declination of target source in degrees, defaults to None
    :type sourceDec: float
    :param frequency: Observation frequency in Hz, defaults to None
    :type frequency: float
    :param simLength: Simulation duration in seconds, defaults to None
    :type simLength: int
    :param timeStep: Simulation time step in seconds, defaults to None
    :type timeStep: int
    :param allsky: Calculate all-sky coverage? Defaults to None
    :type allsky: bool    
    :return: Array of SpaceTelescope objects
    :rtype: list
    :return: Array of GroundTelescope objects
    :rtype: list
    """
    
    # Speed of light
    c = const.c

    # Calculate wavelength
    wavelength = c / frequency
    
    # Form raan and dec arrays if allsky is 1
    if allsky == 1:
        # Define celestial sphere ranges
        raan = range(0,360,30)
        declination = range(-90,120,30)
    else:
        raan = [sourceRa]
        declination = [sourceDec]
        
    # Iterate through right ascension range
    for r in range(len(raan)):
        # Iterate through declination range
        for d in range(len(declination)):
            
            ra = radians(raan[r])
            dec = radians(declination[d])
            
            # Convert source ra and dec to vector in ECI frame
            sourceECI = np.array([cos(dec) * cos(ra), \
                                 cos(dec) * sin(ra), \
                                 sin(dec)])
            sourceVec = (sourceECI/norm(sourceECI))
            # Calculate proj U and V directions
            projU = cross(np.array([0,0,1]),sourceVec)
            projU = projU/norm(projU)
            projV = -cross(projU, sourceVec)

            # Calculate baselines
            # Iterate through spaceTelescopes
            for j in range(len(spaceTelescopes)):
                
                # If spaceTelescope's baseline property has length 0, calculate
                # the required length and initialise zero array
                if len(spaceTelescopes[j].baselines) == 0:
                    rows = len(raan) * len(declination)
                    stationNo = len(spaceTelescopes) + len(groundTelescopes)
                    matrix = []
                    matrixLost = []
                    for _ in range(rows):
                        matrix.append(np.zeros((int(round(simLength/timeStep)),\
                                                int(stationNo)*3)))
                        matrixLost.append(np.zeros((int(round(simLength/timeStep)),\
                                                    int(stationNo)*3)))
                    spaceTelescopes[j].baselines  = matrix
                    spaceTelescopes[j].lostBaselines = matrixLost
                
                baselineCount = 1
                # Extract position vector of telescope
                eciVector1 = spaceTelescopes[j].eciPosition[-1,:]
                
                # Check whether source is insight of space telescope
                sourceVis1 = Constraints.SourceVisibility(raan[r], declination[d], \
                                              eciVector1)
            
                # Iterate through spaceTelescope array. Calculate baselines formed
                # with primary spaceTelescope of this iteration
                for l in range(len(spaceTelescopes)):
                    # Extract position vector of telescope
                    eciVector2 = spaceTelescopes[l].eciPosition[-1,:]
            
                    # Determine whether observation is possible with these two antenna
                    # based on operational limitations
                    if allsky == 0:
                        obsFlag = ObsMask(spaceTelescopes[j], spaceTelescopes[l])
                    else:
                        obsFlag = 1

                    # Check whether source is insight of space telescope
                    sourceVis2 = Constraints.SourceVisibility(raan[r], declination[d], \
                                                  eciVector2)
                    
                    # Form (u,v) points
                    u = dot((eciVector1-eciVector2)/wavelength, projU)/1e9
                    v = dot((eciVector1-eciVector2)/wavelength, projV)/1e9
                    vector = np.array([u.value,v.value,0]).reshape((3,1))
                    
                    if j != l:
                        # If any conditions prohibiting observations are true, zero
                        # baseline, else calculate baseline
                        if obsFlag == 0:
                            spaceTelescopes[j].baselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                            # Add lost baselines to space telescope property
                            spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                        elif obsFlag == 1 and sourceVis1 == 1 and sourceVis2 == 1:
                            # Calculate uvw coordinates
                            spaceTelescopes[j].baselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                            spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        else:
                            spaceTelescopes[j].baselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                            spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    else:
                        spaceTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    
                    # Increase baseline index
                    baselineCount = baselineCount + 1

                # Iterate through ground array
                for l in range(len(groundTelescopes)):
                    # Extract position vector of telescope
                    eciVector2 = groundTelescopes[l].eciPosition[-1,:]
            
                    # Determine whether observation is possible with these two antenna
                    # based on operational limitations
                    if allsky == 0:
                        obsFlag = ObsMask(spaceTelescopes[j], groundTelescopes[l])
                    else:
                        obsFlag = 1
                    
                    # Check whether source is insight of ground telescope
                    el = Constraints.Elevation(eciVector2, raan[r], \
                                            declination[d])
                    if el < 90 and el > groundTelescopes[l].minEl:
                        sourceVis2 = 1
                    else:
                        sourceVis2 = 0
                    
                    # Form (u,v) points
                    u = dot((eciVector1-eciVector2)/wavelength, projU)/1e9
                    v = dot((eciVector1-eciVector2)/wavelength, projV)/1e9
                    vector = np.array([u.value,v.value,0]).reshape((3,1))
                    
                    # If any conditions prohibiting observations are true, zero
                    # baseline, else calculate baseline
                    if obsFlag == 0:
                        spaceTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        # Add lost baselines to space telescope property
                        spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                    elif obsFlag == 1 and sourceVis1 == 1 and sourceVis2 == 1:
                        # Calculate uvw coordinates
                        spaceTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                        spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0];
                    else:
                        spaceTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        spaceTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    
                    # Increase baseline index
                    baselineCount = baselineCount + 1

            # Iterate through groundTelescope array
            for j in range(len(groundTelescopes)):
                
                # If groundTelescope's baseline property has length 0, calculate
                # the required length and initialise zero array
                if len(groundTelescopes[j].baselines) == 0:
                    rows = len(raan) * len(declination)
                    stationNo = len(spaceTelescopes) + len(groundTelescopes)
                    matrix = []
                    matrixLost = []
                    for _ in range(rows):
                        matrix.append(np.zeros((int(round(simLength/timeStep)),int(stationNo)*3)))
                        matrixLost.append(np.zeros((int(round(simLength/timeStep)),int(stationNo)*3)))
                    groundTelescopes[j].baselines  = matrix
                    groundTelescopes[j].lostBaselines = matrixLost
                                                                 
                baselineCount = 1
                # Extract position vector of telescope
                eciVector1 = groundTelescopes[j].eciPosition[-1,:];
                
                # Check whether source is insight of ground telescope
                el = Constraints.Elevation(eciVector1, raan[r], \
                                        declination[d])
                if el < 90 and el > groundTelescopes[j].minEl:
                    sourceVis1 = 1
                else:
                    sourceVis1 = 0
            
                # Iterate through groundTelescopes array
                for l in range(len(groundTelescopes)):
                    # Extract position vector of telescope
                    eciVector2 = groundTelescopes[l].eciPosition[-1,:]
            
                    # Determine whether observation is possible with these two antenna
                    # based on operational limitations
                    if allsky == 0:
                        obsFlag = ObsMask(groundTelescopes[j], groundTelescopes[l])
                    else:
                        obsFlag = 1
                    
                    # Check whether source is insight of ground telescope
                    el = Constraints.Elevation(eciVector2, raan[r], \
                                            declination[d])
                    if el < 90 and el > groundTelescopes[l].minEl:
                        sourceVis2 = 1
                    else:
                        sourceVis2 = 0
                        
                    # Form (u,v) points
                    u = dot((eciVector1-eciVector2)/wavelength, projU)/1e9
                    v = dot((eciVector1-eciVector2)/wavelength, projV)/1e9
                    vector = np.array([u.value,v.value,0]).reshape((3,1))
                        
                    if j != l:
                        # If any conditions prohibiting observations are true, zero
                        # baseline, else calculate baseline
                        if obsFlag == 0:
                            groundTelescopes[j].baselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                            # Add lost baselines to space telescope property
                            groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                        elif obsFlag == 1 and sourceVis1 == 1 and sourceVis2 == 1:
                            # Calculate uvw coordinates
                            groundTelescopes[j].baselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                            groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                                [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        else:
                            groundTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                            groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    else:
                        groundTelescopes[j].baselines[r*len(declination)+d]\
                        [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                        [i,baselineCount*3-3:baselineCount*3] = [0,0,0]

                    # Increase baseline index
                    baselineCount = baselineCount + 1

                # Iterate through spaceStation array
                for l in range(len(spaceTelescopes)):
                    # Extract position vector of telescope
                    eciVector2 = spaceTelescopes[l].eciPosition[-1,:]
                    
                    # Determine whether observation is possible with these two antenna
                    # based on operational limitations
                    if allsky == 0:
                        obsFlag = ObsMask(groundTelescopes[j], spaceTelescopes[l])
                    else:
                        obsFlag = 1
                    
                    # Check whether source is insight of space telescope
                    sourceVis2 = Constraints.SourceVisibility(raan[r], declination[d], \
                                                  eciVector2)

                    # Form (u,v) points
                    u = dot((eciVector1-eciVector2)/wavelength, projU)/1e9
                    v = dot((eciVector1-eciVector2)/wavelength, projV)/1e9
                    vector = np.array([u.value,v.value,0]).reshape((3,1))
            
                    # If any conditions prohibiting observations are true, zero
                    # baseline, else calculate baseline
                    if obsFlag == 0:
                        groundTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        # Add lost baselines to space telescope property
                        groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                    elif obsFlag == 1 and sourceVis1 == 1 and sourceVis2 == 1:
                        # Calculate uvw coordinates
                        groundTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = vector.reshape((1,3))
                        groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    else:
                        groundTelescopes[j].baselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                        groundTelescopes[j].lostBaselines[r*len(declination)+d]\
                            [i,baselineCount*3-3:baselineCount*3] = [0,0,0]
                    
                    baselineCount = baselineCount + 1
           
    return spaceTelescopes, groundTelescopes