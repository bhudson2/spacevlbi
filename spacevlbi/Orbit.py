# -*- coding: utf-8 -*-
"""
Orbit.py

Contains the functions required to model a space telescope's orbit in spacevlbi.

@author: BenHudson - 05/07/2024
"""

import numpy as np
from numpy.linalg import norm
from numpy import degrees, arccos
from astropy.coordinates import SkyCoord,GCRS,ITRS
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation, J3_perturbation

###############################################################################
# OrbitPropagation
###############################################################################

def OrbitPropagation(spaceTelescopes, time, duration, rSun, rMoon):
    """Propagate space telescopes' orbits using poliastro Orbit functionality.
       Propagation currently includes following perturbations: Earth J2, J3
       harmonics. Additional perturbing forces can be added by editing the Force
       function.

       Args:
           spaceTelescopes (obj): Array of SpaceTelescope objects
           time (Time): Current time in simulation
           duration (int): Time across which to propagate orbits
           rSun (float): Sun position vector in ECI frame, metres
           rMoon (float): Moon position vector in ECI frame, metres

       Returns:
           spaceTelescopes (obj): Array of SpaceTelescope objects
    """   

    # Iterate through spaceTelescopes
    for j in range(len(spaceTelescopes)):
        
        # Extract initial orbit from space telescope
        initPosition = spaceTelescopes[j].orbit
        
        # Propagate orbit using Cowell's method and defined perturbations
        currentPosition = initPosition.propagate(duration << u.s, \
                            method=CowellPropagator(f=Force))
            
        # Convert propagated orbit to ECI
        posNow = np.array(currentPosition.r*1000) << u.m
        velNow = np.array(currentPosition.v*1000) << (u.m / u.s)
        # Convert ECI to ECEF
        eci = SkyCoord(x=posNow[0],y=posNow[1],z=posNow[2],\
                       unit='m',frame='gcrs',representation_type='cartesian', \
                           obstime = time)
        ecef = (eci.transform_to(ITRS)).cartesian
        ecefNow = np.array([ecef.x.value, ecef.y.value, ecef.z.value]) << u.m
        
        # Update spaceTelescope object with new position
        spaceTelescopes[j].eciPosition = np.vstack((spaceTelescopes[j].eciPosition, posNow))
        spaceTelescopes[j].eciVelocity = np.vstack((spaceTelescopes[j].eciVelocity, velNow))
        spaceTelescopes[j].ecefPosition = np.vstack((spaceTelescopes[j].ecefPosition, ecefNow))
    return spaceTelescopes


###############################################################################
# Perturbations
###############################################################################

def Force(t0, u_, k):
    """Force model for performing orbit propagation.

       Args:
           t0 (float): Current simulation time, sec
           u_ (float): Current position, velocity state vector, m
           k (float): Gravitational parameter of central body

       Returns:
           accel (float): Acceleration vector imparted by perturbations, m/s**2
    """   

    du_kep = func_twobody(t0, u_, k)
    # J2 acceleration
    ax, ay, az = J2_perturbation(t0, u_, k, J2=Earth.J2.value, \
                    R=Earth.R.to(u.km).value)
    # J3 acceleration
    ax1, ay1, az1 = J3_perturbation(t0, u_, k, J3=Earth.J3.value, \
                    R=Earth.R.to(u.km).value)
    # Form acceleration vectors
    du_ad = np.array([0, 0, 0, ax, ay, az])
    du_ad1 = np.array([0, 0, 0, ax1, ay1, az1])
    # Total acceleration
    accel = du_kep + du_ad + du_ad1
    return accel


###############################################################################
# SatGroundAccess
###############################################################################

def SatGroundAccess(spaceTelescopes, groundStations, time):
    """Calculate range and elevation between space telescopes and groundstations,
       measured in topocentric frame.

       Args:
           spaceTelescopes (obj): Array of SpaceTelescope objects
           groundStations (obj): Array of GroundStation objects
           time (Time): Current time in simulation

       Returns:
           spaceTelescopes (obj): Array of SpaceTelescope objects
           groundStations (obj): Array of GroundStation objects
    """   
           
    # Iterate through ground stations, transform coordinates to ECI and 
    # calculate range and elevation
    for k in range(len(groundStations)):
        # Initialise range and elevation vectors
        satRange = np.zeros((1,len(spaceTelescopes))) << u.m
        elevation = np.zeros((1,len(spaceTelescopes)))
        
        # Iterate through space telescopes
        for j in range(len(spaceTelescopes)):
            # Extract spacecraft position
            spacecraftPos = spaceTelescopes[j].eciPosition[-1,:].reshape((3,1))
            # Extract ground station position
            groundECEF = groundStations[k].ecefPosition.reshape((3,1))
            # Convert ground station position to ECI
            ecef =SkyCoord(x=groundECEF[0],y=groundECEF[1],z=groundECEF[2],\
                            unit='m',frame='itrs',representation_type=\
                            'cartesian', obstime = time)
            eci = (ecef.transform_to(GCRS)).cartesian    
            eciPosition = np.array([eci.x.value, eci.y.value, \
                            eci.z.value]).reshape((3,1)) << u.m
            # Update ground station position
            groundStations[k].eciPosition = np.vstack((groundStations[k].eciPosition, \
                            eciPosition.reshape((1,3))))
            # Calculate satellite - ground station range magnitude
            satGroundRange = norm(spacecraftPos - eciPosition)
            # Calculate satellite - ground station range vector
            satGroundVector = (spacecraftPos - eciPosition).reshape((3,1))

            # Calculate elevation of spacecraft w.r.t ground station
            elRad = np.dot(satGroundVector[:,0], eciPosition[:,0]) / \
                            (norm(satGroundVector)*norm(eciPosition))
            el = degrees((np.pi/2) - arccos(elRad).value)
            
            # Assign to temporary vector
            satRange[0,j] = satGroundRange
            elevation[0,j] = el
            
        
        # Pad satRange and satElev parameters of ground stations to match the 
        # number of spacecraft.
        if len(groundStations[k].satRange) == 1:
            groundStations[k].satRange = np.append(groundStations[k].satRange,\
                            np.zeros((1,len(spaceTelescopes)-1)))
        if len(groundStations[k].satElev) == 1:
            groundStations[k].satElev = np.append(groundStations[k].satElev, \
                            np.zeros((1,len(spaceTelescopes)-1)))
            
        # Update ground stations
        groundStations[k].satRange = np.vstack((groundStations[k].satRange,\
                                                satRange))
        groundStations[k].satElev = np.vstack((groundStations[k].satElev,\
                                                elevation))
            
    return spaceTelescopes, groundStations

