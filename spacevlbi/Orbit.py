# Orbit.py
#
# Contains the functions required to model a space telescope's orbit in spacevlbi.
#
# @author: BenHudson - 28/07/2024

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
    Propagation currently includes following perturbations: Earth J2 and J3
    harmonics. Additional perturbing forces can be added by editing the Force()
    function.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list 
    :param time: Current time in simulation, defaults to None
    :type time: str
    :param duration: Time across which to propagate orbits, defaults to None
    :type duration: int
    :param rSun: Sun position vector in ECI frame in metres, defaults to None
    :type rSun: float
    :param rMoon: Moon position vector in ECI frame in metres, defaults to None
    :type rMoon: float
    :return: Array of spaceTelescope objects
    :rtype: list
    """

    # Iterate through spaceTelescopes
    if spaceTelescopes:
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

    :param t0: Current simulation time in seconds, defaults to None
    :type t0: float 
    :param u_: Current spacecraft state vector in ECI, defaults to None
    :type u_: numpy.ndarray
    :param k: Gravitational parameter of central body, defaults to None
    :type k: float
    :return: Acceleration vector imparted by perturbations in m/s^2
    :rtype: numpy.ndarray
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
    """Calculate range and elevation from ground stations to space telescopes,
    measured in topocentric frame.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundStations: Array of GroundStation objects, defaults to None
    :type groundStations: list
    :param time: Current simulation time, defaults to None
    :type time: str
    :return: Array of SpaceTelescope objects
    :rtype: list
    :return: Array of GroundStation objects
    :rtype: list
    """
           
    # Iterate through ground stations, transform coordinates to ECI and 
    # calculate range and elevation
    if groundStations:
        for k in range(len(groundStations)):
            # Initialise range and elevation vectors
            satRange = np.zeros((1,len(spaceTelescopes))) << u.m
            elevation = np.zeros((1,len(spaceTelescopes)))
            
            # Iterate through space telescopes
            if spaceTelescopes:
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
                    #print(eciPosition)
                    #print(spacecraftPos)
                    # Calculate satellite - ground station range magnitude
                    satGroundRange = norm(spacecraftPos - eciPosition)
                    # Calculate satellite - ground station range vector
                    satGroundVector = (spacecraftPos - eciPosition).reshape((3,1))
        
                    # Calculate elevation of spacecraft w.r.t ground station
                    elRad = np.dot(satGroundVector[:,0], eciPosition[:,0]) / \
                                    (norm(satGroundVector)*norm(eciPosition))
                    el = degrees((np.pi/2) - arccos(elRad.value))
                    
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
                    
                #print(elevation)
            
    return spaceTelescopes, groundStations

