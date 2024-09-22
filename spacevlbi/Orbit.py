"""
Library for simulating space-based VLBI missions (spacevlbi)

Copyright 2024 Ben Hudson

Licensed under the GNU General Public License, Version 3.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.gnu.org/licenses/gpl-3.0.en.html

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

# Orbit.py
#
# Contains the functions required to model a space telescope's orbit in spacevlbi.
#
# @author: BenHudson - 22/09/2024

import numpy as np
from numpy.linalg import norm
from numpy import degrees, arccos
from astropy.coordinates import SkyCoord,GCRS,ITRS
from astropy import units as u
from poliastro.bodies import Earth, Moon, Sun
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation, J3_perturbation, \
    atmospheric_drag_exponential
from poliastro.constants import rho0_earth, H0_earth

###############################################################################
# OrbitPropagation
###############################################################################

def OrbitPropagation(spaceTelescopes, time, rSun, rMoon, i, timeStep):
    """Propagate space telescopes' orbits using poliastro Orbit functionality.
    Propagation currently includes following perturbations: Earth J2, J3
    harmonics, luni-solar third-body, atmospheric drag and solar radiation pressure. 
    Propulsive burns in the form of a low-thrust, long-duration electric 
    propulsion system can also be modelled. Additional perturbing forces can be 
    added by editing the Force() function.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list 
    :param time: Current time in simulation, defaults to None
    :type time: str
    :param rSun: Sun position vector in ECI frame in metres, defaults to None
    :type rSun: list
    :param rMoon: Moon position vector in ECI frame in metres, defaults to None
    :type rMoon: list
    :param i: Current time step, defaults to None
    :type i: int
    :param timeStep: simulation time step in seconds, defaults to None
    :type timeStep: int
    :return: Array of spaceTelescope objects
    :rtype: list
    """

    # Iterate through spaceTelescopes
    if spaceTelescopes:
        for j in range(len(spaceTelescopes)):
            
            # Extract initial orbit from space telescope
            initPosition = spaceTelescopes[j].orbit
            
            # Extract spacecraft properties
            mass = spaceTelescopes[j].mass
            areaDrag = spaceTelescopes[j].areaDrag * 1e-6  # Convert to km^2
            cD = spaceTelescopes[j].cD
            areaSolar = spaceTelescopes[j].areaSolar * 1e-6  # Convert to km^2
            cR = spaceTelescopes[j].cR
            gravityJ2 = spaceTelescopes[j].gravityJ2
            gravityJ3 = spaceTelescopes[j].gravityJ3
            atmosDrag = spaceTelescopes[j].atmosDrag
            solarPress = spaceTelescopes[j].solarPress
            solarFlux = spaceTelescopes[j].solarFlux
            gravityLuniSolar = spaceTelescopes[j].gravityLuniSolar
            
            duration = i * timeStep
            # Propagate orbit using Cowell's method and defined perturbations
            currentPosition = initPosition.propagate(duration << u.s, \
                                method=CowellPropagator(f=Force_Wrapper(rSun,\
                                rMoon,mass,areaDrag,cD,areaSolar,cR, gravityJ2, \
                                gravityJ3, atmosDrag, solarPress, \
                                solarFlux, gravityLuniSolar)))
                
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

def Force_Wrapper(rSun, rMoon, mass, areaDrag, cD, areaSolar, cR, gravityJ2=1, \
                  gravityJ3=1, atmosDrag=1, solarPress=1, solarFlux=1367, \
                  gravityLuniSolar=1):
    def Force(t0, u_, k):
        """Force model for performing orbit propagation.
    
        :param t0: Current simulation time in seconds, defaults to None
        :type t0: float 
        :param u_: Current spacecraft state vector in ECI, defaults to None
        :type u_: numpy.ndarray
        :param k: Gravitational parameter of central body, defaults to None
        :type k: float
        :param rSun: Sun position vector in ECI frame in metres, defaults to None
        :type rSun: list
        :param rMoon: Moon position vector in ECI frame in metres, defaults to None
        :type rMoon: list
        :param mass: Spacecraft mass in kg, defaults to None
        :type mass: float
        :param areaDrag: Spacecraft drag area in km^2, defaults to None
        :type areaDrag: float
        :param cD: Spacecraft drag coefficient, defaults to None
        :type cD: float
        :param areaSolar: Spacecraft solar area in km^2, defaults to None
        :type areaSolar: float
        :param cR: Spacecraft reflectivity coefficient, defaults to None
        :type cR: float
        :param gravityJ2: Model Earth gravity harmonic J2? Defaults to 1
        :type gravityJ2: BOOL
        :param gravityJ3: Model Earth gravity harmonic J3? Defaults to 1
        :type gravityJ3: BOOL
        :param atmosDrag: Model atmospheric drag? Defaults to 1
        :type atmosDrag: BOOL
        :param solarPress: Model solar radiation pressure? Defaults to 1
        :type solarPress: BOOL
        :param solarFlux: Solar flux in W/m^2, defaults to 1367
        :type solarFlux: float
        :param gravityLuniSolar: Model luni-solar point model gravity? Defaults to 1
        :type gravityLuniSolar: BOOL
        :return: Acceleration vector imparted by perturbations in m/s^2
        :rtype: numpy.ndarray
        """
    
        # Two body acceleration
        du_kep = func_twobody(t0, u_, k)
        accel = du_kep
        
        # Sun - Moon parameters
        kMoon = Moon.k.to(u.km**3 / u.s**2).value
        kSun = Sun.k.to(u.km**3 / u.s**2).value
        
        rSat = u_[0:3]
        rEarthMoon = rMoon.value / 1000
        rSatMoon = rEarthMoon - rSat
        
        rEarthSun = rSun.value / 1000
        rSatSun = rEarthSun - rSat
        
        if gravityJ2 == 1:
            # J2 acceleration
            ax, ay, az = J2_perturbation(t0, u_, k, J2=Earth.J2.value, \
                            R=Earth.R.to(u.km).value)
            accel = accel + np.array([0, 0, 0, ax, ay, az])
            
        
        if gravityJ3 == 1:
            # J3 acceleration
            ax1, ay1, az1 = J3_perturbation(t0, u_, k, J3=Earth.J3.value, \
                            R=Earth.R.to(u.km).value)
            accel = accel + np.array([0, 0, 0, ax1, ay1, az1])            
            
        if gravityLuniSolar == 1:
            # Luni-Solar acceleration
            ax2, ay2, az2 = kMoon*((rSatMoon/(rSatMoon**3))-(rEarthMoon/(rEarthMoon**3)))
            ax3, ay3, az3 = kSun*((rSatSun/(rSatSun**3))-(rEarthSun/(rEarthSun**3)))
            
            accel = accel + np.array([0, 0, 0, ax2, ay2, az2]) + \
                np.array([0, 0, 0, ax3, ay3, az3])  
            
        if atmosDrag == 1:
            # Drag acceleration
            R = Earth.R.to(u.km).value
            rho0 = rho0_earth.to(u.kg / u.km**3).value  # kg/km^3
            H0 = H0_earth.to(u.km).value
            A_m = areaDrag / mass
            ax4, ay4, az4 = atmospheric_drag_exponential(t0, u_, k, R=R, C_D=cD, \
                                                          A_over_m=A_m, H0=H0, rho0=rho0)
            accel = accel + np.array([0, 0, 0, ax4, ay4, az4])
            
        if solarPress == 1:
            # Solar pressure acceleration
            flux = solarFlux * 1e-6  # W/km^2
            p = flux / 299792.458
            ax5, ay5, az5 = -p*cR*(areaSolar/mass)*(rSatSun/np.linalg.norm(rSatSun))
            accel = accel + np.array([0, 0, 0, ax5, ay5, az5])
        
        return accel
    
    return Force


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

