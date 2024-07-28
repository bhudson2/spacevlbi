# TimeLoop.py
#
# The primary simulation function of spacevlbi, executing all major functionality
# of the tool.
#
# @author: BenHudson - 28/07/2024

from spacevlbi.Orbit import OrbitPropagation, SatGroundAccess
from spacevlbi.Attitude import AttitudePropagation
from spacevlbi.Constraints import ObsLimits
from spacevlbi.Observation import Baselines
import numpy as np
from astropy.timeseries import TimeSeries
from astropy import units as u
from astropy.coordinates import SkyCoord,GCRS
from astropy.coordinates import solar_system_ephemeris, get_body
solar_system_ephemeris.set("jpl")

###############################################################################
# TimeLoop
###############################################################################

def TimeLoop(initTime, simLength, timeStep, spaceTelescopes, groundTelescopes,\
             groundStations, frequency, sourceRa, sourceDec, dutyCycle, \
             intTime, allsky):
    """Iterate through time loop and perform space VLBI simulation.

    :param initTime: Simulation start datetime, defaults to None
    :type initTime: str
    :param simLength: Simulation duration in seconds, defaults to None
    :type simLength: int
    :param timeStep: Simulation time step in seconds, defaults to None
    :type timeStep: int
    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundTelescopes: Array of GroundTelescope objects, defaults to None
    :type groundTelescopes: list
    :param groundStations: Array of GroundStation objects, defaults to None
    :type groundStations: list
    :param frequency: Observation frequency in Hz, defaults to None
    :type frequency: float
    :param sourceRa: Right ascension of target source in degrees, defaults to None
    :type sourceRa: float
    :param sourceDec: Declination of target source in degrees, defaults to None
    :type sourceDec: float
    :param dutyCycle: Time from start of one integration to the next in seconds, 
        defaults to None
    :type dutyCycle: int
    :param intTime: Integration time of instrument, defaults to None
    :type intTime: int
    :param allsky: Calculate all-sky coverage? Defaults to None
    :type allsky: bool
    :return: Array of SpaceTelescope objects
    :rtype: list
    :return: Array of GroundTelescope objects
    :rtype: list
    :return: Array of GroundStation objects
    :rtype: list
    :return: Timeseries of simulation time, defaults to None
    :rtype: list
    """
    
    # Create time series for simulation
    simTime = TimeSeries(time_start=initTime, time_delta=timeStep * u.s,\
                         n_samples=np.floor(simLength/timeStep))
    
    # Initialise parameters for modelling integration time and duty cycle
    lastObservation = 0;
    obsDuration = 0
    
    # Execute time loop
    for i in range(1,int(np.floor(simLength/timeStep))):
        # Display simulation progress
        print(str(i+1) + " / " + str(int((np.floor(simLength/timeStep)))))
        
        # Calculate ground telescope positions in ECI
        if groundTelescopes:
            for j in range(len(groundTelescopes)):
                posECEF = groundTelescopes[j].ecefPosition 
                # Transform ECEF TO ECI
                ecef = SkyCoord(x=posECEF[0,0],y=posECEF[0,1],z=posECEF[0,2],\
                                unit='m',frame='itrs',representation_type=\
                                'cartesian', obstime = simTime.time[i])
                eci = (ecef.transform_to(GCRS)).cartesian
                # Update GroundTelescope objects
                groundTelescopes[j].eciPosition = np.vstack((groundTelescopes[j]. \
                                eciPosition, np.array([eci.x.value, eci.y.value, \
                                eci.z.value]) << u.m))
                
        # Calculate Sun and Moon positions using SPICE kernels
        sun = get_body('sun', simTime.time[i], ephemeris='de432s')
        sun = (sun.transform_to(GCRS)).cartesian
        rSun = np.array([sun.x.value, sun.y.value, sun.z.value])*1000 \
                            << u.m
        moon = get_body('moon', simTime.time[i], ephemeris='de432s')
        moon = (moon.transform_to(GCRS)).cartesian
        rMoon = np.array([moon.x.value, moon.y.value, moon.z.value])*1000 \
                            << u.m
                
        # Spacecraft Orbit Propagation
        spaceTelescopes = OrbitPropagation(spaceTelescopes, simTime.time[i],\
                                         i*timeStep, rSun, rMoon)
        
        # Spacecraft Attitude Propagation
        spaceTelescopes = AttitudePropagation(spaceTelescopes, rSun, \
                            rMoon, sourceRa, sourceDec, i, timeStep)
        
        # Calculate Range / Elevation from Spacecraft to Ground Stations
        spaceTelescopes, groundStations = SatGroundAccess(spaceTelescopes, \
                            groundStations, simTime.time[i])
        
        # Calculate limits on observations for each station
        spaceTelescopes, groundTelescopes = ObsLimits(spaceTelescopes, \
                            groundTelescopes, groundStations, sourceRa,\
                            sourceDec, rSun, rMoon)
            
        # Calculate VLBI baselines
        if (i*timeStep - lastObservation >= dutyCycle) and (obsDuration <= \
                            intTime):
            spaceTelescopes, groundTelescopes = Baselines(i, spaceTelescopes, \
            groundTelescopes, sourceRa, sourceDec, frequency,simLength, \
            timeStep, allsky)
            obsDuration = obsDuration + timeStep
            if obsDuration >= intTime:
                lastObservation = i * timeStep - intTime
                obsDuration = 0

    return spaceTelescopes, groundTelescopes, groundStations, simTime
