# TimeLoop.py

# The primary simulation function of spacevlbi, executing all major functionality
#  of the tool.

# @author: BenHudson - 05/07/2024

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

       Args:
           initTime (time): Simulation start datetime
           simLength (int): Simulation duration, sec
           timeStep (int): Simulation time step, sec
           spaceTelescopes (obj): Array of SpaceTelescope objects
           groundTelescopes (obj): Array of GroundTelescope objects
           groundStations (obj): Array of GroundStation objects
           frequency (int): Observation frequency, Hz
           sourceRa (float): Right ascension of target source, deg
           sourceDec (float): Declination of target source, deg
           dutyCycle (int): Time from start of one integration to the next, sec
           intTime (int): Integration time of instrument, sec
           allysky (BOOL): Calculate all-sky coverage?

       Returns:
           spaceTelescopes (obj): Array of SpaceTelescope objects
           groundTelescopes (obj): Array of GroundTelescope objects
           groundStations (obj): Array of GroundStation objects
           simTime (Time): Timeseries of simulation time
    """   
    
    # Create time series for simulation
    simTime = TimeSeries(time_start=initTime, time_delta=timeStep * u.s,\
                         n_samples=round(simLength/timeStep))
    
    # Initialise parameters for modelling integration time and duty cycle
    lastObservation = 0;
    obsDuration = 0
    
    # Execute time loop
    for i in range(1,int(np.floor(simLength/timeStep))):
        # Display simulation progress
        print(str(i+1) + " / " + str(int((np.floor(simLength/timeStep)))))
        
        # Calculate ground telescope positions in ECI
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
                            rMoon, sourceRa, sourceDec)
        
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
