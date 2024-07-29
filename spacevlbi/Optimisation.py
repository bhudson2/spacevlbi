# Optimisation.py
#
# Contains the function required to optimise the position of a spacecraft
# component with a specific Sun, Earth and/or Moon relationship. Optimisation
# provides vector along which a given component can be pointed in the 
# spacecraft's body-fixed frame that minimises the amount of times during a 
# simulation that a defined functional constraint is violated.
#
# @author: BenHudson - 28/07/2024

import numpy as np
from numpy.linalg import norm
from numpy import arctan, dot, degrees, arccos
from astropy import constants as const

def Optimisation(spaceTelescopes, sunExcl, earthExcl, \
                 moonExcl, telescopeSelect=0, direction="greaterthan"):
    """Iterate through the entire attitude sphere and determine the positions
    at which the defined Sun, Earth and Moon exclusion angles are violated the
    least during the simulation. This function can be used to find the optimal 
    location of specific spacecraft components (E.g. star trackers, radiators,
    etc.) to minimise the impact of their functional constraint on
    observations. The 'optimal' nature of the component must be able to be
    defined in terms of a Sun, Earth and/or Moon exclusion angle.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param telescopeSelect: Index of spaceTelescope array to optimise, defaults to 0
    :type telescopeSelect: int
    :param sunExcl: Sun exclusion angle in degrees. The direction parameter
        controls whether the Sun angle must be less than or greater than this
        value in order to perform an observation, defaults to None
    :type sunExcl: float
    :param earthExcl: Earth exclusion angle in degrees. The direction parameter
        controls whether the Earth angle must be less than or greater than this
        value in order to perform an observation, defaults to None
    :type earthExcl: float
    :param moonExcl:  Moon exclusion angle in degrees. The direction parameter
        controls whether the Moon angle must be less than or greater than this
        value in order to perform an observation, defaults to None
    :type moonExcl: float
    :param direction:  "lessthan" or "greaterthan". Whether the angle between 
        the unit normal vector and the celestial body should be greater than 
        or less than the earthExcl, sunExcl or moonExcl in order for an 
        observation to take place. E.g. for a star tracker, direction is set 
        to "greaterthan" and the sunExcl parameter is set to the exclusion angle of 
        the star tracker unit, defaults to "greaterthan"
    :type direction: str
    :return: List of unit vectors in the spacecraft body-fixed frame 
        and the associated percentage of the simulation for which the Sun, Earth and/or
        Moon exclusion angles were violated. I.e. smaller number is more optimal.
    :rtype: numpy.ndarray
    """
    
    if spaceTelescopes:
    
        # Extract space telescope and simulation length
        telescope = spaceTelescopes[telescopeSelect]
        simLength = len(spaceTelescopes[telescopeSelect].eciPosition)
        
        # Generate search space of unit vectors
        u = np.linspace(0, 2 * np.pi, 40)
        v = np.linspace(0, np.pi, 40)
        x = np.outer(np.cos(u), np.sin(v))
        y = np.outer(np.sin(u), np.sin(v))
        z = np.outer(np.ones(np.size(u)), np.cos(v))
        [m,n] = x.shape
        
        # Initialise results matrix
        results = np.zeros((m*n, 4))
        resultsCount = 0
        
        # Iterate through attitude sphere positions
        for i in range(m):
            for j in range(n):
                # Generate unit vector
                unitVec = np.array([x[i,j], y[i,j], z[i,j]])
                # Initialise status array
                status = np.zeros((simLength,1))
                for k in range(1,simLength):
                    # Extract Sun, Moon and Earth positions in body-fixed frame
                    # from simulation run
                    sunBody = telescope.sunBody[k,:].value
                    moonBody = telescope.moonBody[k,:].value
                    earthBody = telescope.earthBody[k,:].value
                    
                    # Calculate limb angles
                    rSun = telescope.rSun[k,:].value
                    rMoon = telescope.rMoon[k,:].value
                    moonRad = 1737.4e3
                    position = telescope.eciPosition[k,:].value
                    earthLimb = arctan(const.R_earth.value / norm(position))
                    sunLimb = arctan(const.R_sun.value / norm(rSun))
                    moonLimb = arctan(moonRad / norm(rMoon))
                    
                    # Calculate angles between unit vector and Earth, Sun and Moon
                    # limbs
                    sunAngle = degrees(abs(arccos(dot(unitVec, sunBody)/ \
                        (norm(unitVec)*norm(sunBody))))-sunLimb)
                    moonAngle = degrees(abs(arccos(dot(unitVec, moonBody)/ \
                        (norm(unitVec)*norm(moonBody))))-moonLimb)
                    earthAngle = degrees(abs(arccos(dot(unitVec, earthBody)/ \
                        (norm(unitVec)*norm(earthBody))))-earthLimb)
                        
                    #print(sunAngle)
                    #print(moonAngle)
                    #print(earthAngle)
                        
                    if direction == "lessthan":
                        if (sunAngle < sunExcl) and (earthAngle < earthExcl) and \
                            (moonAngle < moonExcl):
                            status[k] = 0
                        else:
                            status[k] = 1
                    elif direction == "greaterthan":
                        if (sunAngle > sunExcl) and (earthAngle > earthExcl) and \
                            (moonAngle > moonExcl):
                            status[k] = 0
                        else:
                            status[k] = 1
                    else:
                        print("Invalid direction parameter")
                        
                # Calculate for how many timesteps Sun, Earth and Moon
                # exclusion angles were not violated for
                results[resultsCount, 0:3] = unitVec
                results[resultsCount, 3] = ((len(status)-1-\
                                    np.count_nonzero(status))/(len(status)-1))*100
                resultsCount = resultsCount + 1
        
        # Sort results placing optimal configurations at the top
        sorted_indices = np.argsort(-results[:,-1])
        results = results[sorted_indices]
        
        # Return sorted list of results and print optimal position with mimimal
        # impact on observations
        print("Optimal position vector: " + str(results[0,0:3]))
        print("Fitness (% of simulation constraint is not violated): " + \
              str(round(results[0,3],2)) + "%")
        print("See function output 'results' for a full list of optimal positions")
    
    else:
        print("Cannot perform optimisation, no space telescopes are modelled")
        results = 0
        
    return results
        
    
    
    
    

