# Figures.py
#
# Functions to generate figures to show results of a spacevlbi simulation.
#
# @author: BenHudson - 28/07/2024

from matplotlib import pyplot as plt
import numpy as np
from numpy import radians
from astropy import constants as const
import os

###############################################################################
# OrbitPlot
###############################################################################

def OrbitPlot(spaceTelescopes):
    """Plot orbits of space telescopes.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    """
    
    if spaceTelescopes:
    
        # Generate figure axis
        fig = plt.figure(figsize=(8,8))
        print("Generating orbit plot...")
    
        ax = fig.add_subplot(221, projection = '3d')
        earthColour = 'g'
        
        # Generate Earth sphere
        r = const.R_earth
        u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
        x = r * np.cos(u) * np.sin(v)
        y = r * np.sin(u) * np.sin(v)
        z = r * np.cos(v)
        ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)
        
        # Initialise maximum array
        maxPos = 0
        
        # Iterate through space telescopes and plot orbits
        for i in range(len(spaceTelescopes)):
            position = spaceTelescopes[i].eciPosition
            if np.max(position.value) > maxPos:
                maxPos = np.max(position.value)
            name = spaceTelescopes[i].name
            ax.plot(position[:,0], position[:,1], position[:,2], label=name)
        
        # Configure axes
        axisLimit = 1.1 * maxPos
        ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit), \
               zlim=(-axisLimit, axisLimit))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.legend(loc="upper left")
        ax.set_aspect('equal', 'box')  
        
        ax = fig.add_subplot(222)
        # Iterate through space telescopes and plot orbits
        for i in range(len(spaceTelescopes)):
            position = spaceTelescopes[i].eciPosition
            ax.plot(position[:,0], position[:,1])
            
        # Add Earth circumference
        circle1 = plt.Circle((0, 0), const.R_earth.value, color= earthColour, \
                             fill=False)
        ax.add_artist(circle1)
        
        # Configure axes
        axisLimit = 1.1 * maxPos
        ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit))
        ax.set_xlabel('x [km]')
        ax.set_ylabel('y [km]')
        ax.set_aspect('equal', 'box')  
        
        ax = fig.add_subplot(223)
        # Iterate through space telescopes and plot orbits
        for i in range(len(spaceTelescopes)):
            position = spaceTelescopes[i].eciPosition
            ax.plot(position[:,0], position[:,2])
            
        # Add Earth circumference
        circle1 = plt.Circle((0, 0), const.R_earth.value, color= earthColour,\
                             fill=False)
        ax.add_artist(circle1)
        
        # Configure axes
        axisLimit = 1.1 * maxPos
        ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit))
        ax.set_xlabel('x [km]')
        ax.set_ylabel('z [km]')
        ax.set_aspect('equal', 'box')  
        
        ax = fig.add_subplot(224)
        # Iterate through space telescopes and plot orbits
        for i in range(len(spaceTelescopes)):
            position = spaceTelescopes[i].eciPosition
            name = spaceTelescopes[i].name
            ax.plot(position[:,1], position[:,2])
            
        # Add Earth circumference
        circle1 = plt.Circle((0, 0), const.R_earth.value, color= earthColour,\
                             fill=False)
        ax.add_artist(circle1)
        
        # Configure axes
        axisLimit = 1.1 * maxPos
        ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit))
        ax.set_xlabel('y [km]')
        ax.set_ylabel('z [km]')
        ax.set_aspect('equal', 'box')  
        
        plt.show()    
        
        #python program to check if a directory exists
        path = "Outputs"
        # Check whether the specified path exists or not
        isExist = os.path.exists(path)
        if not isExist:
           # Create a new directory because it does not exist
           os.makedirs(path)
           print("Creating Outputs folder")
           
        fig.savefig('Outputs/Orbit.pdf', format='pdf', \
                    bbox_inches='tight')
    else:
        print("Cannot generate Orbit plot, no space telescopes are modelled")

###############################################################################
# UvPlot
###############################################################################

def UvPlot(spaceTelescopes, groundTelescopes, allsky, plotLost=1):
    """Plot (u,v) coverage of a target source(s) generated by the simulated 
    array.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundTelescopes: Array of GroundTelescope objects, defaults to None
    :type groundTelescopes: list
    :param plotLost: Plot baselines lost due to functional constraints? 
        Defaults to 1
    :type plotLost: bool
    :param allysky: Plot all-sky coverage? Should be equal to allsky parameter
        defined in the Timeloop() function
    :type allysky: bool
    """
    
    if allsky == 0:
        # Marker size
        m = 10
        
        # Calculate number of stations
        stationNo = len(spaceTelescopes) + len(groundTelescopes)
        if stationNo > 1:
    
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot()
            print("Generating (u,v) coverage plot...")
            
            # Plot u = v = 0 lines
            ax.plot([-200, 200], [0,0], '--', c='grey')
            ax.plot([0, 0], [-200,200], '--', c='grey')
            
            # # Initialise maximum baseline parameter
            maxSpace = 0
            # Iterate through space telescopes and plot (u,v) coverage
            for i in range(len(spaceTelescopes)):
                for j in range(1, int(stationNo)+1):
                    uv = spaceTelescopes[i].baselines[0][:,j*3-3:j*3]
                    uvLost = spaceTelescopes[i].lostBaselines[0][:,j*3-3:j*3]
                    # Remove zeros indicating no baseline formed
                    uv = uv[~np.all(uv == 0, axis=1)]
                    uvLost = uvLost[~np.all(uvLost == 0, axis=1)]
                    if i == 0 and j == 1:
                        ax.scatter(uv[:,0], uv[:,1], s=m, c='blue', label="Space")
                        if plotLost == 1:
                            ax.scatter(uvLost[:,0], uvLost[:,1], s=m, c='red', \
                                      label="Lost")
                    else:
                        ax.scatter(uv[:,0], uv[:,1], s=m, c='blue')
                    if plotLost == 1:
                        ax.scatter(uvLost[:,0], uvLost[:,1], s=m, c='red') 
                uvMax = np.max(np.array([np.max(spaceTelescopes[i].baselines[0]), \
                                    np.max(spaceTelescopes[i].lostBaselines[0])]))
                if uvMax > maxSpace:
                    maxSpace = uvMax
                    
            # Initialise maximum baseline parameter
            maxGround = 0
            # Iterate through ground telescopes and plot (u,v) coverage
            for i in range(len(groundTelescopes)):
                for j in range(1, int(stationNo)+1):
                    # Vary colour depending on ground or space baseline
                    if j > int(stationNo - len(spaceTelescopes)):
                        colour = 'blue'
                    else:
                        colour = 'orange'                
                    uv = groundTelescopes[i].baselines[0][:,j*3-3:j*3]
                    uvLost = groundTelescopes[i].lostBaselines[0][:,j*3-3:j*3]
                    # Remove zeros indicating no baseline formed
                    uv = uv[~np.all(uv == 0, axis=1)]
                    uvLost = uvLost[~np.all(uvLost == 0, axis=1)]
                    if i == 0 and j == 1:
                        ax.scatter(uv[:,0], uv[:,1], s=m, c=colour, label="Ground")
                        if plotLost == 1:
                            ax.scatter(uvLost[:,0], uvLost[:,1], s=m, c='red') 
                    else:
                        ax.scatter(uv[:,0], uv[:,1], s=m, c=colour)
                        if plotLost == 1:
                            ax.scatter(uvLost[:,0], uvLost[:,1], s=m, c='red') 
                uvMax = np.max(groundTelescopes[i].baselines[0])
                if uvMax > maxGround:
                    maxGround = uvMax
            
            # Configure axes
            axisLimit = 1.3 * np.max(np.array([maxSpace, maxGround]))
            ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit))
            ax.set_xlabel('u [G$ \lambda $]')
            ax.set_ylabel('v [G$ \lambda $]')
            ax.set_aspect('equal', 'box')
            ax.legend(loc="upper right")
            ax.xaxis.set_inverted(True)
            plt.show()
            
            #python program to check if a directory exists
            path = "Outputs"
            # Check whether the specified path exists or not
            isExist = os.path.exists(path)
            if not isExist:
               # Create a new directory because it does not exist
               os.makedirs(path)
               print("Creating Outputs folder")
               
            fig.savefig('Outputs/UV.pdf', format='pdf', \
                        bbox_inches='tight')
        else:
            print("Cannot generate (u,v) plot, only one ground or space telescope has been modelled")
    
    else:
        
        # Marker size
        m = 0.1
        # Scale of individual (u,v) plots
        scale = 3
        
        # Calculate number of stations
        stationNo = len(spaceTelescopes) + len(groundTelescopes)
        
        if stationNo > 1:
    
            fig = plt.figure(figsize=(20,40))
            ax = fig.add_subplot()
            print("Generating all-sky (u,v) coverage plot...")
            
            # Define celestial sphere ranges
            raan = range(0,360,30)
            dec = range(-90,120,30)
            
            # Iterate through right ascension range
            for r in range(len(raan)):
                # Iterate through declination range
                for d in range(len(dec)):
            
                    # Iterate through space telescopes and plot (u,v) coverage
                    for i in range(len(spaceTelescopes)):
                        for j in range(1, int(stationNo)+1):
                            uv = spaceTelescopes[i].baselines[r*len(dec)+d]\
                                [:,j*3-3:j*3]
                            # Remove zeros indicating no baseline formed
                            uv = uv[~np.all(uv == 0, axis=1)]
                            if i == 0 and j == 1:
                                ax.scatter(uv[:,0]+raan[r]*scale, uv[:,1]+dec[d]*\
                                           scale, s=m, c='blue', label="Space")
                            else:
                                ax.scatter(uv[:,0]+raan[r]*scale, uv[:,1]+dec[d]*\
                                           scale, s=m, c='blue')
                            
                    # Iterate through ground telescopes and plot (u,v) coverage
                    for i in range(len(groundTelescopes)):
                        for j in range(1, int(stationNo)+1):
                            # Vary colour depending on ground or space baseline
                            if j >= int(stationNo - len(spaceTelescopes)):
                                colour = 'blue'
                            else:
                                colour = 'orange'                
                            uv = groundTelescopes[i].baselines[r*len(dec)+d]\
                                [:,j*3-3:j*3]
                            # Remove zeros indicating no baseline formed
                            uv = uv[~np.all(uv == 0, axis=1)]
                            if i == 0 and j == 1:
                                ax.scatter(uv[:,0]+raan[r]*scale, uv[:,1]+dec[d]*\
                                           scale, s=m, c=colour, label="Ground")
                            else:
                                ax.scatter(uv[:,0]+raan[r]*scale, uv[:,1]+dec[d]*\
                                           scale, s=m, c=colour)
            
            # Configure axes
            ax.set_xticks(range(0,scale*360,scale*30), ['0','30','60','90','120',\
                                                  '150','180','210','240','270',\
                                                  '300','330'])
            ax.set_yticks(range(scale*-90,scale*120,scale*30), ['-90','-60','-30',\
                                                  '0','30','60','90'])
            ax.set_xlabel('Right Ascension [$ \degree $]')
            ax.set_ylabel('Declination [$ \degree $]')
            ax.set_aspect('equal', 'box')
            plt.show()
            
            #python program to check if a directory exists
            path = "Outputs"
            # Check whether the specified path exists or not
            isExist = os.path.exists(path)
            if not isExist:
               # Create a new directory because it does not exist
               os.makedirs(path)
               print("Creating Outputs folder")
               
            fig.savefig('Outputs/UV_AllSky.pdf', format='pdf', \
                         bbox_inches='tight')    
        else:
            print("Cannot generate (u,v) plot, only one ground or space telescope has been modelled")
    
###############################################################################
# AttitudeSphere
###############################################################################

def AttitudeSphere(spaceTelescopes, telescopeSelect=0, azim=45, elev=30, \
                   plotAntenna=1, plotSTR=1, plotRad=1, plotComms=1, plotPanels=1, \
                   plotEarth=1, plotSun=1, plotMoon=1):
    """Plot attitude sphere of spacecraft in the body-fixed frame. Earth, Sun
    and Moon positions depicted. Note. For systems that have different Sun,
    Earth and Moon exclusion angles (e.g. radiator), only the largest is
    plotted.

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param telescopeSelect: Index of spaceTelescope array to plot attitude 
        sphere of, defaults to 0
    :type telescopeSelect: int
    :param azim: Plot viewing angle, azimuth in X-Y plane, defaults to 45 degrees
    :type azim: float
    :param elev: Plot viewing angle, elevation in Z plane, defaults to 30 degrees
    :type elev: float
    :param plotAntenna: Plot antenna? Defaults to 1
    :type plotAntenna: bool
    :param plotSTR: Plot star trackers? Defaults to 1
    :type plotSTR: bool
    :param plotRad: Plot radiators? Defaults to 1
    :type plotRad: bool
    :param plotComms: Plot comms systems? Defaults to 1
    :type plotComms: bool
    :param plotPanels: Plot solar panels? Defaults to 1
    :type plotPanels: bool
    :param plotEarth: Plot Earth vector and limb? Defaults to 1
    :type plotEarth: bool
    :param plotSun: Plot Sun vector and limb? Defaults to 1
    :type plotSun: bool
    :param plotMoon: Plot Moon vector and limb? Defaults to 1
    :type plotMoon: bool
    """
    
    if spaceTelescopes:
    
        # Generate figure axis
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(projection = '3d')
        print("Generating attitude sphere...")
        
        textScale = 1.1
        axisScale = 1.25
        
        # Extract SpaceTelescope from array and plot attitude sphere
        spaceTelescope = spaceTelescopes[telescopeSelect]
        
        # Plot spacecraft body axis
        ax.plot([0,axisScale], [0,0], [0,0], c='black')
        ax.plot([0,0], [0,axisScale], [0,0], c='black')
        ax.plot([0,0], [0,0], [0,axisScale], c='black')
        ax.text(1.3, 0, 0, "X", color='black')
        ax.text(0, 1.3, 0, "Y", color='black')
        ax.text(0, 0, 1.3, "Z", color='black')
        
        # Plot antenna surface
        if plotAntenna==1:
            vector = spaceTelescope.radioPayloads[0].antBoresight
            [X,Y,Z] = Cone3D([0,0,0], 0.1*vector, 0, 0.1*np.tan(radians(80)), 20)
            ax.plot_surface(X, Y, Z, color='grey', linewidth=0, antialiased=False, \
                            alpha=0.5)
        
        # Plot star trackers
        starTrackers = spaceTelescope.starTrackers
        if spaceTelescope.strModel == 1 and plotSTR == 1 and starTrackers != []:
            for k in range(len(starTrackers)):            
                vector = starTrackers[k].strBoresight
                name = starTrackers[k].name
                # Plot 3D cone showing exclusion angles (only largest exclusion
                # angle is plotted)
                fov = radians(np.max([starTrackers[k].strSunExcl, \
                                      starTrackers[k].strEarthExcl, \
                                          starTrackers[k].strMoonExcl]))            
                ax.plot([0,vector[0]], [0,vector[1]], [0,vector[2]], c='red')
                ax.text(textScale*vector[0], textScale*vector[1], textScale*\
                        vector[2], name, color='red')
                # Plot 3D cone showing FOV
                if fov > radians(45):
                    scale = (1 - fov/radians(90))*2
                else:
                    scale = 1
                [X,Y,Z] = Cone3D([0,0,0], scale*vector, 0, scale*np.tan(fov), 20)
                ax.plot_surface(X, Y, Z, color='red', linewidth=0, antialiased=False, \
                                alpha=0.1)
            
        # Plot radiator
        radiators = spaceTelescope.radiators
        if spaceTelescope.radModel == 1 and plotRad == 1 and radiators != []:
            for k in range(len(radiators)):            
                vector = radiators[k].radNorm
                name = radiators[k].name
                # Plot 3D cone showing exclusion angles (only largest radiator 
                # exclusion angle is plotted)
                fov = radians(np.max([radiators[k].radSunExcl, \
                                      radiators[k].radEarthExcl, \
                                          radiators[k].radMoonExcl]))
                ax.plot([0,vector[0]], [0,vector[1]], [0,vector[2]], c='red')
                ax.text(textScale*vector[0], textScale*vector[1], textScale*vector[2], \
                        name, color='red')
                if fov > radians(45):
                    scale = (1 - fov/radians(90))*2
                else:
                    scale = 1
                [X,Y,Z] = Cone3D([0,0,0], scale*vector, 0, scale*np.tan(fov), 20)
                ax.plot_surface(X, Y, Z, color='red', linewidth=0, antialiased=False, \
                                alpha=0.1)
        
        # Plot solar panel directions
        solarPanels = spaceTelescope.solarPanels
        if spaceTelescope.panelModel == 1 and plotPanels == 1 and \
            solarPanels != []:
            for k in range(len(solarPanels)):
                vector = solarPanels[k].panelNorm
                name = solarPanels[k].name
                ax.plot([0,vector[0]], [0,vector[1]], [0,vector[2]], c='red')
                ax.text(textScale*vector[0], textScale*vector[1], textScale*vector[2], \
                        name, color='red')
                
        # Plot comms system
        commsSystems = spaceTelescope.commsSystems
        if spaceTelescope.commsModel == 1 and plotComms == 1 and commsSystems != []:
            for k in range(len(commsSystems)):
                vector = commsSystems[k].commsNorm
                name = commsSystems[k].name
                fov = radians(commsSystems[k].commsFov)
                ax.plot([0,vector[0]], [0,vector[1]], [0,vector[2]], c='red')
                ax.text(textScale*vector[0], textScale*vector[1], textScale*vector[2], \
                        name, color='red')
                if fov > radians(45):
                    scale = (1 - fov/radians(90))*2
                else:
                    scale = 1
                [X,Y,Z] = Cone3D([0,0,0], scale*vector, 0, scale*np.tan(fov), 20)
                ax.plot_surface(X, Y, Z, color='red', linewidth=0, antialiased=False, \
                                alpha=0.1)
        
        # Plot initial Earth, Sun and Moon positions
        if plotEarth == 1:
            earthBody = spaceTelescope.earthBody[1,:]
            ax.plot([0,earthBody[0]], [0,earthBody[1]], [0,earthBody[2]], c='green', \
                    label="Earth")
            # Plot 3D cone showing Earth, Sun and Moon limbs
            fov = spaceTelescope.earthLimbAngle[1]
            [X,Y,Z] = Cone3D([0,0,0], earthBody, 0, np.tan(fov), 5)
            ax.plot_surface(X, Y, Z, color='green', linewidth=0, antialiased=False, \
                            alpha=0.1)
        if plotSun == 1:
            sunBody = spaceTelescope.sunBody[1,:]
            ax.plot([0,sunBody[0]], [0,sunBody[1]], [0,sunBody[2]], c='orange', \
                    label="Sun") 
            fov = spaceTelescope.sunLimbAngle[1]
            [X,Y,Z] = Cone3D([0,0,0], sunBody, 0, np.tan(fov), 5)
            ax.plot_surface(X, Y, Z, color='orange', linewidth=0, antialiased=False, \
                            alpha=0.1)
        if plotMoon == 1:
            moonBody = spaceTelescope.moonBody[1,:]
            ax.plot([0,moonBody[0]], [0,moonBody[1]], [0,moonBody[2]], c='blue', \
                    label="moon")
            fov = spaceTelescope.moonLimbAngle[1]
            [X,Y,Z] = Cone3D([0,0,0], moonBody, 0, np.tan(fov), 5)
            ax.plot_surface(X, Y, Z, color='blue', linewidth=0, antialiased=False, \
                            alpha=0.1)
            
        # Iterate through simulation time and plot Earth, Sun and Moon positions
        for i in range(2,len(spaceTelescope.earthBody)):
            # Plot initial Earth, Sun and Moon positions
            if plotEarth == 1:
                earthBody = spaceTelescope.earthBody[i,:]
                ax.plot([0,earthBody[0]], [0,earthBody[1]], [0,earthBody[2]], c='green')
                # Plot 3D cone showing Earth, Sun and Moon limbs
                fov = spaceTelescope.earthLimbAngle[1]
                [X,Y,Z] = Cone3D([0,0,0], earthBody, 0, np.tan(fov), 5)
                ax.plot_surface(X, Y, Z, color='green', linewidth=0, antialiased=False, \
                                alpha=0.1)
            if plotSun == 1:
                sunBody = spaceTelescope.sunBody[i,:]
                ax.plot([0,sunBody[0]], [0,sunBody[1]], [0,sunBody[2]], c='orange') 
                fov = spaceTelescope.sunLimbAngle[1]
                [X,Y,Z] = Cone3D([0,0,0], sunBody, 0, np.tan(fov), 5)
                ax.plot_surface(X, Y, Z, color='orange', linewidth=0, antialiased=False, \
                                alpha=0.1)
            if plotMoon == 1:
                moonBody = spaceTelescope.moonBody[i,:]
                ax.plot([0,moonBody[0]], [0,moonBody[1]], [0,moonBody[2]], c='blue')
                fov = spaceTelescope.moonLimbAngle[1]
                [X,Y,Z] = Cone3D([0,0,0], moonBody, 0, np.tan(fov), 5)
                ax.plot_surface(X, Y, Z, color='blue', linewidth=0, antialiased=False, \
                                alpha=0.1)
        
        # Configure axes
        axisLimit = 1.1
        ax.set(xlim=(-axisLimit, axisLimit), ylim=(-axisLimit, axisLimit), \
                zlim=(-axisLimit, axisLimit))
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.legend(loc="upper right")
        ax.set_aspect('equal', 'box')  
        ax.view_init(elev, azim)
        plt.show()
        
        #python program to check if a directory exists
        path = "Outputs"
        # Check whether the specified path exists or not
        isExist = os.path.exists(path)
        if not isExist:
           # Create a new directory because it does not exist
           os.makedirs(path)
           print("Creating Outputs folder")
           
        fig.savefig('Outputs/AttitudeSphere.pdf', format='pdf', \
                    bbox_inches='tight')
    else:
        print("Cannot generate Attitude Sphere plot, no space telescopes are modelled")

###############################################################################
# SolarPanelIncidence
###############################################################################

def SolarPanelIncidence(spaceTelescopes, simTime, telescopeSelect=0):
    """Plot angle of Sun incidence on the solar panels, measured with respect to
    the panel normal. 

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param simTime: Timeseries of simulation time, defaults to None
    :type simTime: list
    :param telescopeSelect: Index of spaceTelescope array to plot solar panel 
        incidence of, defaults to 0
    :type telescopeSelect: int
    """

    if spaceTelescopes:
        # Extract SpaceTelescope object
        spaceTelescope = spaceTelescopes[telescopeSelect]
        
        if spaceTelescope.panelModel == 1:
            
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot()
            print("Generating Solar Panel Angle plot...")
            
            # Extract time
            time = simTime.time.value[1:]
            for j in range(len(time)):
                time[j] = time[j][11:-4]
                
            # Iterate through solar panels and plot angle
            solarPanels = spaceTelescope.solarPanels
            for i in range(len(solarPanels)):
                name = solarPanels[i].name
                betaAngle = solarPanels[i].betaAngle
                ax.plot(time, betaAngle[1:], label=name)
                
            # Configure axes
            ax.set_xlabel('Time')
            ax.set_ylabel('Incidence Angle [$ \degree $]')
            ax.set_xticks(time[0::int(np.floor(len(time)/10))])
            ax.legend(loc="upper right")
            plt.xticks(rotation=90)
            plt.show()
            
            #python program to check if a directory exists
            path = "Outputs"
            # Check whether the specified path exists or not
            isExist = os.path.exists(path)
            if not isExist:
               # Create a new directory because it does not exist
               os.makedirs(path)
               print("Creating Outputs folder")
               
            fig.savefig('Outputs/SolarPanelIncidence.pdf', format='pdf', \
                        bbox_inches='tight')
        else:
            print("Solar panel modelling not turned on for selected space telescope")
    else:
        print("Cannot generate Solar Panel Incidence plot, no space telescopes are modelled")

###############################################################################
# GroundStationElevation
###############################################################################

def GroundStationElevation(spaceTelescopes, groundStations, simTime, \
                           telescopeSelect=0):
    """Plot elevation angle of selected space telescope at each ground station
    as a function of the simulation time. 

    :param spaceTelescopes: Array of SpaceTelescope objects, defaults to None
    :type spaceTelescopes: list
    :param groundStations: Array of GroundStation objects, defaults to None
    :type groundStations: list
    :param simTime: Timeseries of simulation time, defaults to None
    :type simTime: list
    :param telescopeSelect: Index of spaceTelescope array to plot ground 
        station elevation of, defaults to 0
    :type telescopeSelect: int
    """
    
    if spaceTelescopes:
        if len(groundStations) == 0:
            print("No ground stations are modelled")
        else:
            
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot()
            print("Generating Ground Station Access plot...")
            
            # Extract time
            time = simTime.time.value[1:]
            
            # Iterate through solar panels and plot angle
            # solarPanels = spaceTelescope.solarPanels
            for i in range(len(groundStations)):
                name = groundStations[i].name
                elevation = groundStations[i].satElev[1:,telescopeSelect]
                ax.plot(time, elevation, label=name)
                
            # Configure axes
            ylabel = spaceTelescopes[telescopeSelect].name + " Elevation Angle [$ \
                \degree $]"
            ax.set(ylim=(0, 90))
            ax.set_xlabel('Time')
            ax.set_ylabel(ylabel)
            ax.set_xticks(time[0::int(np.floor(len(time)/10))])
            ax.legend(loc="upper right")
            plt.xticks(rotation=90)
            plt.show()
            
            #python program to check if a directory exists
            path = "Outputs"
            # Check whether the specified path exists or not
            isExist = os.path.exists(path)
            if not isExist:
               # Create a new directory because it does not exist
               os.makedirs(path)
               print("Creating Outputs folder")
            
            fig.savefig('Outputs/GroundStationAccess.pdf', format='pdf',\
                        bbox_inches='tight')
    else:
        print("Cannot generate Ground Station Elevation plot, no space telescopes are modelled")

###############################################################################
# Miscellaneous
###############################################################################
    
def Cone3D(p0, p1, R0, R1, n):
    """Plot 3D cone. Used in AttitudeSphere figure.
 
    :param p0: Unit vector of vertex position, defaults to None
    :type p0: list
    :param p1: Unit vector of base position, defaults to None
    :type p1: list
    :param R0: Radius at vertex, defaults to None
    :type R0: float
    :param R1: Radius at base, defaults to None
    :type R1: float
    :param n: Number of radial steps used in cone surface generation, defaults 
        to None
    :type n: int
    """
    # vector in direction of axis
    v = p1 - p0
    # find magnitude of vector
    mag = np.linalg.norm(v)
    # unit vector in direction of axis
    v = v / mag
    # make some vector not in the same direction as v
    not_v = np.array([1, 1, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    # make vector perpendicular to v
    n1 = np.cross(v, not_v)
    # print n1,'\t',norm(n1)
    # normalize n1
    n1 /= np.linalg.norm(n1)
    # make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    # surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, n)
    theta = np.linspace(0, 2 * np.pi, n)
    # use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)
    R = np.linspace(R0, R1, n)
    # generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + R *
               np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in \
              [0, 1, 2]]
    return X, Y, Z