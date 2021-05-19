# Python libraries
import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

#######################################################################
# DO NOT REMOVE THIS SEGMENT
import os
import sys

# Get current working directory
dir1 = os.getcwd()
# Insert .. command to go to the upper directory
dir2 = dir1 + '/..'
# Change directory
os.chdir(dir2)
# Get the new current working directory
dir3 = os.getcwd()
# Setting the package directory path for the modules execution
sys.path.append(dir3)
#######################################################################

# Local project
import src.plottingScripts as plts
import src.preprocessor2D as pre2D
import src.linearElastoStaticsSolver as linElastStat
import src.postprocessor2D as post2D
import src.surfaceRefinements as srfn
import src.debugScripts as dbg_scrpt

####################################################
################### MAIN PROBLEM ###################
####################################################

def mainProgram():
    #Data
    E = 1e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    materialProperties = [E,nu,rho]
    u0 = 0.0
    tv = 10 #Pa
    Rmax = 1.0
    Rmin = 0.7
    numGaussPoints = 4

    # gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[Rmin,0],[Rmin,Rmin],[0,Rmin],[Rmax,0],[Rmax,Rmax],[0,Rmax]])

    winit = np.array([[1],[0.5*np.sqrt(2)],[1],[1],[0.5*np.sqrt(2)],[1]])

    #Isogeometric routines
    Uinit = np.array([0,0,0,1,1,1])
    Vinit = np.array([0,0,1,1])

    pinit = 2
    qinit = 1

    doRefinement = 'Y'

    if doRefinement == 'Y':
        reflist = ['h','h','h','h']
        dirlist = ['U','V','U','V']
        Uinp,Vinp,pinp,qinp,Pinp,winp = srfn.surfaceRefinement(reflist,dirlist,Uinit,Vinit,pinit,qinit,Pinit,winit)
    else:
        Uinp = Uinit
        Vinp = Vinit
        pinp = pinit
        qinp = qinit
        Pinp = Pinit
        winp = winit

    displacementConditions = [[[0.0,Rmin],[0.0,Rmax],"S",0.0],[[Rmin,0.0],[Rmax,0.0],"S",0.0]]
    neumannConditions = [[[0.0,0.0],[1.0,0.0],"normal",tv]]

    parametricNodes,nodesInElement,surfacePreprocessing = pre2D.parametricGrid(Uinp,Vinp,pinp,qinp)
    loadElements,boundaryPreprocessing = pre2D.loadPreprocessing(parametricNodes,nodesInElement,neumannConditions,Uinp,Vinp,pinp,qinp)
    dirichletBCList = pre2D.dirichletBCPreprocessingOnFaces(Pinp,displacementConditions)
    numericalquadrature = pre2D.numericalIntegrationPreprocessing(numGaussPoints)
    
#    pre2D.plotGeometry(Uinp,Vinp,pinp,qinp,Pinp,winp,dirichletBCList,neumannConditions,boundaryPreprocessing)

    K,F = linElastStat.assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,surfacePreprocessing,numericalquadrature,materialProperties,boundaryPreprocessing,neumannConditions)

    Kred,Fred,removedDofs,totalDofs = linElastStat.boundaryConditionsEnforcement(K,F,dirichletBCList)

    dtotal,D = linElastStat.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)
    # print(D)

    post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,dtotal,surfacePreprocessing,materialProperties)

mainProgram()

#import cProfile
#import pstats
#profiler = cProfile.Profile()
#profiler.enable()
#mainProgram()
#profiler.disable()
## stats = pstats.Stats(profiler).sort_stats('ncalls')
#stats = pstats.Stats(profiler).sort_stats('tottime')
#stats.print_stats()