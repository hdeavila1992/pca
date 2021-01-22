import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import preprocessor2D as pre2D
import linearElastoStaticsSolver as linElastStat
import postprocessor2D as post2D
import surfaceRefinements as srfn
import debugScripts as dbg_scrpt

####################################################
################## MAIN PROBLEM ####################
####################################################

def main():
    #Data
    L = 1.0
    H = 0.2
    E = 2e5 #Pa
    nu = 0.31
    rho = 0.0 #kg/m3
    u0 = 0.0
    tv = -1 #Pa
    # uDirichlet = [1,4,7]
    # uAxis = [1,0,1,0,1,0]

    numGaussPoints = 4
    gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    Pinit = np.array([[0,0],
                  	 [0.5*L,0],
                  	 [L,0],
                  	 [0,0.5*H],
                  	 [0.5*L,0.5*H],
                  	 [L,0.5*H],
                  	 [0,H],
                  	 [0.5*L,H],
                  	 [L,H]])

    winit = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1]])

    #Isogeometric routines
    # Uinit = np.array([0,0,0,1,1,1])
    # Vinit = np.array([0,0,0,1,1,1])
    Uinit = np.array([0,0,0.5,1,1])
    Vinit = np.array([0,0,0.5,1,1])

    pinit = 1
    qinit = 1

    doRefinement = 'N'

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

    displacementConditions = [[0.0,0,"C"]]
    neumannConditions = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]

    parametricNodes,nodesInElement = pre2D.parametricGrid(Uinp,Vinp)
    loadElements,loadFaces = pre2D.loadPreprocessingv2(parametricNodes,nodesInElement,neumannConditions)
    dirichletCtrlPts,axisRestrictions = pre2D.dirichletBCPreprocessing(Pinp,displacementConditions)

    # pre2D.plotGeometry(Uinp,Vinp,pinp,qinp,Pinp,winp,dirichletCtrlPts,displacementConditions,neumannConditions,parametricNodes,nodesInElement,loadElements,loadFaces)

    dMat = linElastStat.elasticMatrix(E,nu)
    K,F = linElastStat.assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,dMat,rho,loadElements,loadFaces,neumannConditions)

    Kred,Fred,removedDofs,totalDofs = linElastStat.boundaryConditionsEnforcement(K,F,dirichletCtrlPts,axisRestrictions,u0,displacementConditions[0][2])

    dtotal,D = linElastStat.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)
    # print(D)

    post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,parametricNodes,nodesInElement,dtotal,dMat)

main()

# import cProfile
# import pstats
# profiler = cProfile.Profile()
# profiler.enable()
# main()
# profiler.disable()
# # stats = pstats.Stats(profiler).sort_stats('ncalls')
# stats = pstats.Stats(profiler).sort_stats('tottime')
# stats.print_stats()
