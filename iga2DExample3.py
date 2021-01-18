import numpy as np
from numpy.linalg import inv,det,solve
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

#Data
E = 1e5 #Pa
nu = 0.31
rho = 0.0 #kg/m3
u0 = 0.0
tv = -10 #Pa
# uDirichlet = [1,4,5,8,9,12]
# uAxis = [1,0,1,0,1,0]

numGaussPoints = 4
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

Pinit = np.array([[-1,0],[-1,np.sqrt(2)-1],[1-np.sqrt(2),1],[0,1],[-2.5,0],[-2.5,0.75],
              [-0.75,2.5],[0,2.5],[-4,0],[-4,4],[-4,4],[0,4]])

winit = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
              [1],[1],[1],[1],[1],[1]])

#Isogeometric routines
Uinit = np.array([0,0,0,0.5,1,1,1])
Vinit = np.array([0,0,0,1,1,1])

pinit = 2
qinit = 2

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

displacementConditions = [[0.0,0,"S"],[0.0,1,"S"]]
neumannConditions = [[[0.0,1.0],[0.5,1.0]]]

parametricNodes,nodesInElement = pre2D.parametricGrid(Uinp,Vinp)
loadElements,loadFaces = pre2D.loadPreprocessingv2(parametricNodes,nodesInElement,neumannConditions)
dirichletCtrlPts,axisRestrictions = pre2D.dirichletBCPreprocessing(Pinp,displacementConditions)

# cx,cy = rbs.nurbs2DField(Uinit,Vinit,pinp,qinp,Pinp,winp)
# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])),Pinp)

dMat = linElastStat.elasticMatrix(E,nu)
K,F = linElastStat.assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,dMat,rho,loadElements,loadFaces,tv,0)

Kred,Fred,removedDofs,totalDofs = linElastStat.boundaryConditionsEnforcement(K,F,dirichletCtrlPts,axisRestrictions,u0,displacementConditions[0][2])

dtotal,D = linElastStat.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)
# print(D)

post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,parametricNodes,nodesInElement,dtotal,dMat)
