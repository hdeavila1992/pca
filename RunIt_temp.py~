import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt

import F_BASE as FB
import pca_01 as pca
import nurbs as rbs
import plottingScripts as plts
import preprocessor2D as pre2D
import linearElastoStaticsSolver as linElastStat
import postprocessor2D as post2D
import surfaceRefinements as srfn
import debugScripts as dbg_scrpt
from numpy import random

####################################################
###################MAIN PROBLEM#####################
####################################################


#exec(open("first_optimization.py").read() 
doRefinement='N'


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


gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)


parametricNodes,nodesInElement = pre2D.parametricGrid(Uinp,Vinp)
loadElements,loadFaces = pre2D.loadPreprocessingv2(parametricNodes,nodesInElement,neumannConditions)

#Calculo de area, funcion objetivo

totalArea_temp=dbg_scrpt.calculateAreaAndLength(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,loadElements,loadFaces)

print('Area total!!!<------')
print(totalArea_temp)

np.savetxt('Area_data.dat',np.array([0,totalArea_temp]) )

#with open('Area_data.dat', 'a') as f:
#        f.write(str(totalArea_temp)+'\n')

dirichletCtrlPts,axisRestrictions = pre2D.dirichletBCPreprocessing(Pinp,displacementConditions)

#pre2D.plotGeometry(Uinp,Vinp,pinp,qinp,Pinp,winp,dirichletCtrlPts,displacementConditions,neumannConditions,parametricNodes,nodesInElement,loadElements,loadFaces)


dMat = linElastStat.elasticMatrix(E,nu)
K,F = linElastStat.assemblyWeakForm(Uinp,Vinp,winp,pinp,qinp,Pinp,parametricNodes,nodesInElement,gaussLegendreQuadrature,dMat,rho,loadElements,loadFaces,neumannConditions)

Kred,Fred,removedDofs,totalDofs = linElastStat.boundaryConditionsEnforcement(K,F,dirichletCtrlPts,axisRestrictions,u0,displacementConditions[0][2])

dtotal,D = linElastStat.solveMatrixEquations(Kred,Fred,totalDofs,removedDofs)
# print(D)

sx,sy,sxy,ux,uy=post2D.postProcessing(Uinp,Vinp,pinp,qinp,Pinp,D,winp,parametricNodes,nodesInElement,dtotal,dMat)


#print(np.amax(sx))
#print(np.amax(sy))
#print(np.amax(sxy))

svm_temp= np.sqrt( sx**2 - 2*sx*sy + sy**2 + 3*sxy**2 )


#with open('svm_data.txt', 'a') as f:
#        f.write(str(np.amax(svm_temp))+'\n')


np.savetxt('sy_data.dat',np.array([0,np.amax(sy)]) )
np.savetxt('sx_data.dat',np.array([0,np.amax(sx)]) )
np.savetxt('sxy_data.dat',np.array([0,np.amax(sxy)]) )
np.savetxt('Uy_data.dat',np.array([0,np.amax(uy)]) )
np.savetxt('Ux_data.dat',np.array([0,np.amax(ux)]) )
np.savetxt('svm_data.dat',np.array([0,np.amax(svm)]) )

print(np.amax(svm_temp))
print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('Flag here <--------------------------------')
print (Pinp)

#Rango


