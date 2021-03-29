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

#Data
E = 2.1e5 #Pa
nu = 0.3
rho = 0.0 #kg/m3
u0 = 0.0
tv = 2 #Pa	



#Geometria

Pinit=np.array([[0,0],[3.75,0],[7.5,0],[11.25,0],[15,0],[18.75,0],[22.5,0],[26.25,0],[30,0],[0,3],[3.75,3],[7.5,3],[11.25,3],[15,3],[18.75,3],[22.5,3],[26.25,3],[30,3],[0,6],[3.75,6],[7.5,6],[11.25,6],[15,6],[18.75,6],[22.5,6],[26.25,6],[30,6]])


winit = np.ones(27)

print(Pinit)
VV=np.size(Pinit,0)

Actv=np.ones(VV)

#Condicion para activar 
Actv[0:9]=0
#Actv[18]=0
print(Actv)
#Restricciones 



#Isogeometric routines

Uinit=FB.knot(9-1,1)
Vinit=FB.knot(3-1,1)

print(Uinit)
print(Vinit)

pinit = 1
qinit = 1


numGaussPoints = 4

#Condiciones de borde


#displacementConditions = [[0.0,0,"C"]]
#neumannConditions = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]


displacementConditions = [[1.0,0,"S"]]
neumannConditions = [[[0.8133,0.0],[1.0,0.0],"normal",tv]]

#Esta linea es importante no borrar!!!!			
pos_cond=np.where(Actv==1)


#exec(open("debugScripts.py").read())
#exec(open("RunIt.py").read())

print('Final')
