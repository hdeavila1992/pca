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
E = 1e5 #Pa
nu = 0.31
rho = 0.0 #kg/m3
u0 = 0.0
tv = -10 #Pa	
Rmax = 1.0
Rmin = 0.5


#Geometria

#Pinit = np.array([[Rmin,0],[Rmin,Rmin],[0,Rmin],[Rmax,0],[Rmax,Rmax],[0,Rmax]])

#Pinit=np.array([[0,0],[0,1],[0.2,1],[0.4,1],[0.6,1],[0.8,1],[2,1],[2,0],[0.8,0],[0.6,0],[0.4,0],[0.2,0],[0,0]])

Pinit=np.array([[0,0],[0.5,0.0],[1.0,0],[0,0.5],[0.5,0.5],[1.0,0.5],[0,1.0],[0.5,1.0],[1,1.0]])

winit = np.array([[1],[1],[1],[1],[1],[1],[1],[1],[1]])

print(Pinit)
VV=np.size(Pinit,0)

Actv=np.ones(VV)

#Condicion para activar 
Actv[0]=0
Actv[3]=0
Actv[6:]=0
print(Actv)
#Restricciones 



#Isogeometric routines
Uinit = FB.knot(2,1) #np.array([0,0,0,1,1,1])
Vinit = FB.knot(2,1) #np.array([0,0,1,1])


pinit = 1
qinit = 1


numGaussPoints = 4

#Condiciones de borde

#displacementConditions = [[0.0,0,"S"],[0.0,1,"S"]]
#neumannConditions = [[[0.0,0.0],[0.5,0.0],"normal",tv]]

displacementConditions = [[0.0,0,"C"]]
neumannConditions = [[[1.0,0.0],[1.0,1.0],"tangent",tv]]
#exec(open("debugScripts.py").read())

#exec(open("RunIt.py").read())

#Esta linea es importante no borrar!!!!			
pos_cond=np.where(Actv==1)

print('Final')
