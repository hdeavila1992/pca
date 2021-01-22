#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 19:23:10 2020

@author: hernando
@collaborator: paternina
This code creates B-splines surfaces
"""

import numpy as np
import pca_01 as pca
from numpy import linalg as LA
import bSplines as bs
import plottingScripts as plts

Px = np.array([[0,3,6,9],[0,3,6,9],[0,3,6,9]])
Py = np.array([[0,0,0,0],[2,2,2,2],[4,4,4,4]])
Pz = np.array([[0,3,3,0],[2,5,5,2],[0,3,3,0]])

U = np.array([0,0,0,0.5,1,1,1])
V = np.array([0,0,0,1,1,1])
p = 2
q = 2

# U = np.array([0,0,0,0.4,0.6,1,1,1])
# V = np.array([0,0,0,0.2,0.5,0.8,1,1,1])
# p = 2
# q = 2

# U = [0,0,0,0,0.25,0.5,0.75,1,1,1,1]
# V = [0,0,0,0.2,0.4,0.6,0.6,0.8,1,1,1]
# p = 3
# q = 2

# u = bs.uGenerator(0.0,1.0,0.01)
# v = bs.uGenerator(0.0,1.0,0.01)

cx,cy,cz = bs.bSplineSurface(U,V,p,q,Px,Py,Pz)
# plts.plottingSurface(cx,cy,cz,Px,Py,Pz)

Px2 = np.reshape(Px,(Px.shape[0]*Px.shape[1],1))
Py2 = np.reshape(Py,(Py.shape[0]*Py.shape[1],1))
Pz2 = np.reshape(Pz,(Pz.shape[0]*Pz.shape[1],1))

cx2,cy2,cz2 = bs.bSplineSurface2(U,V,p,q,Px2,Py2,Pz2)
# plts.plottingSurface(cx2,cy2,cz2,Px,Py,Pz)
print(LA.norm(cx-cx2))
print(LA.norm(cx-cx2)/LA.norm(cx))

#Examples with partial derivatives
# cPrimeXu,cPrimeYu,cPrimeZu = bs.bSplineSurfaceDerivativeU(U,V,p,q,Px,Py,Pz)
# plts.plotTangentSurface(cx,cy,cz,cPrimeXu,cPrimeYu,cPrimeZu,Px,Py,Pz,'Partial derivative with respect to u')

# Px2 = np.reshape(Px,(Px.shape[0]*Px.shape[1],1))
# Py2 = np.reshape(Py,(Py.shape[0]*Py.shape[1],1))
# Pz2 = np.reshape(Pz,(Pz.shape[0]*Pz.shape[1],1))

# cPrimeXu2,cPrimeYu2,cPrimeZu2 = bs.bSplineSurfaceDerivativeU2(U,V,p,q,Px2,Py2,Pz2)
# plts.plotTangentSurface(cx,cy,cz,cPrimeXu2,cPrimeYu2,cPrimeZu2,Px,Py,Pz,'Partial derivative with respect to u')
# print(LA.norm(cPrimeXu-cPrimeXu2))
# print(LA.norm(cPrimeXu-cPrimeXu2)/LA.norm(cPrimeXu))

# cPrimeXv,cPrimeYv,cPrimeZv = pca.bSplineSurfaceDerivativeV(U,V,p,q,Px,Py,Pz)
# pca.plotTangentSurface(cx,cy,cz,cPrimeXv,cPrimeYv,cPrimeZv,Px,Py,Pz,'Partial derivative with respect to v')

#The operation gives almost zero matrices
# cPrimeXuv,cPrimeYuv,cPrimeZuv = pca.bSplineSurfaceDerivativeUV(U,V,p,q,Px,Py,Pz)
# pca.plotTangentSurface(cx,cy,cz,cPrimeXuv,cPrimeYuv,cPrimeZuv,Px,Py,Pz,'Partial derivative with respect to u and v')
