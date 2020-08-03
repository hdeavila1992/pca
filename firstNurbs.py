#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 18:16:56 2020
Modified on Sun Aug 02 17:59:20 2020

@author: hernando
@collaborator: paternina
"""
import numpy as np
import pca_01 as pca

Px = np.array([0,1,3,4,5])
Py = np.array([0,1,2,1,-1])
U = np.array([0,0,0,1,2,3,3,3])
w = np.array([1,4,1,1,1])
p = 2

cx,cy = pca.nurbsCurve(U,p,Px,Py,w)
pca.plotCurve2d(cx,cy,Px,Py)


