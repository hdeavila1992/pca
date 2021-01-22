"""
Curve Fitting Interpolation
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg as la
import pca_01 as pca


#Function that creates the parameter vector 't' through the chord method
#Section 9.2.1. The NURBS Book
def tVector(qP):
    tv = np.zeros((len(qP[:,0]),1))

    tv[0] = 0
    tv[-1] = 1

    d = 0
    for i in range(1,len(qP[:,0])):
        d += np.sqrt(sum((qP[i,:] - qP[i-1,:])**2))

    for i in range(1,len(qP[:,0])):
        num = np.sqrt(sum((qP[i,:] - qP[i-1,:])**2))
        tv[i] = tv[i-1] + num/d

    return tv

def interpolateCurve(qP,p,method):

    sizeQ = np.shape(qP)
    n = sizeQ[0] - 1

    #t parameter vector generator
    tv = tVector(qP)

    if method == "uniform":
        U = pca.knotGeneratorUniform(n,p)
    elif method == "chord":
        U = pca.knotGeneratorChord(n,p,tv)
    else:
        print("Wrong knot generator parametrization")

    m = len(U) - 1
    M = np.zeros((n+1,n+1))

    for i in range (0,n+1):
        nb = pca.nFunction(U,p,tv[i])
        M[i,:] = nb

    ctrlpts = la.solve(M,qP)
    return ctrlpts,U

####################################################
###################MAIN PROBLEM#####################
####################################################

# qPoints = np.array([[0,0],[3,4],[-1,4],[-4,0],[-4,-3]])
# p = 3
# method = "chord"
#
# controlPoints,U = interpolateCurve(qPoints,p,method)
# cx,cy = pca.bSplineCurve(U,p,controlPoints.T)
# pca.plotInterpolatedCurve(cx,cy,controlPoints.T,qPoints.T)
