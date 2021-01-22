#First B-spline curve
import numpy as np
import matplotlib.pyplot as plt
import pca_01 as pca
import plottingScripts as plts

# Examples with curves

# P = np.array([[-1,-1,1,1],[0,1,1,0]])
# U = np.array([0,0,0,1/2,1,1,1])
# p = 2

# P = np.array([[0,1,3,5,5,8,9],[0,2,4,2,0,0,3]])
# U = np.array([0,0,0,0,2,3,3,5,5,5,5])
# p = 3

# P = np.array([[1,2,3,4,5,6,7,8],[1,1,1,1,0,0,0,0]])
# p = 2
# U = pca.knotGeneratorUniform(len(P[0])-1,p)

# num = 5
# P = np.array([np.linspace(0,1,num),np.zeros(num)])
# p = 2
# U = np.array([0,0,0,1/3,2/3,1.0,1.0,1.0])

U = np.array([0,0,0.2,0.4,0.6,0.8,1,1])
V = np.array([0,0,0.5,1,1])
p = 1
q = 1

u = 0.0225403330759
v = 0.0563508326896

na = pca.nFunction(U,p,u)
dna = pca.derNFunction(U,p,u)
nb = pca.nFunction(V,q,v)
dnb = pca.derNFunction(V,q,v)

dnadu = nb.T@dna

# print(na)
print(dna)
print(nb)
# print(dnb)

print(dnadu)
print(np.reshape(dnadu,(1,dnadu.shape[0]*dnadu.shape[1])))

#B-Spline curve
# cx,cy = pca.bSplineCurve(U,p,P)
# plts.plotCurve2d(cx,cy,P)

# Examples with derivatives

#Derivative B-Spline

# P = np.array([[1,2,6,7,9],[1,5,5,1,3]])
# U = np.array([0,0,0,0.4,0.6,1,1,1])
# p = 2

# cX,cY = pca.bSplineCurve(U,p,P)
# plts.plotCurve2d(cX,cY,Px,Py,"yes","test1")
# plts.plotCurve2d(cX,cY,P,"no")

# cPrimeX,cPrimeY = pca.bSplineCurveDerivative(U,p,P)
# plts.plotTangentCurve2d(cX,cY,cPrimeX,cPrimeY,P,"yes","test2")
