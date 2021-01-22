import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt

import pca_01 as pca
import bSplines as bs
import curveFitting as cfit
import plottingScripts as plts
import refinements as rfn

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-5, atol=1e-8)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def scalingSegment(ua,ub,gausspt):
    localpts = np.zeros(len(gausspt))
    for i in range(len(localpts)):
        localpts[i] = 0.5*(ub - ua)*gausspt[i] + 0.5*(ub + ua)
    return localpts

def jacobian(U,p,u,PGeom):
    derNB = pca.derNFunction(U,p,u)
    jcb = derNB@PGeom[0].T
    return jcb

################ WEAK FORM INTEGRALS ####################

def localStiffnessMatrix(U,p,intgpt):
    derNB = pca.derNFunction(U,p,intgpt)
    lke = derNB.T*E@derNB
    return lke

def localForceVector(U,p,intgpt,xpt):
    nB = pca.nFunction(U,p,intgpt)
    if xpt < 0.5*L:#Be careful with this condition
        bx = 10.0
    else:
        bx = 0.0
    return nB.T*bx

################ EXACT SOLUTIONS ####################

def exactDisplacement(xpt):
    if (xpt < 0.5*L):
        # ux = xpt*((tx/E) + (0.5*bx*L/E)) - (0.5*bx/E)*xpt**2
        ux = 0.075*xpt - 0.005*xpt**2
    else:
        # ux = (tx/E)*(xpt - 0.5*L) + (0.5*tx*L/E) + ((0.125*bx*L**2)/E)
        ux = 0.025*(xpt - 5.0) + 0.25
    return ux

def exactStress(xpt):
    # bx = 10.0
    # L = 10
    # tx = 25
    if (xpt < 0.5*L):
        # sx = (tx + (0.5*bx)) - bx*xpt
        sx = 75.0 - 10.0*xpt
    else:
        sx = 25.0
    return sx

################ ISOGEOMETRIC ANALYSIS ####################

def preProcessing(numnodes,p,numgausspt):
    U = pca.knotGeneratorUniform(numnodes-1,p)
    Ured = U[p:-p] #this guarantees that u for evaluation goes from [0,1]
    # print(Ured) #Uncomment this line to check that knot vector is fine
    gaussQuadInfo = np.polynomial.legendre.leggauss(numgausspt)

    return U,Ured,gaussQuadInfo

def assemblyWeakForm(numnodes,U,p,Ured,PGeom,gaussquad,tx,u0):
    K = np.zeros((numnodes,numnodes))
    F = np.zeros((numnodes,1))

    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    totalLength = 0.0

    numElems = len(Ured) - 1
    for i in range(0,numElems):
        intgPoints = scalingSegment(Ured[i],Ured[i+1],gaussLegendrePoints)
        Ke = np.zeros((numnodes,numnodes))
        Fe = np.zeros((numnodes,1))
        for j in range(len(intgPoints)):
            jacob = jacobian(U,p,intgPoints[j],PGeom)
            print(jacob)
            totalLength += 0.5*(Ured[i+1] - Ured[i])*jacob*gaussLegendreWeights[j]
            Ke += 0.5*(Ured[i+1] - Ured[i])*(1.0/(jacob**2))*localStiffnessMatrix(U,p,intgPoints[j])*jacob*gaussLegendreWeights[j]

            nB = pca.nFunction(U,p,intgPoints[j])
            xpt = nB@PGeom[0].T
            Fe += 0.5*(Ured[i+1] - Ured[i])*localForceVector(U,p,intgPoints[j],xpt)*jacob*gaussLegendreWeights[j]

        K += Ke
        F += Fe

    #Applying Neumann BC
    F[-1] += tx

    # print(totalLength)
    # print(K)
    # print(F)

    #Applying Dirichlet BC
    dirNB = pca.nFunction(U,p,0.0)
    #Replacing basis funcion evaluated in u=0
    K[0,:] = dirNB

    choice = 1

    if choice == 1:
        #Removing the row and column with the Dirichlet BC
        Kred = K[1:,1:]
        #Modifying F column for suitable operation
        Kcol = np.reshape(K[:,0],(len(K[:,0]),1))
        #Updating F vector
        F -= Kcol*u0
        Fred = F[1:,:]

        #Solving the linear system
        d = la.solve(Kred,Fred)

        #Assembly the total solution
        dtotal = np.zeros((numnodes,1))
        dtotal[0] = u0
        dtotal[1:] = d
    else:
        F[0] = u0
        d = la.solve(K,F)
        dtotal = d

    return dtotal

def postProcessing(U,p,PGeom,dtotal):
    #Building the solution field
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(numpoints)
    ux = np.zeros(numpoints)
    sx = np.zeros(numpoints)
    uExact = np.zeros(numpoints)
    sExact = np.zeros(numpoints)

    for i in range(numpoints):
        nVec = pca.nFunction(U,p,urank[i])
        dnVec = pca.derNFunction(U,p,urank[i])
        cx[i] = nVec@PGeom[0].T
        ux[i] = nVec@dtotal
        uExact[i] = exactDisplacement(cx[i])
        sx[i] = E*dnVec@dtotal
        sExact[i] = exactStress(cx[i])

    return cx,ux,sx,uExact,sExact,PGeom

################ ERROR ANALYSIS ####################

def l2ErrorNorm(U,p,Ured,gaussquad,PGeom,dtotal):
    errorIntegral = 0.0

    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    numElems = len(Ured) - 1
    for i in range(0,numElems):
        intgPoints = scalingSegment(Ured[i],Ured[i+1],gaussLegendrePoints)
        lenSegment = Ured[i+1] - Ured[i]
        for j in range(len(intgPoints)):
            nB = pca.nFunction(U,p,intgPoints[j])
            xPt = nB@PGeom[0].T
            uE = exactDisplacement(xPt)
            ux = nB@dtotal
            errorIntegral += 0.5*lenSegment*((ux - uE)**2)*gaussLegendreWeights[j]
    return np.sqrt(errorIntegral)

def convergenceRate(Uinit,pdegree,numIter):
    # numIter = 5
    hVector = np.zeros((1,numIter))
    l2normVector = np.zeros((1,numIter))

    gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

    U = Uinit
    for i in range(numIter):
        print("Iteration #",i)
        Ured = U[pdegree:-pdegree]

        #Construction of geometry
        numpts = len(U) - pdegree - 1
        x = np.linspace(0,L,numpts)
        PGeom = np.array([x,np.zeros(numpts)])

        dtotal = assemblyWeakForm(numpts,U,pdegree,Ured,PGeom,gaussLegendreQuadrature,tx,u0)

        cx,ux,sx,uExact,sExact,PGeom = postProcessing(U,pdegree,PGeom,dtotal)

        #L2-Norm error
        hVector[0][i] = L/(numpts - 1)
        l2normVector[0][i] = l2ErrorNorm(U,pdegree,Ured,gaussLegendreQuadrature,PGeom,dtotal)
        print('L2-norm error',l2normVector[0][i])
        #Refining the knot vector
        X = 0.5*(Ured[0:-1] + Ured[1:])
        Ubar = np.hstack((U,X))
        Ubar.sort()
        U = Ubar

    xVec = np.log(hVector)
    yVec = np.log(l2normVector)

    xAvg = xVec.mean()
    slope = (yVec*(xVec - xAvg)).sum()/(xVec*(xVec - xAvg)).sum()

    print('Convergence Rate of the numerical solution: ',slope)
    # print(hVector)
    # print(l2normVector)
    # plt.plot(np.log(hVector),np.log(l2normVector))
    # plt.show()

####################################################
###################MAIN PROBLEM#####################
####################################################

#Data
L = 10.0
E = 1000
u0 = 0.0
tx = 25.0

numGaussPoints = 3

#Isogeometric routines
Uinit = np.array([0,0,0.5,1,1])
# Uinit = np.array([0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1])
# Uinit = np.array([ 0.,0.,0.0625,0.125,0.1875,0.25,0.3125,0.375,0.4375,0.5,0.5625,
# 0.625,0.6875,0.75,0.8125,0.875,0.9375,1.,1.])

# Uinit = np.array([0,0,0,0.1,0.3,0.5,1,1,1])
# Uinit = np.array([0,0,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1,1])
p = 1
Ured = Uinit[p:-p]
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)
# U,Ured,gaussLegendreQuadrature = preProcessing(numNodes,p,numGaussPoints)

#Construction of geometry
numPtsGeom = len(Uinit) - p - 1
x = np.linspace(0,L,numPtsGeom)
PGeom = np.array([x,np.zeros(numPtsGeom)])

dtotal = assemblyWeakForm(numPtsGeom,Uinit,p,Ured,PGeom,gaussLegendreQuadrature,tx,u0)
#
# cx,ux,sx,uExact,sExact,PGeom = postProcessing(Uinit,p,PGeom,dtotal)
#
# #L2-Norm error
# l2norm = l2ErrorNorm(Uinit,p,Ured,gaussLegendreQuadrature,PGeom,dtotal)
# print('L2-norm error',l2norm)

# plts.plotComparisonCurves(cx,ux,uExact,'u')
# plts.plotComparisonCurves(cx,sx,sExact,'s')

# plts.plot1DField(cx,ux,'u')
# plts.plot1DField(cx,sx,'s')

# convergenceRate(Uinit,p,5)
