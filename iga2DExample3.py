import numpy as np
from numpy.linalg import inv,det,solve
# from scipy import linalg as la
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import pca_01 as pca
import nurbs as rbs
# import curveFitting as cfit
import plottingScripts as plts
# import refinements as rfn

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkRankMatrix(A):
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        print("The matrix has full rank. It is invertible")
    else:
        print("The matrix hast not full rank. It is not invertible")

def plotSparsity(A):
    fig = plt.figure()
    # plt.spy(A,markersize=5)
    plt.imshow(A,cmap=cm.viridis)
    plt.colorbar()
    plt.show()

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def geometricCoordinate(paramcoor,U,V,w,p,q,px,py):
    ratFunc = rbs.ratFunction(U,V,w,p,q,paramcoor[0][0],paramcoor[0][1])
    geomcoor = np.zeros((1,2))
    geomcoor[0][0] = ratFunc@px
    geomcoor[0][1] = ratFunc@py
    return geomcoor

def jacobian(U,V,w,p,q,pta,ptb,px,py,parentgrad):
    n2 = rbs.ratFunction(U,V,w,p,q,pta,ptb)
    dn2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dn2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

    dXdu = dn2u@px
    dXdv = dn2v@px
    dYdu = dn2u@py
    dYdv = dn2v@py

    jacob = np.zeros((2,2))

    jacob[0][0] = dXdu
    jacob[0][1] = dXdv
    jacob[1][0] = dYdu
    jacob[1][1] = dYdv

    jacob = jacob@parentgrad

    return jacob

def weightedJacobian(jac,gwpts,ipta,iptb):
    return abs(det(jac))*gwpts[ipta]*gwpts[iptb]

def strainDisplacementMatrix(U,V,w,p,q,pta,ptb,jacob):
    dN2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dN2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

    # print(jacob)
    invJac = inv(jacob)
    # print(invJac)
    dN2 = np.vstack((dN2u,dN2v))
    dN2dxi = invJac@dN2

    numpts = dN2dxi.shape[1]
    bMat = np.zeros((3,2*numpts))
    #dNx
    bMat[0,0::2] = dN2dxi[0]
    bMat[2,0::2] = dN2dxi[0]
    #dNy
    bMat[1,1::2] = dN2dxi[1]
    bMat[2,1::2] = dN2dxi[1]
    return bMat

def elasticMatrix(E,nu):
    dmat = np.zeros((3,3))
    dmat[0][0] = 1 - nu
    dmat[1][1] = 1 - nu
    dmat[2][2] = (1 - 2*nu)/2
    dmat[0][1] = nu
    dmat[1][0] = nu
    dmat *= E/((1+nu)*(1-2*nu))
    return dmat

def shapeFunctionMatrix(U,V,w,p,q,pta,ptb):
    N2 = rbs.ratFunction(U,V,w,p,q,pta,ptb)
    nMat = np.zeros((2,2*N2.shape[1]))

    nMat[0,0::2] = N2
    nMat[1,1::2] = N2
    return nMat

################ WEAK FORM INTEGRALS ####################

def localStiffnessMatrix(U,V,w,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad,dmat):
    lke = np.zeros((2*px.shape[0],2*px.shape[0]))
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            # print(coor)
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            bMat = strainDisplacementMatrix(U,V,w,p,q,coor[0][0],coor[0][1],jac)
            lke += (bMat.T@dmat@bMat)*wJac

    # print(lke)
    return lke

def localBodyVector(U,V,w,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad,rho):
    lbe = np.zeros((2*px.shape[0],1))
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
            lbe += (nMat.T@bvec)*wJac

    return lbe

# def appliedLoadVector(U,V,w,uval,vseg,p,q,px,py,gausspoints,gaussweights,load):
def appliedLoadVector(U,V,w,useg,vval,p,q,px,py,gausspoints,gaussweights,load):
    lle = np.zeros((2*px.shape[0],1))
    tvec = np.zeros((2,1))
    tvec[0][0] = -load

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        # coor = parametricCoordinate(uval,uval,vseg[0],vseg[1],gausspoints[qj],gausspoints[qj])
        coor = parametricCoordinate(useg[0],useg[1],vval,vval,gausspoints[qj],gausspoints[qj])
        gcoor = geometricCoordinate(coor,U,V,w,p,q,px,py)
        print("Geometric coor")
        print(gcoor)
        # jac2 = 0.5*(vseg[1] - vseg[0])
        jac2 = 0.5*(useg[1] - useg[0])
        du = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])
        print(du)
        dxdu = du@px
        dydu = du@py
        jac1 = np.sqrt(dxdu**2 + dydu**2)
        # print(dydu)
        nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
        lle += (nMat.T@tvec)*4.0*jac2*gaussweights[qj]

    return lle

def elementArea(U,V,w,useg,vseg,p,q,px,py,gausspoints,gaussweights,parentgrad):
    elemA = 0
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(useg[0],useg[1],vseg[0],vseg[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py,parentgrad)
            wJac = weightedJacobian(jac,gaussweights,qi,qj)
            elemA += 1.0*wJac

    return elemA

def elementLength(U,V,w,useg,vval,p,q,px,py,gausspoints,gaussweights):
    elemL = 0
    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        # coor = parametricCoordinate(uval,uval,vseg[0],vseg[1],gausspoints[qj],gausspoints[qj])
        coor = parametricCoordinate(useg[0],useg[1],vval,vval,gausspoints[qj],gausspoints[qj])
        gcoor = geometricCoordinate(coor,U,V,w,p,q,px,py)
        print("Geometric coor")
        print(gcoor)
        # jac2 = 0.5*(vseg[1] - vseg[0])
        jac2 = 0.5*(useg[1] - useg[0])
        dxdu = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])@px
        dydu = rbs.dRatdU(U,V,w,p,q,coor[0][0],coor[0][1])@py
        jac1 = np.sqrt(dxdu**2 + dydu**2)
        # print(dydu)
        elemL += 1.0*jac1*jac2*gaussweights[qj]

    return elemL

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(U,V,w,Ured,Vred,p,q,px,py,gaussquad,dmat,rho,uneu,vneu,load):
    K = np.zeros((2*px.shape[0],2*px.shape[0]))
    F = np.zeros((2*px.shape[0],1))
    Fb = np.zeros((2*px.shape[0],1))
    Fl = np.zeros((2*px.shape[0],1))
    totalArea = 0
    totalLength = 0
    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    parentElemGrad = np.zeros((2,2))
    numElems = 0

    for ji in range(len(Vred)-1):
        for ii in range(len(Ured)-1):
            #Inside each element
            numElems += 1

            parentElemGrad[0][0] = 0.5*(Ured[ii+1] - Ured[ii])
            parentElemGrad[1][1] = 0.5*(Vred[ji+1] - Vred[ji])

            print("Element #",numElems)
            print("U coor:",np.array([[Ured[ii],Ured[ii+1]]]))
            print("V coor:",np.array([[Vred[ji],Vred[ji+1]]]))
            print("---")
            usegment = np.array([Ured[ii],Ured[ii+1]])
            vsegment = np.array([Vred[ji],Vred[ji+1]])
            K += localStiffnessMatrix(U,V,w,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad,dmat)
            Fb += localBodyVector(U,V,w,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad,rho)
            totalArea += elementArea(U,V,w,usegment,vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,parentElemGrad)
            if Ured[ii] < uneu and abs(Vred[ji+1] - vneu) < 1e-5:
                # Fl += appliedLoadVector(U,V,w,Ured[ii+1],vsegment,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,load)
                Fl += appliedLoadVector(U,V,w,usegment,Vred[ji+1],p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,load)
                # totalLength += elementLength(U,V,w,usegment,Vred[ji+1],p,q,px,py,gaussLegendrePoints,gaussLegendreWeights)

            print("---")

    # print("Total Length")
    # print(totalLength)
    F = Fb + Fl
    return K,F,totalArea

def boundaryConditionsEnforcement(K,F,udisp,axisrestrict,ucond):
    remdofs = 2*(np.array(udisp) - 1)
    restricteddofs = np.array(axisrestrict)
    remdofs = remdofs + restricteddofs

    # remdofs = np.hstack((dofs1,dofs2))
    remdofs.sort()
    # print(remdofs)

    print("First reduction")
    Fred = np.delete(F,remdofs,0)
    Kred = np.delete(K,remdofs,0)

    print("Modification of Freduced")
    for u in udisp:
        Kcol = Kred[:,u]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*ucond

    print("Second reduction")
    Kred = np.delete(Kred,remdofs,1)

    return Kred,Fred,remdofs

def solveMatrixEquations(Kred,Fred,totaldofs,remdofs):
    dred = solve(Kred,Fred)
    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred
    return dtotal,dred

################ POSTPROCESSING ####################

def displacementField(U,V,w,p,q,D):
    ux,uy = rbs.nurbs2DField(U,V,p,q,D,w)
    return ux,uy

def stressField(U,V,w,p,q,ured,vred,P,dtot,dmat):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    parentgrad = np.zeros((2,2))

    sx = np.zeros([len(urank),len(vrank)])
    sy = np.zeros([len(urank),len(vrank)])
    sxy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):

            if urank[i] != 0.5 and vrank[j] != 1.0:
                xpcoor = urank[i]
                ypcoor = vrank[j]
            else:
                xpcoor = 0.52
                ypcoor = vrank[j]

            ii = pca.findKnotInterval(ured,xpcoor)
            ji = pca.findKnotInterval(vred,ypcoor)

            parentgrad[0][0] = 0.5*(ured[ii+1] - ured[ii])
            parentgrad[1][1] = 0.5*(vred[ji+1] - vred[ji])

            jac = jacobian(U,V,w,p,q,xpcoor,ypcoor,P[:,0],P[:,1],parentgrad)
            # print("---")
            # print(parentgrad)
            # print(urank[i])
            # print(vrank[j])
            # print(jac)
            # print("---")
            bmat = strainDisplacementMatrix(U,V,w,p,q,xpcoor,ypcoor,jac)

            svec = dmat@(bmat@dtot)
            sx[i,j] = svec[0]
            sy[i,j] = svec[1]
            sxy[i,j] = svec[2]

    return sx,sy,sxy

def postProcessing(U,V,w,p,q,ured,vred,P,D,dtot,dmat):
    cx,cy = rbs.nurbs2DField(U,V,p,q,P,w)
    ux,uy = displacementField(U,V,w,p,q,D)
    sx,sy,sxy = stressField(U,V,w,p,q,ured,vred,P,dtot,dmat)
    # plts.plotting2DField(cx,cy,ux,P,["Ux Displacement Field","[m]"])
    # plts.plotting2DField(cx,cy,sx,P,["Sx Component Stress Field","[Pa]"])
    print(sx.max())
    print(sx.min())

####################################################
###################MAIN PROBLEM#####################
####################################################

#Data
E = 1e5 #Pa
nu = 0.31
rho = 0 #kg/m3
u0 = 0.0
tv = 10 #Pa
uDirichlet = [1,4,5,8,9,12]
uAxis = [1,0,1,0,1,0]
uNeumann = 0.5
vNeumann = 1

numGaussPoints = 4
gaussLegendreQuadrature = np.polynomial.legendre.leggauss(numGaussPoints)

P = np.array([[-1,0],[-1,np.sqrt(2)-1],[1-np.sqrt(2),1],[0,1],[-2.5,0],[-2.5,0.75],
              [-0.75,2.5],[0,2.5],[-4,0],[-4,4],[-4,4],[0,4]])
w = np.array([[1],[0.5*(1+(1/np.sqrt(2)))],[0.5*(1+(1/np.sqrt(2)))],[1],[1],[1],
              [1],[1],[1],[1],[1],[1]])

#Isogeometric routines
Uinit = np.array([0,0,0,0.5,1,1,1])
Vinit = np.array([0,0,0,1,1,1])

p = 2
q = 2

# cx,cy = rbs.nurbs2DField(Uinit,Vinit,p,q,P,w)
# plts.plotting2DField(cx,cy,np.zeros((cx.shape[0],cx.shape[1])))

Ured = Uinit[p:-p]
Vred = Vinit[q:-q]

dMat = elasticMatrix(E,nu)
K,F,totalArea = assemblyWeakForm(Uinit,Vinit,w,Ured,Vred,p,q,P[:,0],P[:,1],gaussLegendreQuadrature,dMat,rho,uNeumann,vNeumann,tv)
# print(F)
# plotSparsity(K)

Kred,Fred,removedDofs = boundaryConditionsEnforcement(K,F,uDirichlet,uAxis,u0)

totaldofs = np.arange(2*P.shape[0])
dtotal,dred = solveMatrixEquations(Kred,Fred,totaldofs,removedDofs)

dx = dtotal[0::2]
dy = dtotal[1::2]
D = np.hstack((dx,dy))
# print(P)
# print(D)

postProcessing(Uinit,Vinit,w,p,q,Ured,Vred,P,D,dtotal,dMat)