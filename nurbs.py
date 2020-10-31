import numpy as np
import pca_01 as pca

####################################################
######################NURBS#########################
####################################################

def ratFunction(U,V,w,p,q,u,v):
    nVecU = pca.nFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)

    n2 = nVecV.T@nVecU
    n2 = np.reshape(n2,(1,n2.shape[0]*n2.shape[1]))

    ratfunc = (n2*w.T)/(n2@w)
    return ratfunc

def dRatdU(U,V,w,p,q,u,v):
    dnVecU = pca.derNFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)

    dn2u = nVecV.T@dnVecU
    dn2u = np.reshape(dn2u,(1,dn2u.shape[0]*dn2u.shape[1]))

    ratfunc = (dn2u*w.T)/(dn2u@w)
    return ratfunc

def dRatdV(U,V,w,p,q,u,v):
    nVecU = pca.nFunction(U,p,u)
    dnVecV = pca.derNFunction(V,q,v)

    dn2v = dnVecV.T@nVecU
    dn2v = np.reshape(dn2v,(1,dn2v.shape[0]*dn2v.shape[1]))

    ratfunc = (dn2v*w.T)/(dn2v@w)
    return ratfunc

##################NURBS CURVES######################

def nurbsCurve(U,p,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(len(P[0]),1))
    py = np.reshape(P[1],(len(P[1]),1))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        ratFunc = (nVec*w)/(nVec@w)
        cx[i] = ratFunc@px
        cy[i] = ratFunc@py
    return cx,cy

#################NURBS SURFACES#####################

def nurbs2DField(U,V,p,q,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            ratFunc = ratFunction(U,V,w,p,q,urank[i],vrank[j])

            cx[i,j] = ratFunc@px
            cy[i,j] = ratFunc@py
    return cx,cy

def displacementField(U,V,p,q,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            ratFunc = ratFunction(U,V,w,p,q,urank[i],vrank[j])

            cx[i,j] = ratFunc@px
            cy[i,j] = ratFunc@py
    return cx,cy

def stressField(U,V,p,q,P,w,dMat):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            ratFunc = ratFunction(U,V,w,p,q,urank[i],vrank[j])

            cx[i,j] = ratFunc@px
            cy[i,j] = ratFunc@py
    return cx,cy

def nurbsSurface(U,V,p,q,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))
    pz = np.reshape(P[:,2],(len(P[:,2]),1))

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            ratFunc = ratFunction(U,V,w,p,q,urank[i],vrank[j])
            # ratFunc = (n2Func*w.T)/(n2Func@w)

            cx[i,j] = ratFunc@px
            cy[i,j] = ratFunc@py
            cz[i,j] = ratFunc@pz
    return cx,cy,cz
