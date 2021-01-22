import numpy as np
import pca_01 as pca

####################################################
####################B-SPLINES#######################
####################################################

################B-SPLINES CURVES####################

def bSpline1DCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
    return cx

def bSplineCurve(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        cx[i] = nVec@px.T
        cy[i] = nVec@py.T
    return cx,cy

def bSplineCurveDerivative(U,p,P):
    numpoints = 101
    urank = np.linspace(U.min(),U.max(),numpoints)

    cprimex = np.zeros(len(urank))
    cprimey = np.zeros(len(urank))

    px = np.reshape(P[0],(1,len(P[0])))
    py = np.reshape(P[1],(1,len(P[1])))

    for i in range(len(urank)):
        dnVec = pca.derNFunction(U,p,urank[i])
        cprimex[i] = dnVec@px.T
        cprimey[i] = dnVec@py.T
    return cprimex,cprimey

################B-SPLINES SURFACES###################

def bivariateShapeFunction(U,V,p,q,u,v):
    nVecU = pca.nFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)

    n2 = nVecV.T@nVecU
    n2 = np.reshape(n2,(1,n2.shape[0]*n2.shape[1]))
    return n2

def bivariateDerivativeU(U,V,p,q,u,v):
    dnVecU = pca.derNFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)

    dn2u = nVecV.T@dnVecU
    dn2u = np.reshape(dn2u,(1,dn2u.shape[0]*dn2u.shape[1]))
    return dn2u

def bivariateDerivativeV(U,V,p,q,u,v):
    nVecU = pca.nFunction(U,p,u)
    dnVecV = pca.derNFunction(V,q,v)

    dn2v = dnVecV.T@nVecU
    dn2v = np.reshape(dn2v,(1,dn2v.shape[0]*dn2v.shape[1]))
    return dn2v

def bSpline2DField(U,V,p,q,px,py):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            n2 = bivariateShapeFunction(U,V,p,q,urank[i],vrank[j])

            cx[i,j] = n2@px
            cy[i,j] = n2@py
    return cx,cy

def bSplineSurface(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = pca.nFunction(U,p,urank[i])
            nVecj = pca.nFunction(V,q,vrank[j])

            cx[i,j] = nVecj@px@nVeci.T
            cy[i,j] = nVecj@py@nVeci.T
            cz[i,j] = nVecj@pz@nVeci.T
    return cx,cy,cz

def bSplineSurface2(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cx = np.zeros([len(urank),len(vrank)])
    cy = np.zeros([len(urank),len(vrank)])
    cz = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            n2 = bivariateShapeFunction(U,V,p,q,urank[i],vrank[j])

            cx[i,j] = n2@px
            cy[i,j] = n2@py
            cz[i,j] = n2@pz
    return cx,cy,cz

def bSplineSurfaceDerivativeU(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = pca.derNFunction(U,p,urank[i])
            nVecj = pca.nFunction(V,q,vrank[j])

            cprimex[i,j] = nVecj@px@dnVeci.T
            cprimey[i,j] = nVecj@py@dnVeci.T
            cprimez[i,j] = nVecj@pz@dnVeci.T
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeU2(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            dn2u = bivariateDerivativeU(U,V,p,q,urank[i],vrank[j])

            cprimex[i,j] = dn2u@px
            cprimey[i,j] = dn2u@py
            cprimez[i,j] = dn2u@pz
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            nVeci = pca.nFunction(U,p,urank[i])
            dnVecj = pca.derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@nVeci
            cprimey[i,j] = dnVecj.T@py@nVeci
            cprimez[i,j] = dnVecj.T@pz@nVeci
    return cprimex,cprimey,cprimez

def bSplineSurfaceDerivativeUV(U,V,p,q,px,py,pz):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    cprimex = np.zeros([len(urank),len(vrank)])
    cprimey = np.zeros([len(urank),len(vrank)])
    cprimez = np.zeros([len(urank),len(vrank)])
    for i in range(len(urank)):
        for j in range(len(vrank)):
            dnVeci = pca.derNFunction(U,p,urank[j])
            dnVecj = pca.derNFunction(V,q,vrank[j])

            cprimex[i,j] = dnVecj.T@px@dnVeci
            cprimey[i,j] = dnVecj.T@py@dnVeci
            cprimez[i,j] = dnVecj.T@pz@dnVeci
    return cprimex,cprimey,cprimez
