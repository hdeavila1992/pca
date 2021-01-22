import numpy as np
import pca_01 as pca

#####################################################
##################CONTROL POINTS#####################
#####################################################

#Convert from real control points to homogeneous ones
def weightedControlPoints(P,w):
    Pw = np.hstack((P,np.ones((P.shape[0],1))))
    Pw *= w
    return Pw

#Convert from homogeneous to real projection
def geometricControlPoints(Pw):
    w = Pw[:,-1]
    w = np.reshape(w,(len(w),1))
    P = Pw[:,:-1]/w
    return P,w

#Convert list of control points to a spatial grid
def listToGridControlPoints(Pl,U,V,p,q):
    #Number of control points in the U direction
    NU = len(U) - p - 1

    #Number of control points in the V direction
    NV = len(V) - q - 1

    Pg = np.zeros((Pl.shape[1],NU,NV))

    Pg[0,:,:] = np.reshape(Pl[:,0],(NU,NV),order='F')
    Pg[1,:,:] = np.reshape(Pl[:,1],(NU,NV),order='F')
    Pg[2,:,:] = np.reshape(Pl[:,2],(NU,NV),order='F')

    return Pg

def gridToListControlPoints(Pg):
    Pl = np.zeros((Pg.shape[1]*Pg.shape[2],Pg.shape[0]))

    Pl[:,0] = np.reshape(Pg[0,:,:],(Pg.shape[1]*Pg.shape[2]),order='F')
    Pl[:,1] = np.reshape(Pg[1,:,:],(Pg.shape[1]*Pg.shape[2]),order='F')
    Pl[:,2] = np.reshape(Pg[2,:,:],(Pg.shape[1]*Pg.shape[2]),order='F')

    return Pl

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
    nVecU = pca.nFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)
    dnVecU = pca.derNFunction(U,p,u)

    nunv = nVecV.T@nVecU
    nunv = np.reshape(nunv,(1,nunv.shape[0]*nunv.shape[1]))
    dnunv = nVecV.T@dnVecU
    dnunv = np.reshape(dnunv,(1,dnunv.shape[0]*dnunv.shape[1]))

    W = nunv@w
    dW = dnunv@w

    dratfunc = (w.T)*(dnunv*W - nunv*dW)/(W**2)
    return dratfunc

def dRatdV(U,V,w,p,q,u,v):
    nVecU = pca.nFunction(U,p,u)
    nVecV = pca.nFunction(V,q,v)
    dnVecV = pca.derNFunction(V,q,v)

    nunv = nVecV.T@nVecU
    nunv = np.reshape(nunv,(1,nunv.shape[0]*nunv.shape[1]))
    nudnv = dnVecV.T@nVecU
    nudnv = np.reshape(nudnv,(1,nudnv.shape[0]*nudnv.shape[1]))

    W = nunv@w
    dW = nudnv@w

    dratfunc = (w.T)*(nudnv*W - nunv*dW)/(W**2)
    return dratfunc

##################NURBS CURVES######################

def nurbsCurve(U,p,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)

    cx = np.zeros(len(urank))
    cy = np.zeros(len(urank))

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))

    for i in range(len(urank)):
        nVec = pca.nFunction(U,p,urank[i])
        ratFunc = (nVec*w.T)/(nVec@w)
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

def nurbs2DBoundary(U,V,p,q,P,w):

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))

    boundarycoor = []

    boundlimits = [[np.array([0.0,0.0]),np.array([1.0,0.0])],
                  [np.array([1.0,0.0]),np.array([1.0,1.0])],
                  [np.array([1.0,1.0]),np.array([0.0,1.0])],
                  [np.array([0.0,1.0]),np.array([0.0,0.0])]]

    numpt = 11
    for bndlim in boundlimits:
        parampath = np.linspace(bndlim[0],bndlim[1],numpt,endpoint=True)
        coor = np.zeros((parampath.shape[0],2))
        ipath = 0
        for ppath in parampath:
            ratFunc = ratFunction(U,V,w,p,q,ppath[0],ppath[1])

            coor[ipath,0] = ratFunc@px
            coor[ipath,1] = ratFunc@py
            ipath += 1

        boundarycoor.append(coor)

    for bc in range(len(boundarycoor)):
        if bc == 0:
            boundarycoor1 = boundarycoor[bc]
        else:
            boundarycoor1 = np.vstack((boundarycoor1,boundarycoor[bc]))

    cx = boundarycoor1[:,0]
    cy = boundarycoor1[:,1]
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

            cx[i,j] = ratFunc@px
            cy[i,j] = ratFunc@py
            cz[i,j] = ratFunc@pz
    return cx,cy,cz

def nurbsSurfaceTangent(U,V,p,q,P,w):
    numpoints = 21
    urank = np.linspace(U.min(),U.max(),numpoints)
    vrank = np.linspace(V.min(),V.max(),numpoints)

    px = np.reshape(P[:,0],(len(P[:,0]),1))
    py = np.reshape(P[:,1],(len(P[:,1]),1))
    pz = np.reshape(P[:,2],(len(P[:,2]),1))

    cpx = np.zeros([len(urank),len(vrank)])
    cpy = np.zeros([len(urank),len(vrank)])
    cpz = np.zeros([len(urank),len(vrank)])
    for j in range(len(vrank)):
        for i in range(len(urank)):
            dratFunc = dRatdU(U,V,w,p,q,urank[i],vrank[j])

            cpx[i,j] = dratFunc@px
            cpy[i,j] = dratFunc@py
            cpz[i,j] = dratFunc@pz
    return cpx,cpy,cpz
