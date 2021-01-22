import numpy as np
import numpy.linalg
import nurbs as rbs

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

def jacobian(U,V,w,p,q,pta,ptb,px,py):
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

    return jacob

def weightedJacobian(jac,paramgrad,gwpts,ipta,iptb):
    return abs(np.linalg.det(jac))*abs(np.linalg.det(paramgrad))*gwpts[ipta]*gwpts[iptb]

def strainDisplacementMatrix(U,V,w,p,q,pta,ptb,jacob):
    dN2u = rbs.dRatdU(U,V,w,p,q,pta,ptb)
    dN2v = rbs.dRatdV(U,V,w,p,q,pta,ptb)

    invJac = np.linalg.inv(jacob)
    dN2 = np.vstack((dN2u,dN2v))
    dN2dxi = invJac.T@dN2

    numpts = dN2dxi.shape[1]
    bMat = np.zeros((3,2*numpts))
    #dNx
    bMat[0,0::2] = dN2dxi[0,:]
    bMat[2,0::2] = dN2dxi[1,:]
    #dNy
    bMat[1,1::2] = dN2dxi[1,:]
    bMat[2,1::2] = dN2dxi[0,:]
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

def localStiffnessMatrix(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt,dmat):
    lke = np.zeros((2*px.shape[0],2*px.shape[0]))
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            bMat = strainDisplacementMatrix(U,V,w,p,q,coor[0][0],coor[0][1],jac)
            lke += (bMat.T@dmat@bMat)*wJac

    return lke

def localBodyVector(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt,rho):
    lbe = np.zeros((2*px.shape[0],1))
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8

    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
            lbe += (nMat.T@bvec)*wJac

    return lbe

def appliedLoadVector(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt,loadvalue,loadtype,paramaxis,rotmat):
    lle = np.zeros((2*px.shape[0],1))

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)

        Jac = jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
        jac1 = np.linalg.norm(Jac[:,paramaxis])
        jac2 = 0.5*np.sum(bpt-apt)

        if jac1 > 1e-6:
            unitTangetVec = Jac[:,paramaxis]/jac1
        else:
            unitTangetVec = np.zeros((2,1))

        if loadtype == "tangent":
            tvec = (loadvalue/abs(loadvalue))*unitTangetVec
        elif loadtype == "normal":
            unitNormalVec = rotmat@unitTangetVec
            tvec = (loadvalue/abs(loadvalue))*unitNormalVec
        else:
            print("Wrong load configuration")

        tvec = np.reshape(tvec,(2,1))
        nMat = shapeFunctionMatrix(U,V,w,p,q,coor[0][0],coor[0][1])
        lle += (nMat.T@tvec)*jac1*jac2*gaussweights[qj]

    return lle

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(U,V,w,p,q,P,paramnodes,nodeselem,gaussquad,dmat,rho,loadelems,loadfaces,neumannconditions):
    K = np.zeros((2*P.shape[0],2*P.shape[0]))
    F = np.zeros((2*P.shape[0],1))
    Fb = np.zeros((2*P.shape[0],1))
    Fl = np.zeros((2*P.shape[0],1))
    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    paramGrad = np.zeros((2,2))
    numElems = nodeselem.shape[0]

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    loadtype = neumannconditions[0][2]
    loadvalue = neumannconditions[0][3]

    for ielem in range(0,numElems):
        """
        - paramGrad is conformed by the difference between u values and v values
        - The difference in u values is determined by paramnodes[C][0] - paramnodes[A][0]
          where C and A stand for the corners of the parametric element ABCD and
          0 stands for the u component
        - C can be obtained as the 3rd element in a given row of the matrix nodeselem:
        ielem -> [A B C D]
                  0 1 2 3
        - Likewise, A is the 1st element in a ielem row of nodeselem
        - It means that uC = paramnodes[nodeselem[ielem][2]][0]
                        uA = paramnodes[nodeselem[ielem][0]][0]
                        vC = paramnodes[nodeselem[ielem][2]][1]
                        vA = paramnodes[nodeselem[ielem][0]][1]
        """
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        paramGrad[0][0] = 0.5*(uC - uA)
        paramGrad[1][1] = 0.5*(vC - vA)

        aPoint = np.array([uA,vA])
        cPoint = np.array([uC,vC])

        print("---")
        print("Element #",ielem)
        K += localStiffnessMatrix(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,dmat)
        Fb += localBodyVector(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint,rho)

        if ielem in loadelems:
            print('Loaded element')

            iside = loadelems.index(ielem)
            paramside = loadfaces[iside]

            if paramside == 0 or paramside == 1:
                startindex = paramside
                endindex = paramside + 1
            elif paramside == 2:
                startindex = paramside + 1
                endindex = paramside
            else:
                startindex = 3
                endindex = 0

            if paramside == 1 or paramside == 3:
                paramaxis = 1
            else:
                paramaxis = 0

            uB = paramnodes[nodeselem[ielem][endindex]][0]
            uA = paramnodes[nodeselem[ielem][startindex]][0]
            vB = paramnodes[nodeselem[ielem][endindex]][1]
            vA = paramnodes[nodeselem[ielem][startindex]][1]

            aPoint = np.array([uA,vA])
            bPoint = np.array([uB,vB])

            Fl += appliedLoadVector(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,aPoint,bPoint,loadvalue,loadtype,paramaxis,rotMat)

        print("---")

    F = Fb + Fl
    return K,F

################ MATRIX EQUATION SOLUTION ####################

def boundaryConditionsEnforcement(K,F,udisp,axisrestrict,ucond,restriction):
    if restriction == "C":
        print("Clamped condition")
        remdofs = np.array(udisp) - 1
        remdofs = np.hstack((2*remdofs,2*remdofs))
        remdofs.sort()
    elif restriction == "S":
        print("Supported condition")
        remdofs = 2*(np.array(udisp) - 1)
    else:
        print("Wrong restriction")

    restricteddofs = np.array(axisrestrict)
    remdofs = remdofs + restricteddofs
    remdofs.sort()

    print("First reduction")
    Fred = np.delete(F,remdofs,0)
    Kred = np.delete(K,remdofs,0)

    print("Modification of Freduced")
    for rdof in remdofs:
        Kcol = Kred[:,rdof]
        Kcol = np.reshape(Kcol,(Kcol.shape[0],1))
        Fred -= Kcol*ucond

    print("Second reduction")
    Kred = np.delete(Kred,remdofs,1)

    totaldofs = np.arange(F.shape[0])

    return Kred,Fred,remdofs,totaldofs

def solveMatrixEquations(Kred,Fred,totaldofs,remdofs):
    # Checking full rank in matrix
    mRank = np.linalg.matrix_rank(Kred)
    mRows = Kred.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        fullRank = True
    else:
        fullRank = False

    # fullRank = True
    if fullRank:
        print("The matrix has full rank. It is invertible")
        dred = np.linalg.solve(Kred,Fred)
    else:
        print("The matrix has not full rank. It is not invertible")
        dred = np.linalg.lstsq(Kred,Fred,rcond=None)[0]

    reduceddofs = np.setdiff1d(totaldofs,remdofs)
    dtotal = np.zeros((totaldofs.shape[0],1))
    dtotal[reduceddofs,:] = dred

    dx = dtotal[0::2]
    dy = dtotal[1::2]
    D = np.hstack((dx,dy))

    return dtotal,D
