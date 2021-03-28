import numpy as np
import nurbs as rbs
import plottingScripts as plts

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

################ POSTPROCESSING ####################

def nurbs2DField(numpoints,U,V,p,q,P,w,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

    cx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    cy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                ratFunc = rbs.ratFunction(U,V,w,p,q,urank[i],vrank[j])

                cx[iu*numpoints + i,jv*numpoints + j] = ratFunc@px
                cy[iu*numpoints + i,jv*numpoints + j] = ratFunc@py

    return cx,cy

def displacementField(numpoints,U,V,p,q,D,w,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    dx = np.reshape(D[:,0],(D.shape[0],1))
    dy = np.reshape(D[:,1],(D.shape[0],1))

    ux = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    uy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):
                ratFunc = rbs.ratFunction(U,V,w,p,q,urank[i],vrank[j])

                ux[iu*numpoints + i,jv*numpoints + j] = ratFunc@dx
                uy[iu*numpoints + i,jv*numpoints + j] = ratFunc@dy

    return ux,uy

# Improve this function
def stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem):
    numelemsu = len(np.unique(U)) - 1
    numelemsv = len(np.unique(V)) - 1
    numElems = nodeselem.shape[0]

    paramgrad = np.zeros((2,2))

    sx = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))
    sxy = np.zeros((numelemsu*numpoints,numelemsv*numpoints))

    maxSx = -1e20
    maxu = 0
    maxv = 0
    maxelem = 0

    for ielem in range(0,numElems):
        uC = paramnodes[nodeselem[ielem][2]][0]
        uA = paramnodes[nodeselem[ielem][0]][0]
        vC = paramnodes[nodeselem[ielem][2]][1]
        vA = paramnodes[nodeselem[ielem][0]][1]

        urank = np.linspace(uA,uC,numpoints)
        vrank = np.linspace(vA,vC,numpoints)

        jv = ielem//numelemsu
        iu = ielem%numelemsu

        for j in range(0,numpoints):
            for i in range(0,numpoints):

                xpcoor = urank[i]
                ypcoor = vrank[j]

                jac = jacobian(U,V,w,p,q,xpcoor,ypcoor,P[:,0],P[:,1])
                detJac = np.linalg.det(jac)

                if abs(detJac) > 1e-5:
                    bmat = strainDisplacementMatrix(U,V,w,p,q,xpcoor,ypcoor,jac)
                    svec = dmat@(bmat@dtot)
                else:
                    print("Singularity")
                    print(detJac)
                    # print(urank[i])
                    # print(vrank[j])
                    # xpcoor = 1.15*urank[i]
                    # xpcoor = urank[i]
                    # ypcoor = vrank[j]
                    xleft = urank[i] - 0.05
                    xright = urank[i] + 0.05
                    ydown = vrank[j] - 0.05
                    yup = vrank[j] + 0.05

                    if xleft < 1e-5:
                        xleft = 0.0

                    if xright > (1.0 + 1e-5):
                        xright = 1.0

                    if ydown < 1e-5:
                        ydown = 0.0

                    if yup > (1.0 + 1e-5):
                        yup = 1.0

                    apt = np.array([xleft,yup])
                    bpt = np.array([xright,yup])
                    cpt = np.array([xleft,ydown])
                    dpt = np.array([xright,ydown])

                    jaca = jacobian(U,V,w,p,q,apt[0],apt[1],P[:,0],P[:,1])
                    bmata = strainDisplacementMatrix(U,V,w,p,q,apt[0],apt[1],jaca)
                    sveca = dmat@(bmata@dtot)

                    jacb = jacobian(U,V,w,p,q,bpt[0],bpt[1],P[:,0],P[:,1])
                    bmatb = strainDisplacementMatrix(U,V,w,p,q,bpt[0],bpt[1],jacb)
                    svecb = dmat@(bmatb@dtot)

                    jacc = jacobian(U,V,w,p,q,cpt[0],cpt[1],P[:,0],P[:,1])
                    bmatc = strainDisplacementMatrix(U,V,w,p,q,cpt[0],cpt[1],jacc)
                    svecc = dmat@(bmatc@dtot)

                    jacd = jacobian(U,V,w,p,q,dpt[0],dpt[1],P[:,0],P[:,1])
                    bmatd = strainDisplacementMatrix(U,V,w,p,q,dpt[0],dpt[1],jacd)
                    svecd = dmat@(bmatd@dtot)

                    svec = 0.25*(sveca + svecb + svecc + svecd)

                sx[iu*numpoints + i,jv*numpoints + j] = svec[0]
                sy[iu*numpoints + i,jv*numpoints + j] = svec[1]
                sxy[iu*numpoints + i,jv*numpoints + j] = svec[2]

                if svec[0] > maxSx:
                    maxu = xpcoor
                    maxv = ypcoor
                    maxelem = ielem
                    maxSx = svec[0]
                    # print('----')
                    # print(maxu)
                    # print(maxv)
                    # print(maxSx)
                    # print('----')

    # print('Element with maximum value',maxelem)
    # print('u: ',maxu)
    # print('v: ',maxv)

    # Maximum stress value
     #maxSx = np.amax(sx)
     #maxSy = np.amax(sy)
     #maxSxy = np.amax(sxy)
     #print('Maximum stress value',maxSx)
    # # Index of the minimum stress value
     #imaxSx = np.where(sx == np.amax(sx))
     #listOfCordinates = list(zip(imaxSx[0], imaxSx[1]))
     #for cord in listOfCordinates:
      #   print(cord)

    return sx,sy,sxy

def stressTensor(uval,vval,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem):
    numElems = nodeselem.shape[0]

    paramgrad = np.zeros((2,2))

    # Finding the element
    iu = pca.findKnotInterval(np.unique(U),uval)
    jv = pca.findKnotInterval(np.unique(V),vval)

    ielem = jv*len(np.unique(U)) + iu
    print('# of element')
    print(ielem)

    if abs(uval - 0.5) > 1e-5:
        xpcoor = uval
        ypcoor = vval
    else:
        xpcoor = 1.05*uval
        ypcoor = vval

    print('u: ',xpcoor)
    print('v: ',ypcoor)

    jac = jacobian(U,V,w,p,q,xpcoor,ypcoor,P[:,0],P[:,1])
    bmat = strainDisplacementMatrix(U,V,w,p,q,xpcoor,ypcoor,jac)
    svec = dmat@(bmat@dtot)

    return svec[0]

def postProcessing(U,V,p,q,P,D,w,paramnodes,nodeselem,dtot,dmat):
    numpoints = 11
    cx,cy = nurbs2DField(numpoints,U,V,p,q,P,w,paramnodes,nodeselem)
    ux,uy = displacementField(numpoints,U,V,p,q,D,w,paramnodes,nodeselem)
    sx,sy,sxy = stressField(numpoints,U,V,p,q,P,w,dtot,dmat,paramnodes,nodeselem)
    # plts.plotting2DField(cx,cy,ux,P,["Ux Displacement Field","[m]"])
    #plts.plotting2DField(cx,cy,sx,P,["Sx Stress Field","[Pa]"])
    return sx,sy,sxy
