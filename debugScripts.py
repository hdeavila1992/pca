import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import linearElastoStaticsSolver as linElastStat


print('debugScripts.py RUN!!')

def checkingSymmetricMatrix(A):
    check = np.allclose(A, A.T, rtol=1e-2, atol=1e-3)
    if check:
        print("The given matrix is symmetric")
    else:
        print("The given matrix is not symmetric")

def checkFullRankMatrix(A):
    mRank = np.linalg.matrix_rank(A)
    mRows = A.shape[0]
    print("Number of rows: ",mRows)
    print("Rank of matrix: ",mRank)
    if mRank == mRows:
        return True
    else:
        return False

def plotSparsity(A):
    fig = plt.figure()
    plt.spy(A,markersize=5)
    # plt.imshow(A,cmap=cm.viridis)
    # plt.colorbar()
    plt.show()

################ AREA AND LENGTH INTEGRALS ####################

def elementArea(U,V,w,p,q,px,py,gausspoints,gaussweights,paramgrad,apt,cpt):
    elemA = 0
    for qj in range(len(gausspoints)):
        for qi in range(len(gausspoints)):
            coor = linElastStat.parametricCoordinate(apt[0],cpt[0],apt[1],cpt[1],gausspoints[qi],gausspoints[qj])
            jac = linElastStat.jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
            wJac = linElastStat.weightedJacobian(jac,paramgrad,gaussweights,qi,qj)
            elemA += 1.0*wJac

    return elemA

def elementLength(U,V,w,p,q,px,py,gausspoints,gaussweights,apt,bpt,paramside):
    elemL = 0

    if paramside == 1 or paramside == 3:
        paramaxis = 1
    else:
        paramaxis = 0

    for qj in range(len(gausspoints)):
        #The first gausspoints does not influence in the output due to uval
        coor = linElastStat.parametricCoordinate(apt[0],bpt[0],apt[1],bpt[1],gausspoints[qj],gausspoints[qj])
        gcoor = linElastStat.geometricCoordinate(coor,U,V,w,p,q,px,py)
        # print("Geometric coor")
        # print(gcoor)

        Jac = linElastStat.jacobian(U,V,w,p,q,coor[0][0],coor[0][1],px,py)
        jac1 = np.linalg.norm(Jac[:,paramaxis])
        jac2 = 0.5*np.sum(bpt-apt)

        elemL += 1.0*jac1*jac2*gaussweights[qj]

    return elemL

def calculateAreaAndLength(U,V,w,p,q,P,paramnodes,nodeselem,gaussquad,loadelems,loadfaces):
    totalArea = 0
    totalLength = 0
    gaussLegendrePoints = gaussquad[0]
    gaussLegendreWeights = gaussquad[1]

    paramGrad = np.zeros((2,2))
    numElems = nodeselem.shape[0]

    px = np.reshape(P[:,0],(P.shape[0],1))
    py = np.reshape(P[:,1],(P.shape[0],1))

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
        totalArea += elementArea(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,paramGrad,aPoint,cPoint)
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

            uB = paramnodes[nodeselem[ielem][endindex]][0]
            uA = paramnodes[nodeselem[ielem][startindex]][0]
            vB = paramnodes[nodeselem[ielem][endindex]][1]
            vA = paramnodes[nodeselem[ielem][startindex]][1]

            aPoint = np.array([uA,vA])
            bPoint = np.array([uB,vB])

            #totalLength +=    elementLength(U,V,w,p,q,px,py,gaussLegendrePoints,gaussLegendreWeights,aPoint,bPoint,paramside)

        print("---")

    # print("Total Length")
    # print(totalLength)
    # print("Total Area")
    # print(totalArea)

    return totalArea #,totalLength
