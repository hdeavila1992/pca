# Python libraries
import numpy as np
import numpy.linalg

# Local project
import src.basisFunctions as bfunc
import src.nurbs as rbs

def parametricCoordinate(ua,ub,va,vb,gausspta,gaussptb):
    localpts = np.zeros((1,2))
    localpts[0][0] = 0.5*(ub - ua)*gausspta + 0.5*(ub + ua)
    localpts[0][1] = 0.5*(vb - va)*gaussptb + 0.5*(vb + va)
    return localpts

def elasticMatrix(E,nu):
    dmat = np.zeros((3,3))
    dmat[0][0] = 1 - nu
    dmat[1][1] = 1 - nu
    dmat[2][2] = (1 - 2*nu)/2
    dmat[0][1] = nu
    dmat[1][0] = nu
    dmat *= E/((1+nu)*(1-2*nu))
    return dmat

################ ISOGEOMETRIC ANALYSIS ####################

def assemblyWeakForm(U,V,w,p,q,P,surfaceprep,numquad,matprop,boundaryprep,neumannconditions):
    K = np.zeros((2*P.shape[0],2*P.shape[0]))
    F = np.zeros((2*P.shape[0],1))
    Fb = np.zeros((2*P.shape[0],1))
    Fl = np.zeros((2*P.shape[0],1))
    
    numquad2d = numquad[0]
    numquad1d = numquad[1]
    
    # Extraction of surface preprocessing
    nonzeroctrlpts = surfaceprep[0]
    surfacespan = surfaceprep[1]
    elementcorners = surfaceprep[2]
    
    # Extraction of boundary preprocessing
    nonzeroctrlptsload = boundaryprep[0]
    boundaryspan = boundaryprep[1]
    boundarycorners = boundaryprep[2]
    axisselector = boundaryprep[3]

    paramGrad = np.zeros((2,2))
    numElems = len(elementcorners)
    numLoadedElems = len(boundarycorners)

    Pwl = rbs.weightedControlPoints(P,w)
    Pw = rbs.listToGridControlPoints(Pwl,U,V,p,q)

    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])

    loadtype = neumannconditions[0][2]
    loadvalue = neumannconditions[0][3]
    
    # Definition of the material matrix
    E = matprop[0]
    nu = matprop[1]
    rho = matprop[2]
    dMat = elasticMatrix(E,nu)
    
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8
    
    # Precomputing info for the nurbs derivatives
    mU = len(U) - 1
    mV = len(V) - 1
    nU = mU - p - 1
    nV = mV - q - 1
    
    # Strain-energy and body force integrals
    print('Computing the strain-energy and body forces integrals')
    for ielem in range(0,numElems):
        # Extracting the indices of the non-zero control points
        idR = nonzeroctrlpts[ielem]
        # Extracting the indices for the location of the parametric element
        uspan = surfacespan[ielem][0]
        vspan = surfacespan[ielem][1]
        # Extracting the corners of the parametric element
        aPoint = elementcorners[ielem][0]
        cPoint = elementcorners[ielem][1]
        
        # Computing the parametric gradient and its determinant
        paramGrad[0][0] = 0.5*(cPoint[0] - aPoint[0])
        paramGrad[1][1] = 0.5*(cPoint[1] - aPoint[1])
        detJac2 = abs(np.linalg.det(paramGrad))
        
        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1
            
        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')

        # K stiffness matrix
        for iquad in range(numquad2d.shape[0]):
            coor = parametricCoordinate(aPoint[0],cPoint[0],aPoint[1],cPoint[1],numquad2d[iquad][0],numquad2d[iquad][1])
            
            # NURBS gradient
            biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)
            
            # Jacobian
            jac = (biRatGrad[1:3,:]@P[idR,:]).T
            wJac = abs(np.linalg.det(jac))*detJac2*numquad2d[iquad][2]
            
            # Strain displacement matrix
            invJac = np.linalg.inv(jac)
            dN2 = biRatGrad[1:3,:]
            dN2dxi = invJac.T@dN2
            
            bMat = np.zeros((3,2*dN2dxi.shape[1]))
            #dNx
            bMat[0,0::2] = dN2dxi[0,:]
            bMat[2,0::2] = dN2dxi[1,:]
            #dNy
            bMat[1,1::2] = dN2dxi[1,:]
            bMat[2,1::2] = dN2dxi[0,:]
            
            # Global indexing
            K[globalDOFx,globalDOFy] += (bMat.T@dMat@bMat)*wJac
        
            # Body forces integral
            if abs(rho) > 1e-5:
                nMat = np.zeros((2,2*biRatGrad.shape[1]))
                nMat[0,0::2] = biRatGrad[0,:]
                nMat[1,1::2] = biRatGrad[0,:]
                
                Fb[globalDOF] += (nMat.T@bvec)*wJac
    
    # Load integrals
    print('Computing the load integrals')
    for iload in range(0,numLoadedElems):
        # Extracting the indices of the non-zero control points of the loaded elements
        idR = nonzeroctrlptsload[iload]
        # Extracting the indices for the location of the parametric segment
        uspan = boundaryspan[iload][0]
        vspan = boundaryspan[iload][1]
        # Extracting the corners of the parametric segment
        aPoint = boundarycorners[iload][0]
        bPoint = boundarycorners[iload][1]
        # Extracting the non-zero column index of the boundary jacobian
        paramaxis = axisselector[iload]
        
        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1
        
        for iquad in range(numquad1d.shape[0]):
            coor = parametricCoordinate(aPoint[0],bPoint[0],aPoint[1],bPoint[1],numquad1d[iquad][0],numquad1d[iquad][0])

            biRatGrad = rbs.bivariateRationalGradient(mU,mV,p,q,uspan,vspan,coor[0][0],coor[0][1],U,V,Pw)
            Jac = (biRatGrad[1:3,:]@P[idR,:]).T
            jac1 = np.linalg.norm(Jac[:,paramaxis])
            jac2 = 0.5*np.sum(bPoint-aPoint)

            if jac1 > 1e-6:
                unitTangetVec = Jac[:,paramaxis]/jac1
            else:
                unitTangetVec = np.zeros((2,1))

            if loadtype == "tangent":
#                tvec = (loadvalue/abs(loadvalue))*unitTangetVec
                tvec = loadvalue*unitTangetVec
            elif loadtype == "normal":
                unitNormalVec = rotMat@unitTangetVec
#                tvec = (loadvalue/abs(loadvalue))*unitNormalVec
                tvec = loadvalue*unitNormalVec
            else:
                print("Wrong load configuration")

            tvec = np.reshape(tvec,(2,1))
            nMat = np.zeros((2,2*biRatGrad.shape[1]))

            nMat[0,0::2] = biRatGrad[0,:]
            nMat[1,1::2] = biRatGrad[0,:]
            
            Fl[globalDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]


    F = Fb + Fl
    return K,F

def assemblyMultipatchWeakForm(mulU,mulV,fullw,mulp,mulq,fullP,idctrlpts,surfaceprep,numquad,matprop,boundaryprep):
    Ktotal = np.zeros((2*fullP.shape[0],2*fullP.shape[0]))
    Ftotal = np.zeros((2*fullP.shape[0],1))
    Fbtotal = np.zeros((2*fullP.shape[0],1))
    Fltotal = np.zeros((2*fullP.shape[0],1))
    
    numquad2d = numquad[0]
    numquad1d = numquad[1]
    
    # Definition of the material matrix
    E = matprop[0]
    nu = matprop[1]
    rho = matprop[2]
    dMat = elasticMatrix(E,nu)
    
    # Rotation matrix for -pi/2
    rotMat = np.array([[0.0,1.0],[-1.0,0.0]])
    
    paramGrad = np.zeros((2,2))
    
    bvec = np.zeros((2,1))
    bvec[1][0] = -rho*9.8
    
    numpatches = len(mulU)
    
    # Patch loop
    print('Computing the strain-energy and body forces integrals')
    for ipatch in range(0,numpatches):
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]
        
        pi = mulp[ipatch]
        qi = mulq[ipatch]
        
        Pi = fullP[idctrlpts[ipatch],:]
        wi = fullw[idctrlpts[ipatch],:]
        
        Kpatch = np.zeros((2*Pi.shape[0],2*Pi.shape[0]))
        Fbpatch = np.zeros((2*Pi.shape[0],1))
        
        # Extraction of surface preprocessing
        nonzeroctrlpts = surfaceprep[ipatch][0]
        surfacespan = surfaceprep[ipatch][1]
        elementcorners = surfaceprep[ipatch][2]
        numElems = len(elementcorners)
        
        Pwl = rbs.weightedControlPoints(Pi,wi)
        Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)
        
        # Precomputing info for the nurbs derivatives
        mU = len(Ui) - 1
        mV = len(Vi) - 1
        nU = mU - pi - 1
        nV = mV - qi - 1
        
        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idctrlpts[ipatch]),dtype=int)
        dof0 = 2*np.array(idctrlpts[ipatch])
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1
        
        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')
        
        # Strain-energy and body force integrals
        for ielem in range(0,numElems):
            # Extracting the indices of the non-zero control points
            idR = nonzeroctrlpts[ielem]
            # Extracting the indices for the location of the parametric element
            uspan = surfacespan[ielem][0]
            vspan = surfacespan[ielem][1]
            # Extracting the corners of the parametric element
            aPoint = elementcorners[ielem][0]
            cPoint = elementcorners[ielem][1]
            
            # Computing the parametric gradient and its determinant
            paramGrad[0][0] = 0.5*(cPoint[0] - aPoint[0])
            paramGrad[1][1] = 0.5*(cPoint[1] - aPoint[1])
            detJac2 = abs(np.linalg.det(paramGrad))
            
            # Patch degrees of freedom
            patchDOF = np.zeros(2*len(idR),dtype=int)
            dof0 = 2*np.array(idR)
            dof1 = dof0 + 1
            patchDOF[0::2] = dof0
            patchDOF[1::2] = dof1
                
            patchDOFx,patchDOFy = np.meshgrid(patchDOF,patchDOF,indexing='xy')

            # K stiffness matrix
            for iquad in range(numquad2d.shape[0]):
                coor = parametricCoordinate(aPoint[0],cPoint[0],aPoint[1],cPoint[1],numquad2d[iquad][0],numquad2d[iquad][1])
                
                # NURBS gradient
                biRatGrad = rbs.bivariateRationalGradient(mU,mV,pi,qi,uspan,vspan,coor[0][0],coor[0][1],Ui,Vi,Pwi)
                
                # Jacobian
                jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
                wJac = abs(np.linalg.det(jac))*detJac2*numquad2d[iquad][2]
                
                # Strain displacement matrix
                invJac = np.linalg.inv(jac)
                dN2 = biRatGrad[1:3,:]
                dN2dxi = invJac.T@dN2
                
                bMat = np.zeros((3,2*dN2dxi.shape[1]))
                #dNx
                bMat[0,0::2] = dN2dxi[0,:]
                bMat[2,0::2] = dN2dxi[1,:]
                #dNy
                bMat[1,1::2] = dN2dxi[1,:]
                bMat[2,1::2] = dN2dxi[0,:]
                
                # Patch indexing
                Kpatch[patchDOFx,patchDOFy] += (bMat.T@dMat@bMat)*wJac
            
                # Body forces integral
                if abs(rho) > 1e-5:
                    nMat = np.zeros((2,2*biRatGrad.shape[1]))
                    nMat[0,0::2] = biRatGrad[0,:]
                    nMat[1,1::2] = biRatGrad[0,:]
                    
                    Fbpatch[patchDOF] += (nMat.T@bvec)*wJac
            # End of quadrature loop
        # End of element loop
        Ktotal[globalDOFx,globalDOFy] += Kpatch
        Fbtotal[globalDOF] += Fbpatch
    # End of patch loop
    
    # Load conditions loop
    # Extraction of boundary preprocessing
    loadedpatches = boundaryprep[0]
    nonzeroctrlptsload = boundaryprep[1]
    boundaryspan = boundaryprep[2]
    boundarycorners = boundaryprep[3]
    axisselector = boundaryprep[4]
    valuesload = boundaryprep[5]
    loadtype = boundaryprep[6]

    numLoadedElems = len(boundarycorners)

    # Load integrals
    print('Computing the load integrals')
    for iload in range(0,numLoadedElems):
        # Extracting the indices of the patch
        ipatch = loadedpatches[iload]
        # Extracting the indices of the non-zero control points of the loaded elements
        idR = nonzeroctrlptsload[iload]
        # Extracting the indices for the location of the parametric segment
        uspan = boundaryspan[iload][0]
        vspan = boundaryspan[iload][1]
        # Extracting the corners of the parametric segment
        aPoint = boundarycorners[iload][0]
        bPoint = boundarycorners[iload][1]
        # Extracting the non-zero column index of the boundary jacobian
        paramaxis = axisselector[iload]
        # Extracting the value of the load in the face
        load = valuesload[iload]
        
        Ui = mulU[ipatch]
        Vi = mulV[ipatch]
        
        pi = mulp[ipatch]
        qi = mulq[ipatch]
        
        Pi = fullP[idctrlpts[ipatch],:]
        wi = fullw[idctrlpts[ipatch],:]
        
        Flpatch = np.zeros((2*Pi.shape[0],1))
        
        mu = len(Ui) - 1
        mv = len(Vi) - 1
        
        Pwl = rbs.weightedControlPoints(Pi,wi)
        Pwi = rbs.listToGridControlPoints(Pwl,Ui,Vi,pi,qi)
        
        # Global degrees of freedom
        globalDOF = np.zeros(2*len(idctrlpts[ipatch]),dtype=int)
        dof0 = 2*np.array(idctrlpts[ipatch])
        dof1 = dof0 + 1
        globalDOF[0::2] = dof0
        globalDOF[1::2] = dof1
        
        globalDOFx,globalDOFy = np.meshgrid(globalDOF,globalDOF,indexing='xy')
        
        # Patch degrees of freedom
        patchDOF = np.zeros(2*len(idR),dtype=int)
        dof0 = 2*np.array(idR)
        dof1 = dof0 + 1
        patchDOF[0::2] = dof0
        patchDOF[1::2] = dof1
        
        patchDOFx,patchDOFy = np.meshgrid(patchDOF,patchDOF,indexing='xy')
        
        for iquad in range(numquad1d.shape[0]):
            coor = parametricCoordinate(aPoint[0],bPoint[0],aPoint[1],bPoint[1],numquad1d[iquad][0],numquad1d[iquad][0])

            biRatGrad = rbs.bivariateRationalGradient(mU,mV,pi,qi,uspan,vspan,coor[0][0],coor[0][1],Ui,Vi,Pwi)
            Jac = (biRatGrad[1:3,:]@Pi[idR,:]).T
            jac1 = np.linalg.norm(Jac[:,paramaxis])
            jac2 = 0.5*np.sum(bPoint-aPoint)

            if jac1 > 1e-6:
                unitTangetVec = Jac[:,paramaxis]/jac1
            else:
                unitTangetVec = np.zeros((2,1))

            if loadtype[iload] == "tangent":
#                tvec = (loadvalue/abs(loadvalue))*unitTangetVec
                tvec = load*unitTangetVec
            elif loadtype[iload] == "normal":
                unitNormalVec = rotMat@unitTangetVec
#                tvec = (loadvalue/abs(loadvalue))*unitNormalVec
                tvec = load*unitNormalVec
            else:
                print("Wrong load configuration")

            tvec = np.reshape(tvec,(2,1))
            nMat = np.zeros((2,2*biRatGrad.shape[1]))

            nMat[0,0::2] = biRatGrad[0,:]
            nMat[1,1::2] = biRatGrad[0,:]
            
            Flpatch[patchDOF] += (nMat.T@tvec)*jac1*jac2*numquad1d[iquad][1]
        # End quadrature loop
        Fltotal[globalDOF] += Flpatch
    # End load loop

    Ftotal = Fbtotal + Fltotal
    return Ktotal,Ftotal