import numpy as np
import nurbs as rbs

tol = 1e-5

def findKnotIndex(U,u):
    for k in range(len(U)):
        if U[k] < (u + tol) and (u + tol) < (U[k+1]):
            kindex = k

        if abs(u - U.max())<tol:
            if U[k] < (u - tol) and (u - tol) < U[k+1]:
                kindex = k
    return kindex

############# SINGLE KNOT INSERTION FOR SURFACE #################

def knotInsertion(unew,vnew,dir,U,V,p,q,Pwg):
    if dir == "UDIR":
        Unew = np.zeros(len(U) + 1)
        Vnew = np.zeros(len(V) + 1)

        k = findKnotIndex(U,unew)
        r = 1 #Single insertion
        s = 0 #Multiplicity of the knots

        # nP is the number of rows of control points minus 1
        # before the knot insertion
        # Np is the number of rows of control points
        # before the knot insertion

        # mp is the number of columns of control points minus 1
        # before the knot insertion

        nP = Pwg.shape[1] - 1
        Np = nP + 1
        mP = nP + p + 1
        mp = Pwg.shape[2] - 1

        # nq is the number of rows of control points minus one
        # after the knot insertion

        nq = nP + r
        Nq = nq + 1

        Qwg = np.zeros((Pwg.shape[0],Nq,Pwg.shape[2]))
        Rw = np.zeros((p+1,Pwg.shape[0]))
        alpha = np.zeros((p-r-s+1,r+1))

        #Load new vector
        for i in range(0,k+1):
            Unew[i] = U[i]

        for i in range(1,r+1):
            Unew[k+i] = unew

        for i in range(k+1,mP+1):
            Unew[i+r] = U[i]

        #Copy old vector in new vector
        Vnew = V

        #Save the alphas
        for j in range(1,r+1):
            L = k - p + j

            for i in range(0,p-j-s+1):
                alpha[i][j] = (unew - U[L+i])/(U[i+k+1] - U[L+i])

        #For each row do
        for row in range(0,mp+1):

            #Save unaltered control points
            for i in range(0,k-p+1):
                Qwg[:,i,row] = Pwg[:,i,row]

            for i in range(k-s,nP+1):
                Qwg[:,i+r,row] = Pwg[:,i,row]

            #Load auxiliary control points
            for i in range(0,p-s+1):
                Rw[i,:] = Pwg[:,k-p+i,row]

            #Insert the knot r times
            for j in range(1,r+1):
                L = k-p+j

                for i in range(0,p-j-s+1):
                    Rw[i,:] = alpha[i][j]*Rw[i+1,:] + (1.0 - alpha[i][j])*Rw[i,:]

                Qwg[:,L,row] = Rw[0,:]
                Qwg[:,k+r-j-s,row] = Rw[p-j-s,:]

            #Load remaining control points
            for i in range(L+1,k-s):
                Qwg[:,i,row] = Rw[i-L,:]

    if dir == "VDIR":
        Unew = np.zeros(len(U) + 1)
        Vnew = np.zeros(len(V) + 1)

        k = findKnotIndex(V,vnew)
        r = 1 #Single insertion
        s = 0 #Multiplicity of the knots

        # nP is the number of columns of control points minus 1
        # before the knot insertion
        # Np is the number of colums of control points
        # before the knot insertion

        # mp is the number of rows of control points minus 1
        # before the knot insertion

        nP = Pwg.shape[2] - 1
        Np = nP + 1
        mP = nP + q + 1
        mp = Pwg.shape[1] - 1

        # nq is the number of columns of control points minus one
        # after the knot insertion

        nq = nP + r
        Nq = nq + 1

        Qwg = np.zeros((Pwg.shape[0],Pwg.shape[1],Nq))
        Rw = np.zeros((q+1,Pwg.shape[0]))
        alpha = np.zeros((q-r-s+1,r+1))

        #Load new vector
        for i in range(0,k+1):
            Vnew[i] = V[i]

        for i in range(1,r+1):
            Vnew[k+i] = vnew

        for i in range(k+1,mP+1):
            Vnew[i+r] = V[i]

        #Copy old vector in new vector
        Unew = U

        #Save the alphas
        for j in range(1,r+1):
            L = k - q + j

            for i in range(0,q-j-s+1):
                alpha[i][j] = (vnew - V[L+i])/(V[i+k+1] - V[L+i])

        #For each row do
        for col in range(0,mp+1):

            #Save unaltered control points
            for i in range(0,k-q+1):
                Qwg[:,col,i] = Pwg[:,col,i]

            for i in range(k-s,nP+1):
                Qwg[:,col,i+r] = Pwg[:,col,i]

            #Load auxiliary control points
            for i in range(0,q-s+1):
                Rw[i,:] = Pwg[:,col,k-q+i]

            #Insert the knot r times
            for j in range(1,r+1):
                L = k-q+j

                for i in range(0,q-j-s+1):
                    Rw[i,:] = alpha[i][j]*Rw[i+1,:] + (1.0 - alpha[i][j])*Rw[i,:]

                Qwg[:,col,L] = Rw[0,:]
                Qwg[:,col,k+r-j-s] = Rw[q-j-s,:]

            #Load remaining control points
            for i in range(L+1,k-s):
                Qwg[:,col,i] = Rw[i-L,:]

    return Unew,Vnew,Qwg

################ KNOT REFINEMENT FOR SURFACE ####################

def knotRefinement(X,dir,U,V,p,q,Pwg):
    if dir == "UDIR":
        a = findKnotIndex(U,X[0])
        b = findKnotIndex(U,X[-1])
        b += 1

        n = Pwg.shape[1] - 1
        r = len(X) - 1
        m = n + p + 1
        mcol = Pwg.shape[2] - 1
        Qwg = np.zeros((Pwg.shape[0],n + 1 + r + 1,Pwg.shape[2]))
        Ubar = np.zeros(len(U) + r + 1)

        #Initializing Ubar
        for j in range(0,a+1):
            Ubar[j] = U[j]

        for j in range(b+p,m+1):
            Ubar[j+r+1] = U[j]

        #Copying V into Vbar
        Vbar = V

        #Save unaltered control points
        for row in range(0,mcol+1):

            for k in range(0,a-p+1):
                Qwg[:,k,row] = Pwg[:,k,row]

            for k in range(b-1,n+1):
                Qwg[:,k+r+1,row] = Pwg[:,k,row]

        i = b + p - 1
        k = b + p + r
        for j in range(r,-1,-1):

            while X[j] <= U[i] and i > a:
                Ubar[k] = U[i]
                for row in range(0,mcol+1):
                    Qwg[:,k-p-1,row] = Pwg[:,i-p-1,row]

                k -= 1
                i -= 1

            for row in range(0,mcol+1):
                Qwg[:,k-p-1,row] = Qwg[:,k-p,row]

            for l in range(1,p+1):
                ind = k - p + l
                alpha = Ubar[k+l] - X[j]
                if abs(alpha) < 1e-5:
                    for row in range(0,mcol+1):
                        Qwg[:,ind-1,row] = Qwg[:,ind,row]
                else:
                    alpha /= (Ubar[k+l] - U[i-p+l])
                    for row in range(0,mcol+1):
                        Qwg[:,ind-1,row] = alpha*Qwg[:,ind-1,row] + (1.0 - alpha)*Qwg[:,ind,row]

            Ubar[k] = X[j]
            k -= 1

    if dir == "VDIR":
        a = findKnotIndex(V,X[0])
        b = findKnotIndex(V,X[-1])
        b += 1

        n = Pwg.shape[2] - 1
        r = len(X) - 1
        m = n + q + 1
        nrow = Pwg.shape[1] - 1
        Qwg = np.zeros((Pwg.shape[0],Pwg.shape[1],n + 1 + r + 1))
        Vbar = np.zeros(len(V) + r + 1)

        #Initializing Ubar
        for j in range(0,a+1):
            Vbar[j] = V[j]

        for j in range(b+q,m+1):
            Vbar[j+r+1] = V[j]

        #Copying U into Ubar
        Ubar = U

        #Save unaltered control points
        for col in range(0,nrow+1):

            for k in range(0,a-q+1):
                Qwg[:,col,k] = Pwg[:,col,k]

            for k in range(b-1,n+1):
                Qwg[:,col,k+r+1] = Pwg[:,col,k]

        i = b + q - 1
        k = b + q + r
        for j in range(r,-1,-1):

            while X[j] <= V[i] and i > a:
                Vbar[k] = V[i]
                for col in range(0,nrow+1):
                    Qwg[:,col,k-q-1] = Pwg[:,col,i-q-1]

                k -= 1
                i -= 1

            for col in range(0,nrow+1):
                Qwg[:,col,k-q-1] = Qwg[:,col,k-q]

            for l in range(1,q+1):
                ind = k - q + l
                alpha = Vbar[k+l] - X[j]
                if abs(alpha) < 1e-5:
                    for col in range(0,nrow+1):
                        Qwg[:,col,ind-1] = Qwg[:,col,ind]
                else:
                    alpha /= (Vbar[k+l] - V[i-q+l])
                    for col in range(0,nrow+1):
                        Qwg[:,col,ind-1] = alpha*Qwg[:,col,ind-1] + (1.0 - alpha)*Qwg[:,col,ind]

            Vbar[k] = X[j]
            k -= 1

    return Ubar,Vbar,Qwg

################ Pre-SPLINE DECOMPOSITION FOR SURFACE ####################

def preSplineDecomposition(refdir,U,V,p,q,Pw):
    X = np.array([])
    multiVec = np.ones(p)
    for u in U:
        if abs(u-U.min()) > 1e-5 and abs(u-U.max()) > 1e-5:
            X = np.concatenate([X,u*multiVec])

    if len(X) > 0:
        # refdir = "UDIR"
        Usplit,Vsplit,Qsplit = knotRefinement(X,refdir,U,V,p,q,Pw)
    else:
        Usplit = U
        Vsplit = V
        Qsplit = Pw
        print("The spline only has one element")

    return Usplit,Vsplit,Qsplit

################ SPLINE DECOMPOSITION FOR CURVE ####################

################ DEGREE ELEVATION FOR SURFACE ####################

def binomialCoefficient(a,b):
    bc = 1.0
    for j in range(1,b+1):
        bc *= ((a+1-j)/j)

    return bc

def surfaceDegreeElevation(dir,U,V,p,q,Pwg,tp,tq):
    if dir == "UDIR":
        n = Pwg.shape[1] - 1 #n is the number of rows of Pwg minus 1
        mN = n + p + 1 #mN is the m value from the U knot vector

        ph = p + tp
        ph2 = ph//2

        s0 = len(np.unique(U))-2

        bezalfs = np.zeros((p+tp+1,p+1))
        bpts = np.zeros((Pwg.shape[0],p+1,Pwg.shape[2]))
        ebpts = np.zeros((Pwg.shape[0],p+tp+1,Pwg.shape[2]))
        Nextbpts = np.zeros((Pwg.shape[0],p-1,Pwg.shape[2]))
        alfs = np.zeros(p-1)
        Uh = np.zeros(len(U) + len(np.unique(U))*tp)
        Qwg = np.zeros((Pwg.shape[0],len(Uh) - (p+tp) - 1,Pwg.shape[2]))

        Vh = V
        qh = q

        #Compute bezier degree elevation coefficients
        bezalfs[0][0] = 1.0
        bezalfs[ph][p] = 1.0

        for i in range(1,ph2+1):
            inv = 1.0/binomialCoefficient(ph,i)
            mpi = min(p,i)
            for j in range(max(0,i-tp),mpi+1):
                bezalfs[i][j] = inv*binomialCoefficient(p,j)*binomialCoefficient(tp,i-j)

        for i in range(ph2+1,ph):
            mpi = min(p,i)
            for j in range(max(0,i-tp),mpi+1):
                bezalfs[i][j] = bezalfs[ph-i][p-j]

        mh = ph
        kind = ph + 1
        r = -1
        a = p
        b = p + 1
        cind = 1
        ua = Uh[0]
        # Qw[0,:] = Pw[0,:]
        Qwg[:,0,:] = Pwg[:,0,:]

        for i in range(0,ph+1):
            Uh[i] = ua

        #Initialize first bezier segment
        for i in range(0,p+1):
            # bpts[i,:] = Pw[i,:]
            bpts[:,i,:] = Pwg[:,i,:]

        #Big loop through knot vector
        while b < mN:
            i = b
            while b < mN and abs(U[b+1] - U[b]) < 1e-5:
                b += 1

            mul = b - i + 1
            mh += mul + tp
            ub = U[b]
            oldr = r
            r = p - mul
            #Insert knot u(b) r times
            if oldr > 0:
                lbz = (oldr + 2)//2
            else:
                lbz = 1

            if r > 0:
                rbz = ph - (r + 1)//2
            else:
                rbz = ph

            #Insert knot to get bezier segment
            if r > 0:
                numer = ub - ua
                for k in range(p,mul,-1):
                    alfs[k-mul-1] = numer/(U[a+k] - ua)

                for j in range(1,r+1):
                    save = r - j
                    s = mul + j
                    for k in range(p,s-1,-1):
                        bpts[:,k,:] = alfs[k-s]*bpts[:,k,:] + (1.0 - alfs[k-s])*bpts[:,k-1,:]

                    Nextbpts[:,save,:] = bpts[:,p,:]

            #End of insert knot
            #Degree elevate bezier
            #Only points lbz,...,ph are used below
            # ebpts = np.zeros((2,p+t+1))
            for i in range(lbz,ph+1):
                ebpts[:,i,:] = np.zeros((Pwg.shape[0],Pwg.shape[2]))
                mpi = min(p,i)
                for j in range(max(0,i-tp),mpi+1):
                    ebpts[:,i,:] += bezalfs[i][j]*bpts[:,j,:]

            # print('Bezier + t control points of current segment')
            # print(ebpts)
            #End of degree elevating bezier
            #Must remove knot u = U[a] oldr times
            if oldr > 1:
                first = kind - 2
                last = kind
                den = ub - ua
                bet = (ub - Uh[kind-1])/den
                #Knot removal loop
                for tr in range(1,oldr):
                    i = first
                    j = last
                    kj = j - kind + 1
                    #Loop and compute the new
                    #Control points for one removal step
                    while (j - i) > tr:
                        if i < cind:
                            alf = (ub - Uh[i])/(ua - Uh[i])
                            Qwg[:,i,:] = alf*Qwg[:,i,:] + (1.0 - alf)*Qwg[:,i-1,:]

                        if j >= lbz:
                            if (j - tr) <= (kind - ph + oldr):
                                gam = (ub - Uh[j - tr])/den
                                ebpts[:,kj,:] = gam*ebpts[:,kj,:] + (1.0 - gam)*ebpts[:,kj+1,:]
                            else:
                                ebpts[:,kj,:] = bet*ebpts[:,kj,:] + (1.0 - bet)*ebpts[:,kj+1,:]

                        i += 1
                        j -= 1
                        kj -= 1

                    first -= 1
                    last += 1

            #End of removing knot u = U[a]
            #Load the knot ua
            if a != p:
                for i in range(0,ph-oldr):
                    Uh[kind] = ua
                    kind += 1

            #Load control points into Qw
            for j in range(lbz,rbz+1):
                Qwg[:,cind,:] = ebpts[:,j,:]
                cind += 1

            if b < mN:
                #Set up for the next loop through
                for j in range(0,r):
                    bpts[:,j,:] = Nextbpts[:,j,:]

                for j in range(r,p+1):
                    bpts[:,j,:] = Pwg[:,b-p+j,:]

                a = b
                b += 1
                ua = ub
            else:
                #End knot
                for i in range(0,ph+1):
                    Uh[kind+i] = ub

        #End of while loop (b < m)
        nh = mh - ph - 1

    if dir == "VDIR":
        m = Pwg.shape[2] - 1 #m is the number of columns of Pwg minus 1
        mM = m + q + 1 #mM is the m value from the V knot vector

        qh = q + tq
        qh2 = qh//2

        s0 = len(np.unique(V))-2

        bezalfs = np.zeros((q+tq+1,q+1))
        bpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q+1))
        ebpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q+tq+1))
        Nextbpts = np.zeros((Pwg.shape[0],Pwg.shape[1],q-1))
        alfs = np.zeros(q-1)
        Vh = np.zeros(len(V) + len(np.unique(V))*tq)
        Qwg = np.zeros((Pwg.shape[0],Pwg.shape[1],len(Vh) - (q+tq) - 1))

        Uh = U
        ph = p

        #Compute bezier degree elevation coefficients
        bezalfs[0][0] = 1.0
        bezalfs[qh][q] = 1.0

        for i in range(1,qh2+1):
            inv = 1.0/binomialCoefficient(qh,i)
            mqi = min(q,i)
            for j in range(max(0,i-tq),mqi+1):
                bezalfs[i][j] = inv*binomialCoefficient(q,j)*binomialCoefficient(tq,i-j)

        for i in range(qh2+1,qh):
            mqi = min(q,i)
            for j in range(max(0,i-tq),mqi+1):
                bezalfs[i][j] = bezalfs[qh-i][q-j]

        mh = qh
        kind = qh + 1
        r = -1
        a = q
        b = q + 1
        cind = 1
        va = Vh[0]
        # Qw[0,:] = Pw[0,:]
        Qwg[:,:,0] = Pwg[:,:,0]

        for i in range(0,qh+1):
            Vh[i] = va

        #Initialize first bezier segment
        for i in range(0,q+1):
            # bpts[i,:] = Pw[i,:]
            bpts[:,:,i] = Pwg[:,:,i]

        #Big loop through knot vector
        while b < mM:
            i = b
            while b < mM and abs(V[b+1] - V[b]) < 1e-5:
                b += 1

            mul = b - i + 1
            mh += mul + tq
            vb = V[b]
            oldr = r
            r = q - mul
            #Insert knot u(b) r times
            if oldr > 0:
                lbz = (oldr + 2)//2
            else:
                lbz = 1

            if r > 0:
                rbz = qh - (r + 1)//2
            else:
                rbz = qh

            #Insert knot to get bezier segment
            if r > 0:
                numer = vb - va
                for k in range(q,mul,-1):
                    alfs[k-mul-1] = numer/(V[a+k] - va)

                for j in range(1,r+1):
                    save = r - j
                    s = mul + j
                    for k in range(q,s-1,-1):
                        bpts[:,:,k] = alfs[k-s]*bpts[:,:,k] + (1.0 - alfs[k-s])*bpts[:,:,k-1]

                    Nextbpts[:,:,save] = bpts[:,:,q]

            #End of insert knot
            #Degree elevate bezier
            #Only points lbz,...,ph are used below
            # ebpts = np.zeros((2,p+t+1))
            for i in range(lbz,qh+1):
                ebpts[:,:,i] = np.zeros((Pwg.shape[0],Pwg.shape[1]))
                mqi = min(q,i)
                for j in range(max(0,i-tq),mqi+1):
                    ebpts[:,:,i] += bezalfs[i][j]*bpts[:,:,j]

            # print('Bezier + t control points of current segment')
            # print(ebpts)
            #End of degree elevating bezier
            #Must remove knot u = U[a] oldr times
            if oldr > 1:
                first = kind - 2
                last = kind
                den = vb - va
                bet = (vb - Vh[kind-1])/den
                #Knot removal loop
                for tr in range(1,oldr):
                    i = first
                    j = last
                    kj = j - kind + 1
                    #Loop and compute the new
                    #Control points for one removal step
                    while (j - i) > tr:
                        if i < cind:
                            alf = (vb - Vh[i])/(va - Vh[i])
                            Qwg[:,:,i] = alf*Qwg[:,:,i] + (1.0 - alf)*Qwg[:,:,i-1]

                        if j >= lbz:
                            if (j - tr) <= (kind - qh + oldr):
                                gam = (vb - Vh[j - tr])/den
                                ebpts[:,:,kj] = gam*ebpts[:,:,kj] + (1.0 - gam)*ebpts[:,:,kj+1]
                            else:
                                ebpts[:,:,kj] = bet*ebpts[:,:,kj] + (1.0 - bet)*ebpts[:,:,kj+1]

                        i += 1
                        j -= 1
                        kj -= 1

                    first -= 1
                    last += 1

            #End of removing knot v = V[a]
            #Load the knot va
            if a != q:
                for i in range(0,qh-oldr):
                    Vh[kind] = va
                    kind += 1

            #Load control points into Qw
            for j in range(lbz,rbz+1):
                Qwg[:,:,cind] = ebpts[:,:,j]
                cind += 1

            if b < mM:
                #Set up for the next loop through
                for j in range(0,r):
                    bpts[:,:,j] = Nextbpts[:,:,j]

                for j in range(r,q+1):
                    bpts[:,:,j] = Pwg[:,:,b-q+j]

                a = b
                b += 1
                va = vb
            else:
                #End knot
                for i in range(0,qh+1):
                    Vh[kind+i] = vb

        #End of while loop (b < m)
        nh = mh - qh - 1

    return Uh,Vh,ph,qh,Qwg

################ GENERAL REFINEMENT ####################

#Knot refinement
def hRefinement(rfndir,U,V,p,q,P,w):
    if rfndir == "U":
        dir = "UDIR"
        Ured = U[p:-p]
        X = 0.5*(Ured[0:-1] + Ured[1:])
    elif rfndir == "V":
        dir = "VDIR"
        Ured = V[q:-q]
        X = 0.5*(Ured[0:-1] + Ured[1:])
    else:
        print("Wrong parameter direction")

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    Uh,Vh,Qwhgrid = knotRefinement(X,dir,U,V,p,q,Pwgrid)
    Qw = rbs.gridToListControlPoints(Qwhgrid)
    Ph,wh = rbs.geometricControlPoints(Qw)

    ph = p
    qh = q

    return Uh,Vh,ph,qh,Ph,wh

#Degree Elevation
def pRefinement(rfndir,U,V,p,q,P,w):
    if rfndir == "U":
        dir = "UDIR"
    elif rfndir == "V":
        dir = "VDIR"
    else:
        print("Wrong parameter direction")

    Pw = rbs.weightedControlPoints(P,w)
    Pwgrid = rbs.listToGridControlPoints(Pw,U,V,p,q)

    tp = 1
    tq = 1
    Up,Vp,pp,qp,Qwpgrid = surfaceDegreeElevation(dir,U,V,p,q,Pwgrid,tp,tq)
    Qwp = rbs.gridToListControlPoints(Qwpgrid)
    Pp,wp = rbs.geometricControlPoints(Qwp)

    return Up,Vp,pp,qp,Pp,wp

def kRefinement(rfndir,U,V,p,q,P,w):

    Up,Vp,pp,qp,Pp,wp = pRefinement(rfndir,U,V,p,q,P,w)
    Uk,Vk,pk,qk,Pk,wk = hRefinement(rfndir,Up,Vp,pp,qp,Pp,wp)

    return Uk,Vk,pk,qk,Pk,wk

def surfaceRefinement(refinementlist,directionlist,U,V,p,q,P,w):
    href = 0
    pref = 0
    kref = 0

    Uin = U
    Vin = V
    pin = p
    qin = q
    Pin = P
    win = w

    numindex = 0

    for rfnlist in refinementlist:
        dirlist = directionlist[numindex]
        # print(dirlist)
        if rfnlist == 'h':
            Uout,Vout,pout,qout,Pout,wout = hRefinement(dirlist,Uin,Vin,pin,qin,Pin,win)
            href += 1
        elif rfnlist == 'p':
            Uout,Vout,pout,qout,Pout,wout = pRefinement(dirlist,Uin,Vin,pin,qin,Pin,win)
            pref += 1
        elif rfnlist == 'k':
            Uout,Vout,pout,qout,Pout,wout = kRefinement(dirlist,Uin,Vin,pin,qin,Pin,win)
            kref += 1
        else:
            print("Invalid option")

        Uin = Uout
        Vin = Vout
        pin = pout
        qin = qout
        Pin = Pout
        win = wout

        numindex += 1

    print('Number of h-refinements')
    print(href)
    print('Number of p-refinements')
    print(pref)
    print('Number of k-refinements')
    print(kref)

    return Uout,Vout,pout,qout,Pout,wout
