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

############# SINGLE KNOT INSERTION FOR CURVE #################

def knotInsertion(U,p,Pw,unew):

    k = findKnotIndex(U,unew)
    r = 1 #Single insertion
    s = 0 #Multiplicity of the knots

    nP = (len(U) - 1) - p - 1
    mp = nP + p + 1
    nq = nP + r

    Qw = np.zeros((Pw.shape[0] + 1,Pw.shape[1]))
    Unew = np.zeros(len(U) + 1)
    Rw = np.zeros((p+1,Pw.shape[1]))

    #Load new knot vector
    for i in range(0,k+1):
        Unew[i] = U[i]

    for i in range(1,r+1):
        Unew[k+i] = unew

    for i in range(k+1,mp+1):
        Unew[i+r] = U[i]

    #Save unaltered control points
    for i in range(0,k-p+1):
        Qw[i,:] = Pw[i,:]

    for i in range(k-s,nP+1):
        Qw[i+r,:] = Pw[i,:]

    for i in range(0,p-s+1):
        Rw[i,:] = Pw[k-p+i,:]

    #Insert the knot r times
    for j in range(1,r+1):
        L = k-p+j

        for i in range(0,p-j-s+1):
            alpha = (unew - U[L+i])/(U[i+k+1]- U[L+i])
            Rw[i,:] = alpha*Rw[i+1,:] + (1.0 - alpha)*Rw[i,:]

        Qw[L,:] = Rw[0,:]
        Qw[k+r-j-s,:] = Rw[p-j-s,:]

    for i in range(L+1,k-s):
        Qw[i,:] = Rw[i-L,:]

    return Unew,Qw

################ KNOT REFINEMENT FOR CURVE ####################

def knotRefinement(U,X,p,Pw):
    a = findKnotIndex(U,X[0])
    b = findKnotIndex(U,X[-1])
    b += 1

    n = Pw.shape[0] - 1
    r = len(X) - 1
    m = n + p + 1
    Qw = np.zeros((n + 1 + r + 1,Pw.shape[1]))
    Ubar = np.zeros(len(U) + r + 1)

    for j in range(0,a-p+1):
        Qw[j,:] = Pw[j,:]

    for j in range(b-1,n+1):
        Qw[j+r+1,:] = Pw[j,:]

    for j in range(0,a+1):
        Ubar[j] = U[j]

    for j in range(b+p,m+1):
        Ubar[j+r+1] = U[j]

    i = b + p - 1
    k = b + p + r
    for j in range(r,0-1,-1):

        while X[j] <= U[i] and i > a:
            Qw[k-p-1,:] = Pw[i-p-1,:]
            Ubar[k] = U[i]
            k -= 1
            i -= 1

        Qw[k-p-1,:] = Qw[k-p,:]
        for l in range(1,p+1):
            ind = k - p + l
            alpha = Ubar[k+l] - X[j]
            if abs(alpha) < 1e-5:
                Qw[ind-1,:] = Qw[ind,:]
            else:
                alpha /= (Ubar[k+l] - U[i-p+l])
                Qw[ind-1,:] = alpha*Qw[ind-1,:] + (1.0 - alpha)*Qw[ind,:]

        Ubar[k] = X[j]
        k -= 1

    return Qw,Ubar

################ Pre-SPLINE DECOMPOSITION FOR CURVE ####################

def preSplineDecomposition(U,p,Pw):
    X = np.array([])
    multiVec = np.ones(p)
    for u in U:
        if abs(u-U.min()) > 1e-5 and abs(u-U.max()) > 1e-5:
            X = np.concatenate([X,u*multiVec])

    if len(X) > 0:
        Qsplit,Usplit = knotRefinement(U,X,p,Pw)
    else:
        Qsplit = Pw
        Usplit = U
        print("The spline only has one element")

    return Qsplit,Usplit

################ SPLINE DECOMPOSITION FOR CURVE ####################

def splineSplitting(U,p,Pw):
    n = Pw.shape[0] - 1
    m = n + p + 1
    a = p
    b = p + 1
    nb = 0
    qList = []
    numSegments = len(np.unique(U)) - 1
    alphas = np.zeros(p-1)

    for i in range(0,numSegments):
        qList.append(np.zeros((p+1,Pw.shape[1])))

    for i in range(0,p+1):
        qList[nb][i,:] = Pw[i,:]

    while b < m:
        i = b
        while b < m and abs(U[b+1] - U[b]) < 1e-5:
            b += 1

        mult = b - i + 1
        if mult < p:
            numer = U[b] - U[a]

            for j in range(p,mult,-1):
                alphas[j-mult-1] = numer/(U[a+j] - U[a])

            r = p - mult
            for j in range(1,r+1):
                save = r - j
                s = mult + j

                for k in range(p,s-1,-1):
                    alpha = alphas[k-s]
                    qList[nb][k,:] = alpha*qList[nb][k,:] + (1.0 - alpha)*qList[nb][k-1,:]

                if b < m:
                    qList[nb+1][save,:] = qList[nb][p,:]

        nb += 1
        if b < m:
            for i in range(p-mult,p+1):
                qList[nb][i,:] = Pw[b-p+i,:]

            a = b
            b += 1

    Qx = np.array([])
    Qy = np.array([])
    for q in qList:
        Qx = np.concatenate([Qx,q[:,0]])
        Qy = np.concatenate([Qy,q[:,1]])
    Q = np.vstack((Qx,Qy))
    return Q

def splineSplittingv2(U,p,Pw):
    n = len(Pw.shape[0]) - 1
    m = n + p + 1
    a = p
    b = p + 1
    alphas = np.zeros(p-1)

    Qw = np.zeros((p+1,Pw.shape[1]))
    NextQw = np.zeros((p-1,Pw.shape[1]))

    for i in range(0,p+1):
        Qw[i,:] = Pw[i,:]

    while b < m:
        i = b
        while b < m and abs(U[b+1] - U[b]) < 1e-5:
            b += 1

        mult = b - i + 1
        if mult < p:
            #numerator of alpha
            numer = U[b] - U[a]

            #Compute and store alpha
            for j in range(p,mult,-1):
                alphas[j-mult-1] = numer/(U[a+j] - U[a])

            r = p - mult
            #Insert knot r times
            for j in range(1,r+1):
                save = r - j
                s = mult + j #This many new points

                for k in range(p,s-1,-1):
                    alpha = alphas[k-s]
                    Qw[k,:] = alpha*Qw[k,:] + (1.0 - alpha)*Qw[k-1,:]

                if b < m:
                    #Control point of the next segment
                    NextQw[save,:] = Qw[p,:]

        print('Bezier control points of the segment')
        print(Qw)
        if b < m:
            #Initialize for next segment
            for i in range(0,p-1):
                Qw[i,:] = NextQw[i,:]

            for i in range(p-mult,p+1):
                Qw[i,:] = Pw[b-p+i,:]

            a = b
            b += 1

################ DEGREE ELEVATION FOR CURVE ####################

def binomialCoefficient(a,b):
    bc = 1.0
    for j in range(1,b+1):
        bc *= ((a+1-j)/j)

    return bc

def degreeElevation(U,p,Pw,t):
    n = Pw.shape[0] - 1
    m = n + p + 1
    ph = p + t
    ph2 = ph//2

    s0 = len(np.unique(U))-2

    bezalfs = np.zeros((p+t+1,p+1))
    bpts = np.zeros((p+1,Pw.shape[1]))
    ebpts = np.zeros((p+t+1,Pw.shape[1]))
    Nextbpts = np.zeros((p-1,Pw.shape[1]))
    alfs = np.zeros(p-1)
    Uh = np.zeros(len(U) + len(np.unique(U))*t)
    Qw = np.zeros((len(Uh) - (p+t) - 1,Pw.shape[1]))

    #Compute bezier degree elevation coefficients
    bezalfs[0][0] = 1.0
    bezalfs[ph][p] = 1.0

    for i in range(1,ph2+1):
        inv = 1.0/binomialCoefficient(ph,i)
        mpi = min(p,i)
        for j in range(max(0,i-t),mpi+1):
            bezalfs[i][j] = inv*binomialCoefficient(p,j)*binomialCoefficient(t,i-j)

    for i in range(ph2+1,ph):
        mpi = min(p,i)
        for j in range(max(0,i-t),mpi+1):
            bezalfs[i][j] = bezalfs[ph-i][p-j]

    mh = ph
    kind = ph + 1
    r = -1
    a = p
    b = p + 1
    cind = 1
    ua = Uh[0]
    Qw[0,:] = Pw[0,:]

    for i in range(0,ph+1):
        Uh[i] = ua

    #Initialize first bezier segment
    for i in range(0,p+1):
        bpts[i,:] = Pw[i,:]

    #Big loop through knot vector
    while b < m:
        i = b
        while b < m and abs(U[b+1] - U[b]) < 1e-5:
            b += 1

        mul = b - i + 1
        mh += mul + t
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
                    bpts[k,:] = alfs[k-s]*bpts[k,:] + (1.0 - alfs[k-s])*bpts[k-1,:]

                Nextbpts[save,:] = bpts[p,:]

        #End of insert knot
        #Degree elevate bezier
        #Only points lbz,...,ph are used below
        # ebpts = np.zeros((2,p+t+1))
        for i in range(lbz,ph+1):
            ebpts[i,:] = np.zeros(Pw.shape[1])
            mpi = min(p,i)
            for j in range(max(0,i-t),mpi+1):
                ebpts[i,:] += bezalfs[i][j]*bpts[j,:]

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
                        Qw[i,:] = alf*Qw[i,:] + (1.0 - alf)*Qw[i-1,:]

                    if j >= lbz:
                        if (j - tr) <= (kind - ph + oldr):
                            gam = (ub - Uh[j - tr])/den
                            ebpts[kj,:] = gam*ebpts[kj,:] + (1.0 - gam)*ebpts[kj+1,:]
                        else:
                            ebpts[kj,:] = bet*ebpts[kj,:] + (1.0 - bet)*ebpts[kj+1,:]

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
            Qw[cind,:] = ebpts[j,:]
            cind += 1

        if b < m:
            #Set up for the next loop through
            for j in range(0,r):
                bpts[j,:] = Nextbpts[j,:]

            for j in range(r,p+1):
                bpts[j,:] = Pw[b-p+j,:]

            a = b
            b += 1
            ua = ub
        else:
            #End knot
            for i in range(0,ph+1):
                Uh[kind+i] = ub

    #End of while loop (b < m)
    nh = mh - ph - 1

    return Qw,Uh,ph

################ GENERAL REFINEMENT ####################

#Knot refinement
def hRefinement(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)

    Ured = U[p:-p]
    X = 0.5*(Ured[0:-1] + Ured[1:])

    Qwh,Uh = knotRefinement(U,X,p,Pw)
    Ph,wh = rbs.geometricControlPoints(Qwh)

    ph = p

    return Uh,ph,Ph,wh

#Degree Elevation
def pRefinement(U,p,P,w):
    Pw = rbs.weightedControlPoints(P,w)
    t = 1
    Qp,Up,pp = degreeElevation(U,p,Pw,t)
    Pp,wp = rbs.geometricControlPoints(Qp)

    return Up,pp,Pp,wp

def kRefinement(U,p,P,w):
    Up,pp,Pp,wp = pRefinement(U,p,P,w)
    Uk,pk,Pk,wk = hRefinement(Up,pp,Pp,wp)

    return Uk,pk,Pk,wk

def curveRefinement(refinementlist,U,p,P,w):
    href = 0
    pref = 0
    kref = 0

    Uin = U
    pin = p
    Pin = P
    win = w

    for rfnlist in refinementlist:
        if rfnlist == 'h':
            Uout,pout,Pout,wout = hRefinement(Uin,pin,Pin,win)
            href += 1
        elif rfnlist == 'p':
            Uout,pout,Pout,wout = pRefinement(Uin,pin,Pin,win)
            pref += 1
        elif rfnlist == 'k':
            Uout,pout,Pout,wout = kRefinement(Uin,pin,Pin,win)
            kref += 1
        else:
            print("Invalid option")

        Uin = Uout
        pin = pout
        Pin = Pout
        win = wout

    print('Number of h-refinements')
    print(href)
    print('Number of p-refinements')
    print(pref)
    print('Number of k-refinements')
    print(kref)

    return Uout,pout,Pout,wout
