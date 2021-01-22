import numpy as np
import pca_01 as pca
import bSplines as bs
import plottingScripts as plts
import matplotlib.pyplot as plt
import refinements as rfn

####################################################
###################MAIN PROGRAM#####################
####################################################

################ KNOT INSERTION EXAMPLE ####################

# P = np.array([[0.0,0.5,1.0],[0.0,1.0,0.0]])
# p = 2
# U = np.array([0,0,0,1,1,1])
# X = np.array([0.5])
#
# Uh,Ph = knotInsertion(U,p,P,0.5)
#
# cx,cy = pca.bSplineCurve(U,p,P)
# cxh,cyh = pca.bSplineCurve(Uh,p,Ph)

# plts.plotCurve2d(cx,cy,P)
# plts.plotCurve2d(cxnew,cynew,Pbar)
# plts.plotKnotInsertion(cx,cy,cxh,cyh,P,Ph)

################ KNOT REFINEMENT EXAMPLE ####################

P = np.array([[0.0,1.0,0.75,4.5,3.0,5.0],[0.0,0.1,2.0,3.0,0.2,0.0]])
p = 3
U = np.array([0,0,0,0,0.3,0.7,1,1,1,1])
X = np.array([0.15,0.5,0.85])

Ured = U[p:-p]
X = 0.5*(Ured[0:-1] + Ured[1:])

Q,Ubar = rfn.knotRefinement(U,X,p,P)

cx,cy = bs.bSplineCurve(U,p,P)
cxref,cyref = bs.bSplineCurve(Ubar,p,Q)
# plts.plotCurve2d(cx,cy,P)
# plts.plotCurve2d(cxref,cyref,Q)
# plts.plotKnotRefinement(cx,cy,cxref,cyref,P,Q)

################ Pre-SPLINE DECOMPOSITION EXAMPLE ####################

# P = np.array([[0,0.1,1.8,2.3,3.8,4.5,3.5],[0.1,1.0,1.1,0.0,0.0,0.7,1.25]])
# p = 3
# U = np.array([0,0,0,0,1,2,3,4,4,4,4])
#
# cx,cy = pca.bSplineCurve(U,p,P)
# # plts.plotCurve2d(cx,cy,P)
#
# Qdecmp,Udecmp = preSplineDecomposition(U,p,P)
# print(Qdecmp)
# print(Udecmp)
#
# cxdecmp,cydecmp = pca.bSplineCurve(Udecmp,p,Qdecmp)
# plts.plotCurve2d(cxdecmp,cydecmp,Qdecmp)

################ SPLINE DECOMPOSITION EXAMPLE ####################
# P = np.array([[0.0,1.0,0.5,2.5,2.0,3.0],[0.0,0.2,1.0,1.0,0.25,0.0]])
# p = 3
# U = np.array([0,0,0,0,0.4,0.7,1,1,1,1])

# Qsplit,Usplit = rfn.preSplineDecomposition(U,p,P)
# Qsplit2 = rfn.splineSplitting(U,p,P)
# rfn.splineSplittingv2(U,p,P)

# cx,cy = pca.bSplineCurve(U,p,P)
# cxs,cys = pca.bSplineCurve(Usplit,p,Qsplit2)
# plt.plot(cx,cy)
# plt.plot(cxs,cys)
# plt.legend(['Initial curve','Segmented curve'])
# plt.plot(Qsplit2[0],Qsplit2[1],'g')
# plt.plot(Qsplit2[0],Qsplit2[1],'r+')
# plt.show()

################ DEGREE ELEVATION EXAMPLE ####################
# P = np.array([[0.0,1.0,0.5,2.5,2.0,3.0],[0.0,0.2,1.0,1.0,0.25,0.0]])
# p = 3
# U = np.array([0,0,0,0,0.4,0.7,1,1,1,1])

# P = np.array([[0.0,0.5,1.0],[0.0,1.0,0.0]])
# p = 2
# U = np.array([0,0,0,1,1,1])

# t = 1

# Qe,Ue,pe = rfn.degreeElevation(U,p,P,t)

# cx,cy = bs.bSplineCurve(U,p,P)
# cxs,cys = bs.bSplineCurve(Ue,pe,Qe)

# plts.plotDegreeElevation(cx,cy,cxs,cys,P,Qe)
# plt.plot(cx,cy)
# plt.plot(cxs,cys)
# plt.legend(['Initial curve','Degree augmented curve'])
# plt.plot(P[0],P[1],'k')
# plt.plot(P[0],P[1],'ro')
# plt.plot(Qe[0],Qe[1],'g')
# plt.plot(Qe[0],Qe[1],'r+')
# plt.show()
