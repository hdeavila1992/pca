import numpy as np
import matplotlib.pyplot as plt

factorialList = dict([
(0,1),
(1,1),
(2,2),
(3,6),
(4,24),
(5,120),
(6,720),
(7,5040),
(8,40320),
(9,362880),
(10,3628800)
])

def bezierBasisFunction(i,n,u):
    # bn = (np.math.factorial(n)/(np.math.factorial(i)*np.math.factorial(n-i)))*(u**i)*((1.0 - u)**(n-i))
    bn = (factorialList[n]/(factorialList[i]*factorialList[n-i]))*(u**i)*((1.0 - u)**(n-i))
    return bn

def bezierCurve(ubz,P,n):
    cx = np.zeros(len(ubz))
    cy = np.zeros(len(ubz))

    for j in range(len(ubz)):
        for i in range(len(P[0].T)):
            bn = bezierBasisFunction(i,n,ubz[j])
            cx[j] += bn*P[0][i]
            cy[j] += bn*P[1][i]

    return cx,cy

def plotBezierCurve(cx,cy,P):
    fig = plt.figure()
    plt.plot(cx,cy)
    plt.plot(P[0],P[1],'ro')
    plt.plot(P[0],P[1],'k')
    plt.show()

# P = np.array([[0,1],[0,1]])
# n = 1

P = np.array([[0.0,0.5,1.0],[0.0,1.0,0.0]])
n = 2

numpoints = 21
Ubz = np.linspace(0.0,1.0,numpoints)

cx,cy = bezierCurve(Ubz,P,n)
plotBezierCurve(cx,cy,P)
