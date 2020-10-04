import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys
import NURBS as bs
import NURBS2 as bs2
import wcg
import zcg

#RESOLVIENDO EL PRIMER EJEMPLO DEL LIBRO DE FEM

#1.STEP: DISCRETIZAR EL MODELO

Px=np.array([0,2.5,5.0,7.5,10.0])
Py=np.array([0.0,0.0,0.0,0.0,0.0])

#Numero de puntos de la cuadratura Gaussiana "npg"
npg=8

print (Px)
print (Py)
#Grado
p=1

U=bs.knot(len(Px)-1,p)

u=U[p:(len(U)-p)]#bs.ugen(0.0,1.0,0.1)
#Llama los "zita" de la cuadratura dependiendo del numero de puntos "npg"
zita=np.array(zcg.Z[npg])
w=np.array(wcg.W[npg])

#Inserta los puntos de Gauss dentro de un intervalo cualquiera A y B
#A: inferior
#B: Superior
def ins_gauss(A,B,zita):
    z1=-1.0
    z2=1.0
    a=(B-A)/(z2-z1)
    b=A-z1*(B-A)/(z2-z1)
    zita_ab=zita*a+b
    return zita_ab
rango=0
px=np.array([ ])
#=np.array.zeros(16,4)
NZ=np.zeros((npg,len(u)))#Tamaño Marticial: [npgxN] 
NZWE=np.zeros((len(u),npg))#Tamaño Matricial: [Nxnpg]
#npg: Numero de puntos de Gauss
#N: Numero de funciones de forma
Gu=np.array([0.0]*(len(u)-1)*npg)
I=0
suma=0
for i in range(len(u)-1):
    Gu[I:(I+npg)]=ins_gauss(u[i],u[i+1],zita)
    u_knot=Gu[I:(I+npg)]
    #[cx_knot,cy_knot]=bs2.bSplineCurveDerivative(U,u_knot,p,Px,Py)
    #suma=suma+(((cx_knot**2+cy_knot**2)**0.5)@w)*(u[i+1]-u[i])*0.5
    #w*np.identity(8)@(u[i+1]-u[i])*0.5*np.identity(8)
    for j in range(len(u_knot)):
        #print (bs2.derNFunction(U,p,u_knot[j])*w[j])
        NZ[j,:]=bs2.derNFunction(U,p,u_knot[j])
        NZWE[:,j]=np.transpose(bs2.derNFunction(U,p,u_knot[j])*w[j])*(u[i+1]-u[i])*0.5
    M=np.matmul(NZWE,NZ)     
        #print (bs2.derNFunction(U,p,u_knot[j]))
    #suma=suma+cx_local@w*(-u[i]+u[i+1])*0.5
    #I=I+npg




cx=bs.b_spline(Gu,Px,p)
cy=bs.b_spline(Gu,Py,p)
cz=bs.b_spline(Gu,Py,p)

#IGABEM

#Propiedades del material
E=1000
#CARGAS
bx=10.0
tx=[0.0]







#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot3D(cx, cy, cz, 'gray')

#plt.show()


