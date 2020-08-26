#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:11:32 2020

@author: hernando
"""
"NURBS"

import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys


def facorial(n):
    factorial=1
    if n>=1:
        for i in range(1,int(n)+1):
            factorial=factorial*i
    return factorial

def cbinomial(n,k):
    cb=facorial(n)/(facorial(k)*facorial(n-k))
    return cb

def B(t,m):
    #m: La cantidad de puntos totales 
    #n: el grado del polinomio n=m-1
    n=m-1
    b=np.zeros(m)
    for i in range(m):
        b[i]=cbinomial(n, i)*(t**i)*((1-t)**(n-i))
    return b


def B2d(u,v,px):
    T=px.shape
    m=T[0]
    n=T[1]
    cx=np.zeros([len(u),len(v)])
    for i in range(len(u)):
        for j in range(len(v)):
            b1=B(u[i],m)
            b2=B(v[j],n)
            cx[i,j]=b1@px@b2
    return cx

def knot (n,p):
    u0=np.array([0]*(p+1))
    if n<p:
        print ("Fatal error")
        sys.exit("no such vector curve exist, pls try n greater than p")
    else:
        j=np.arange(1,n-p+1)
        ujp=j/(n-p+1)
        um=np.array([1]*(p+1))
        k=np.concatenate([u0,ujp,um])
        
    return k

def N0(U,u):
    #U: el vector de nudos
    #u: valor del parametro ente 0 y 1
    N0=((U[0:-1]<=u)*(u<U[1:]))*1
    return N0

def N_1(N0,U,p,u):
    I=0
    N_1=np.array([0.0]*(len(N0)-1))
    for i in range(len(N0[:-1])):
        num1=u-U[i]
        num2=U[i+p+1]-u
        den1=U[i+p]-U[i]
        den2=U[i+p+1]-U[i+1]
        if den1==0:
            A=0
        else:
            A=num1/den1
        if den2==0:
            B=0
        else:
            B=num2/den2
        N_1[I]=A*N0[i]+B*N0[i+1]# 
        I=I+1
        if u==1:
            N_1[-1:]=1
    return N_1

def N(N0,U,p,u):
    for i in range(p):
        N0=N_1(N0,U,(i+1),u)
    return N0    

def ugen(a,b,step):
    if 0<=a and b<=1:
        u=rango=np.arange(a,b,step,dtype=np.float64)
    else:
        print('Los rangos no pertenecen al dominio de la funcion')
        sys.exit("Error message")
    return u  
        
        
def b_spline(u,px,p):
    U=knot(len(px)-1,p)
    cx=np.array([0.0]*len(u))
    for i in range(len(u)):
        Nn0=N0(U,u[i])
        Nn=N(Nn0,U,p,u[i])
        cx[i]=px@Nn
    return cx

def b_spline2d(u,v,px,p,q):
    T=px.shape
    Ui=knot(T[0]-1,p)
    Uj=knot(T[1]-1,q)
    cx=np.zeros([len(u),len(v)])
    for i in range(len(u)):
        for j in range(len(v)):
            Nn0_i=N0(Ui,u[i])
            Nn0_j=N0(Uj,v[j])
            Nn_i=N(Nn0_i,Ui,p,u[i])
            Nn_j=N(Nn0_j,Uj,q,v[j])
            cx[i,j]=Nn_i@px@Nn_j 
    return cx
        
        
def NURBS(u,px,w,p):
    U=knot(len(px)-1,p)
    cx=np.array([0.0]*len(u))
    for i in range(len(u)):
        Nn0=N0(U,u[i])
        Nn=N(Nn0,U,p,u[i])
        cx[i]=(px@(Nn*w))/(Nn@w)
    return cx


def NURBS2d(u,v,px,w,p,q):
    T=px.shape
    Ui=knot(T[0]-1,p)
    Uj=knot(T[1]-1,q)
    cx=np.zeros([len(u),len(v)])
    for i in range(len(u)):
        for j in range(len(v)):
            Nn0_i=N0(Ui,u[i])
            Nn0_j=N0(Uj,v[j])
            Nn_i=N(Nn0_i,Ui,p,u[i])
            Nn_j=N(Nn0_j,Uj,q,v[j])
            cx[i,j]=(Nn_i@(px*w)@Nn_j)/(Nn_i@(w)@Nn_j) 
    return cx

def ploting2d(cx,cy):
    plt.plot(cx,cy)
    return

def ploting3d(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(cx, cy, cz, 'gray')
    if argv != ():
        if argv[0]=='on':
            if len(argv)==4:
                ax.plot3D(argv[1],argv[2],argv[3], 'red')
            elif len(argv)> 4:
                sys.exit("too much argumets, please delet one or more")
            else:
                sys.exit("missing arguments to plot control points")
                
def plotingsurf(cx,cy,cz,*argv):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(cx, cy, cz, 50, cmap='binary')
    if len(argv)==3:
        ax.plot_wireframe(argv[0], argv[1], argv[2], color='red')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');
            
    
    


