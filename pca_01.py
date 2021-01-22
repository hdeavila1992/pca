#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 15:11:32 2020
Modified on Sun Aug 02 15:29:17 2020

@author: hernando
@collaborator: paternina
"""
#"NURBS"

import numpy as np
import sys

#tolerance for knot vector
tol = 1e-5

def knotGeneratorUniform(n,p):
    if n<p:
        print ("Fatal error")
        sys.exit("No such vector curve exists, please try n greater than p")
        return np.zeros(p+1)
    else:
        ustart = np.zeros(p+1)
        umid = np.zeros(n-p)
        for j in range(1,n-p+1):
            umid[j-1] = j/(n-p+1)
        uend = np.ones(p+1)
        return np.concatenate([ustart,umid,uend])

def knotGeneratorChord(n,p,tv):
    if n<p:
        print ("Fatal error")
        sys.exit("No such vector curve exists, please try n greater than p")
        return np.zeros(p+1)
    else:
        ustart = np.zeros(p+1)
        uend = np.ones(p+1)
        umid = np.zeros(n-p)
        for j in range(1,n-p+1):
            ui = 0
            for i in range(j,j+p):
                ui += tv[i]
            ui *= (1.0/p)
            umid[j-1] = ui
        return np.concatenate([ustart,umid,uend])

#U: knot vector
#u: parameter between 0 and 1
def findKnotInterval(U,u):
    for i in range(len(U)-1):
        if U[i] < (u + tol) and (u + tol) < (U[i+1]):
            # nbasis[i] = 1.0
            index = i

        if abs(u - U.max())<tol:
            if U[i] < (u - tol) and (u - tol) < U[i+1]:
                # nbasis[i] = 1.0
                index = i
    return index

####################################################
#################BASIS FUNCTIONS####################
####################################################

def nFunction(U,p,u):
    m = len(U) - 1
    for pi in range(p+1):
        # print("p-Order: ",pi)
        if pi != 0:
            nbas = np.zeros((1,len(nbasis[0])-1))
            for j in range(m-pi):

                num1 = u - U[j]
                den1 = U[j+pi] - U[j]

                num2 = U[j+pi+1] - u
                den2 = U[j+pi+1] - U[j+1]

                if abs(den1)<tol:
                    A = 0.0
                else:
                    A = num1/den1

                if abs(den2)<tol:
                    B = 0.0
                else:
                    B = num2/den2

                nbas[0][j] = A*nbasis[0][j] + B*nbasis[0][j+1]

            nbasis = nbas
        else:
            nbasis = np.zeros((1,len(U)-1))
            for i in range(len(nbasis[0])):
                if U[i] < (u + tol) and (u + tol) < (U[i+1]):
                    nbasis[0][i] = 1.0

                if abs(u - U.max())<tol:
                    if U[i] < (u - tol) and (u - tol) < U[i+1]:
                        nbasis[0][i] = 1.0
    return nbasis

def derNFunction(U,p,u):
    m = len(U) - 1
    dnbasis = nFunction(U,p-1,u)
    dnbas = np.zeros((1,len(dnbasis[0])-1))
    for j in range(m-p):

        num1 = 1.0
        den1 = U[j+p] - U[j]

        num2 = 1.0
        den2 = U[j+p+1] - U[j+1]

        if abs(den1)<tol:
            A = 0.0
        else:
            A = num1/den1

        if abs(den2)<tol:
            B = 0.0
        else:
            B = num2/den2

        dnbas[0][j] = p*(A*dnbasis[0][j] - B*dnbasis[0][j+1])

    dnbasis = dnbas
    return dnbasis

def secDerNFunction(U,p,u):
    m = len(U) - 1
    d2nbasis = derNFunction(U,p-1,u)
    d2nbas = np.zeros((1,len(d2nbasis[0])-1))
    for j in range(m-p):

        num1 = 1.0
        den1 = U[j+p] - U[j]

        num2 = 1.0
        den2 = U[j+p+1] - U[j+1]

        if abs(den1)<tol:
            A = 0.0
        else:
            A = num1/den1

        if abs(den2)<tol:
            B = 0.0
        else:
            B = num2/den2

        d2nbas[0][j] = p*(A*d2nbasis[0][j] - B*d2nbasis[0][j+1])

    d2nbasis = d2nbas
    return d2nbasis
