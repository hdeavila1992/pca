import numpy as np

def F_BASE(U,u,p):
 vect=range(0,len(U)-1)
 N0=np.array(vect)*0.0

 for i in vect:  
  N0[i]=(U[i]<=u)*(u<U[i+1])
 I=len(N0)-1
 Np_1=np.array(range(0,I))*0.0
 for pn in range(1,p+1):
  for i in range(0,I):
   #range(0,I):
   den1=U[i+pn]-U[i]
   den2=U[i+pn+1]-U[i+1]
   #Primera condicion
   if den1 != 0:
    A=(u-U[i])/den1
   else:
    A=0.0
   #Segunda condicion
   if den2 != 0:
    B=(U[i+pn+1]-u)/den2
   else:
    B=0.0
   Np_1[i]=A*N0[i]+B*N0[i+1]
  I=len(Np_1)-1
  #print(I)
  N0=Np_1
  #print(N0)
  Np_1=np.array(range(0,I))*0.0
 if u==1:
   #N0[-1:]=1
   N0=1  
 return N0
    
 
   
     

#        %Primera condicion
#        if den1~=0
#A=(u-U(i))/den1;
#        else
#            A=0;
#        end
#        %Segunda condicion
#        if den2~=0
#B=(U(i+pn+1)-u)/den2;
#        else
#            B=0;
#        end
#%================================


def knot(n,p):
 #Esta funcion Genera el vector de nudos
 u0=np.zeros(p+1)
 j=range(1,(n))
 ujp=np.array(j)/(n-p+1.0)
 um=np.ones(p+1)
 u0=np.append(u0,ujp)
 k=np.append(u0,um)
 return k




#function [ k ] = knot(n,p)
#%%Generador de nodos 
#u0=0*ones(1,p+1);
#j=1:(n-p);
#ujp=j/(n-p+1);
#um=ones(1,p+1);
#k=[u0,ujp,um];
#end



