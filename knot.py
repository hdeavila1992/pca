import numpy as np

def knot (n,p):
 u0=np.array([0]*(p+1))
 if n==p or n<p:
  print ("no such vector curve exist, pls try n greater than p") 
  u0="NaN"
 else: 
  #range(1,n-p) 
  j=np.array([1]*(n-p-1))
  ujp=j/(n-p)
  um=np.array([1]*(p+1))
  k=np.concatenate([u0,ujp,um])
 return k

#u0=knot(2,1)
#print(u0)
