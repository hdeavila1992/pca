import sys
from PyQt5 import QtGui, QtCore
import numpy as np
from numpy import random	


#Ejemplo 2 del libro de optimización


#def readFile(fileName):



def input_inf(pos_cond,Pinip,x,j):
	Long=np.size(pos_cond,1)	
	for i in range(Long): 
		Pinip[pos_cond[0][i],0]=x[j,i]#x[j,(i+1)]  print(x[j,(0+Long)])
		Pinip[pos_cond[0][i],1]=x[j,i+Long]
	return Pinip 	#np.size(pos_cond)

def Extrac_inf(x):
	
	Area_max=np.zeros(np.size(x,0))
	#exec(open("first_optimization.py").read())
	#exec(open("RunIt.py").read())
	for j in range(np.size(x,0)):
		print('Puntos de control iniciales')
		print (Pinp)
		Pinp=input_inf(pos_cond,Pinp,x,j)
		print('Puntos de control aleatorios')
		print (Pinip)
		exec(open("RunIt_temp.py").read())
		Von_max[j]=np.amax(svm_temp)
		Area_max[j]=totalArea_temp
	return Area_max,Von_max

def funcion(x):
	#[Area_max,Von_max]=Extrac_inf(x)
	Von_max=np.zeros(np.size(x,0))
	Area_max=np.zeros(np.size(x,0))
	for j in range(np.size(x,0)):
		print('Puntos de control iniciales')
		print (Pinp)
		temp=Pinp
		temp=input_inf(pos_cond,temp,x,j)
		print('Puntos de control aleatorios')
		print (Pinp)
		exec(open("RunIt_temp.py").read())
		Vect_temp=np.loadtxt('Area_data.dat')
		Area_max[j]=Vect_temp[1]
		print("save successful")		
	#Vn=open('svm_data.txt').read()
	print(Area_max)
	#Von_max[j]=eval (Vn.readline())
	#Area_max[j]=eval (A.readline())
	#print(Von_max)
	#print(Area_max)
	return Area_max

#Constrain function

def g1(x):
	return 26-( (x[:,0]-5)**2) -(x[:,1]**2)	
	
def g2(x):
	return 20-4*x[:,0]-x[:,1]
		
		
def Teaching(x,Fx,ps,Tf,Penal):
#Penalties: Penal
	Minimun=np.amin(abs(Fx(x)))
	pos_Minimun=np.argmin(abs(Fx(x)))
	print('El valor minimo de la funcion es'+str(Minimun))
	print('Se encuentra en '+str(pos_Minimun))
	#=========================
	colums=np.size(x,1)
	Best_l_x=x[pos_Minimun,:]
	#=========================
	mean=(np.ones(ps)@x)/ps
	#========================= 1.1<---------
	Diff_mean_x=random.rand()*(Best_l_x-Tf*mean)
	Base=np.zeros((ps,2))
	Base[:,0]=1.0
	Diff_mean_x_matrix=Base@np.array([Diff_mean_x,Diff_mean_x])
	#=========================1.2<----------
	xp=x+Diff_mean_x_matrix
	print('Posibles valores futuros x_prima=')
	print(xp)
	Fp=Fx(xp)+Penal
	print('Posibles valores futuros de la funcion F(x)')
	print(Fp)
	#=========================
	(x,Fx)=compared(x,xp,Fx(x),Fp)
	return (x,Fx)
	
def checker(P,Q,ps):
	while P==Q:
		P=np.random.randint(0,high=ps)
	return(P,Q)
	
	
def compared(x,xp,F,F_prima):
	pos_cond=np.where(abs(F_prima)<abs(F))
	x[pos_cond,:]=xp[pos_cond,:]+0.0
	F[pos_cond]=F_prima[pos_cond]+0.0
	if any(abs(F_prima)<abs(F))==True:
		print('Se cambia')
		print('F['+str(pos_cond)+']='+str(F[pos_cond])+' por F_prima['+str(pos_cond)+']='+str(F_prima[pos_cond]))
	else:
		print('No se cambia ninguna posicion')
		print(x)
		print(F)	
	return(x,F)

def learning(x,Fx,ps):
	#=============
	P_pos=np.random.randint(0,high=ps)
	Q_pos=np.random.randint(0,high=ps)
	if P_pos==Q_pos:
		(P_pos,Q_pos)=checker(P_pos,Q_pos,ps)
	#=============
	print('Valor de x antes de learning x='+str(x))
	xprima=Learning_fase(Fx(x),x,Q_pos,P_pos)#ojo
	print('Valor de x despues de learning x1='+str(x)+'<no tiene que cambiar')
	print(x==xprima)
	print('Valor de x_prima despues de learning x1_prima='+str(xprima))
	#=============
	(xf,Fxf)=compared(x,xprima,Fx(x),Fx(xprima))
	return (xf,Fxf)
	
def Learning_fase(Fx,x,Q_pos,P_pos):
#Fx Vector de la Función a evaluar
#x Vector de Estudiantes de la población de un tema
#Q_pos Posición aleatoria 1
#Q_pos Posición aleatoria 2
#Ps: Ta maño de la población 
	x_prima=x+0.0
	if abs(Fx[P_pos])<abs(Fx[Q_pos]):
#===============================
		if x[Q_pos,0]<x[P_pos,0]:
			x_prima[Q_pos,:]=x[Q_pos,:] +random.rand()*(x[Q_pos,:]-x[P_pos,:])
		else:
			x_prima[Q_pos,:]=x[Q_pos,:] +random.rand()*(x[P_pos,:]-x[Q_pos,:]) 		
#=====================
	else:
		if x[Q_pos,0]<x[P_pos,0]:
			x_prima[P_pos,:]=x[P_pos,:] +random.rand()*(x[Q_pos,:]-x[P_pos,:])
		else:
			x_prima[P_pos,:]=x[P_pos,:] +random.rand()*(x[P_pos,:]-x[Q_pos,:])
			
	return (x_prima)
	
#Datos iniciales

exec(open("first_optimization.py").read())
exec(open("RunIt.py").read())



rango=0.2	
ps=np.size(Pinp[pos_cond,0],1)
#Extrae la información de los Los puntos de control a estudiar 
print(pos_cond)
print(Pinp[pos_cond,0])
print(Pinp[pos_cond,1])
Pinp_xx=Pinp[pos_cond,0]+0.0
Pinp_yy=Pinp[pos_cond,1]+0.0
Pinp_x=np.append(Pinp_xx,Pinp_yy)+0.0	
	

Base=np.zeros((2*ps,2))
Base[:,0]=1.0

print('La poblacionde puntos x,y iniciales, ORGANIZADO')
print(Pinp_x)

x=np.zeros( (ps,np.size(Pinp_x,0)) )
print(x)

for i in range(ps):
	x[i,:]=Pinp_x+(random.rand(1,2*ps)-0.5)*2.0*rango

print('La poblacion de estudiantes para el algoritmo TLBO es')
print(x)


F=funcion(x)
print('Los valores iniciales de la funcion son:'+str(F))
	
#c1=g1(x)<=0#funcion(x1,x2)
#c2=g2(x)<=0	
	




#print('Valor inicial de F')
#print(F)	


#print('Este es el valor de x1 antes de la fase Teaching<==========')
#print(x1)
#print('Este es el valor de x2 antes de la fase Teaching<===QQQQQQQ')
#print(x2)
	
TerC=abs(F)<=0.01
print('Termination criteria='+str(TerC))	

for i in range(5):
#while any(TerC)==False:
	for j in range(1):
		Penal=0.0#(10*(g1(x)**2)*c1+10*(g2(x)**2)*c2)
		print('Penal por trolazo '+str(Penal))
		(x,F)=Teaching(x,funcion,ps,1,Penal)
		c1=(g1(x)<=0)*1
		c2=(g2(x)<=0)*1
		print('Valor de F despues de '+str(j)+' iteraciones de Teaching')
		print(F) 
	x_l=0.0+x
	#F_l=0.0+F
	for k in range(5):
		(x_l,F_l)=learning(x_l,funcion,ps)	
		print('Valor de F despues de '+str(k)+' iteraciones de Learning')
		print(F_l)
	#print('Este es el valor de x1 Despues de '+str(i)+' iteracion de la fase Learning')
	#print(x1)
	#print('Este es el valor de x2 Despues de '+str(i)+' iteracion de la fase Learning')
	#print(x2)
	(x,F)=compared(x,x_l,funcion(x),funcion(x_l))	
	print('Mejores Valores de F <====================')
	print(F)
	TerC=abs(F)<=0.1
	print(TerC)
	#if any(TerC)==True:
	#	print('Look at this dude !!!')
	#	print(x)
		#print(x2)
	#	sys.exit()
		
	
#coding completed	
	
	
	
