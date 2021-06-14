import numpy as np
import math 

def legendre(l,x,indic):
	P=np.zeros((l+1,2*l+1),dtype=object) 
	P[0,0]=np.ones(x.shape)
	if (l>=1):
		P[1,0],P[1,1]=x,-np.sqrt(1-x**2)
		P[1,-1]=-P[1,1]/2
	if (l>=2):
		for L in range (2,l+1):
			P[L,L]=-(2*L-1)*np.sqrt(1-x**2)*P[L-1,L-1]
			P[L,0]=((2*L-1)/L)*x*P[L-1,0]-(L-1)*P[L-2,0]/L
			P[L,-L]=((-1)**L)*math.factorial(L-L)/math.factorial(L+L)*P[L,L]
			for m in range(1,L):
					P[L,m]=((2*L-1)/(L-m))*x*P[L-1,m]-(L+m-1)*P[L-2,m]/(L-m)
					P[L,-m]=((-1)**m)*math.factorial(L-m)/math.factorial(L+m)*P[L,m]
	if (indic=='tab'):
		return P
	else :
		return P[l,indic]


def harmonique(l,theta,indic):
	Y=legendre(l,np.cos(theta),indic)
	print(Y.shape)
	Y[0,0]=np.sqrt(1/(4*math.pi))*Y[0,0]
	if (l>=1):
		Y[1,0]=np.sqrt((2*1+1)/(4*math.pi))*Y[1,0]
		Y[1,1]=np.sqrt((2*1+1)/(4*math.pi*math.factorial(2*1)))*Y[1,1]
		Y[1,-1]=np.sqrt((2*1+1)*math.factorial(1+1)/(4*math.pi*math.factorial(1-1)))*Y[1,-1]
	if (l>=2):
		for L in range (2,l+1):
			Y[L,L]=np.sqrt((2*L+1)/(4*math.pi*math.factorial(2*L)))*Y[L,L]
			Y[L,0]=np.sqrt((2*L+1)/(4*math.pi))*Y[L,0]
			Y[L,-L]=np.sqrt((2*L+1)*math.factorial(L+L)/(4*math.pi))*Y[L,-L]
			for m in range(1,L):
				pos=np.sqrt((2*L+1)*math.factorial(L-m)/(4*math.pi*math.factorial(L+m)))
				neg=np.sqrt((2*L+1)*math.factorial(L+m)/(4*math.pi*math.factorial(L-m)))
				Y[L,m]=pos*Y[L,m]
				Y[L,-m]=neg*Y[L,-m]
	if (indic=='tab'):
		return Y
	else :
		return Y[l,indic]


