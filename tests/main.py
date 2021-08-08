import matplotlib.pyplot as plt
import scipy
from scipy import special
import numpy as np
import harmo_sphe as hs
import math
import time

"""
Paramètres
"""
m=3
l=10
h=0.001
x=np.arange(0,math.pi+h,h)

"""
Calcul des harmoniques sphériques 
"""
#Ma fonction 
n=np.size(x)
y1=np.zeros(n)
y2=np.zeros(n)
i=0
for k in x:
	Y=hs.harmonique(l,k,'tab')
	y1[i]=Y[l,m]
	i+=1
#Python 
y2=scipy.special.sph_harm(m, l, 0, x)

hihi=hs.harmonique(l,x,'tab',0)
#print(hihi)

"""
Plot 
"""
#Mon truc
plt.plot(x, y1, color='red', linewidth = 3)
#Python 
plt.plot(x, y2, color='blue', linestyle='dashed', linewidth = 3) 
#Légendes et affichages
plt.title('Harmoniques sphériques')
plt.grid()
plt.show()


"""
Comparaison du temps
"""
x0=math.pi/2
#Ma fonction
ti1=time.time()
Y=hs.harmonique(l,x,'tab',0)
tf1=time.time()-ti1

#Celle de python 

ti2=time.time()
Ypyt=np.zeros((l+1,2*l+1),dtype=object)
for L in range (0,l+1):
	for m in range(0,L+1):
		Ypyt[L,m]=scipy.special.sph_harm(m, L ,0, x)
		Ypyt[L,-m]=scipy.special.sph_harm(-m, L ,0, x)
tf2=time.time()-ti2

print(tf1)
print(tf2)
#je suis plus rapide hehee efin sauf pour les grands l D:

