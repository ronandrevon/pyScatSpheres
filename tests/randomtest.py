import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
import random

#-------------------------------------------------------
# Paramètres
#-------------------------------------------------------
ListN=[5,10,20] #5,10,20
nrefp=1.1
h=0.05
EPS=np.arange(0,0.5+h,h)
SEFF=[]
fig,ax = plt.subplots()
for Np in ListN :
	kap=np.ones(Np)
	
	for eps in EPS:
		kdp=eps*np.random.rand(Np)+2.5*kap*np.arange(Np)
		#print(kdp)
		qdot1 = qsa.QdotSphereArray(N=Np,ka=kap,kd=kdp,kp=nrefp)
		
		qdot1.show_ff(ax=ax)



		#-------------------------------------------------------
		# Tracer champ lointain
		#-------------------------------------------------------
		"""
		t = np.linspace(-np.pi,np.pi,100)
		f=qdot1.get_ff(t)
		print(eps)
		F.append(f)
		"""

		#-------------------------------------------------------
		# Calcul de la section efficace
		#-------------------------------------------------------
		seff=qdot1.get_s()
		SEFF.append(seff/(kap**2))

	#-------------------------------------------------------
	# Section efficace
	#-------------------------------------------------------
	print(len(SEFF))
"""
fig, ax = plt.subplots()
ax.plot(EPS, SEFF[0:int(len(SEFF)/3)], 'r',label='N=5',linewidth=3)
ax.plot(EPS, SEFF[int(len(SEFF)/3):int(2*len(SEFF)/3)],'g',label='N=10',linewidth=3)
ax.plot(EPS, SEFF[int(2*len(SEFF)/3):int(len(SEFF))],'b',label='N=20',linewidth=3)
ax.legend(['N=5', 'N=10','N=20'],loc='upper right', frameon=False)
#legend =ax.legend( shadow=True)
#legend.get_frame().set_facecolor('C0')
plt.title('Section efficace en fonction de epsilon')
plt.grid()
plt.show()
"""