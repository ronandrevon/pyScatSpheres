import importlib as imp
import numpy as np,matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
import random

#-------------------------------------------------------
# Paramètres
#-------------------------------------------------------
ListN=[5,10,20] #5,10,20
kap=1
nrefp=1.1
h=0.05

for Np in ListN :
	EPS=np.arange(0,0.5+h,h)
	SEFF=[]
	fig,ax = plt.subplots()
	for eps in EPS:

		kdp=eps*np.random.rand(Np)+2.5*kap*np.arange(Np)
		print(kdp)

		qdot1 = qsa.QdotSphereArray(N=Np,ka=kap,kd=kdp,kp=nrefp)

		#-------------------------------------------------------
		# Tracer champ lointain
		#-------------------------------------------------------
		F=qdot1.get_ff()
		print(eps)
		

		#-------------------------------------------------------
		# Calcul de la section efficace
		#-------------------------------------------------------
		seff=qdot1.get_s()
		SEFF.append(seff/(kap**2))

	#-------------------------------------------------------
	# Section efficace
	#-------------------------------------------------------
	print(SEFF)
	plt.figure()
	plt.plot(EPS, SEFF, linewidth = 3)
plt.title('Section efficace en fonction de epsilon')
plt.grid()
fig.show()
