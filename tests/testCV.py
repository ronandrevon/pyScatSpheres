import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
from pyScatSpheres import harmo_sphe as hs
from pyScatSpheres import spherical_utils as spu

#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1.5,2]),kd_z=np.array([2,6]),kd_y=np.array([0,7]),kp=1.25,nmax=2)
"""
r=np.array([3,3])
theta=np.array([0,1])
fs=qdot1.compute_f(r,theta,0,ftype='t')
"""

#fs=qdot1.compute_f(r_p,theta_p,phi_p,ftype='t')


### Pilotage -----------------------------------------------------
lmax=40
eps=0.0001
N=2
h=5

### Initialisation -----------------------------------------------
X=np.arange(0,lmax+h,h)
Y=np.zeros(X.size)
indice=0

### Boucle -------------------------------------------------------
for l in X:
	qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1.5,2]),kd_z=np.array([2,6]),kd_y=np.array([0,7]),kp=1.25,nmax=l)
	ka=qdot1.ka
	# sphère 1
	r1_p=np.array([ka[0]+eps]*360)
	r1_m=np.array([ka[0]-eps]*360)
	theta1=np.arange(0,2*np.pi,2*np.pi/360)

	x,y,z = spu.sphere2cart(r1_p,theta1,phi=np.pi/2)
	r_p1_p,theta_p1_p,phi_p1_p = spu.cart2sphere(x,y+qdot1.d_py[0],z+qdot1.d_pz[0])

	x,y,z = spu.sphere2cart(r1_m,theta1,phi=np.pi/2)
	r_p1_m,theta_p1_m,phi_p1_m = spu.cart2sphere(x,y+qdot1.d_py[0],z+qdot1.d_pz[0])


	# sphère 2
	r2_p=np.array([ka[1]+eps]*360)
	r2_m=np.array([ka[1]-eps]*360)

	x,y,z = spu.sphere2cart(r2_p,theta1,phi=np.pi/2)
	r_p2_p,theta_p2_p,phi_p2_p = spu.cart2sphere(x,y+qdot1.d_py[1],z+qdot1.d_pz[1])

	x,y,z = spu.sphere2cart(r2_m,theta1,phi=np.pi/2)
	r_p2_m,theta_p2_m,phi_p2_m = spu.cart2sphere(x,y+qdot1.d_py[1],z+qdot1.d_pz[1])

	#total
	r_tot_p=np.hstack((r_p1_p,r_p2_p))
	theta_tot_p=np.hstack((theta_p1_p,theta_p2_p))
	phi_tot_p=np.hstack((phi_p1_p,phi_p2_p))

	r_tot_m=np.hstack((r_p1_m,r_p2_m))
	theta_tot_m=np.hstack((theta_p1_m,theta_p2_m))
	phi_tot_m=np.hstack((phi_p1_m,phi_p2_m))

	# calcul différence
	f_plus=qdot1.compute_f(r_tot_p,theta_tot_p,phi_tot_p,ftype='t')
	f_moins=qdot1.compute_f(r_tot_m,theta_tot_m,phi_tot_m,ftype='t')
	f_tot=abs(f_plus-f_moins)
	s=np.sum(f_tot)

	Y[indice]=s
	indice+=1

### plot ---------------------------------------------------------------------------
plt.figure()
plt.plot(X, np.log(Y), linewidth = 3)
plt.grid()
plt.show()