import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa			;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb	;imp.reload(qsa_arb)
from pyScatSpheres import harmo_sphe as hs					;imp.reload(hs)
from pyScatSpheres import spherical_utils as spu			;imp.reload(spu)
from pyScatSpheres import displayStandards as dsp
from pyScatSpheres import glob_colors as colors

#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1.5,2]),kd_z=np.array([2,6]),kd_y=np.array([0,7]),kp=1.25,nmax=2)
"""
r=np.array([3,3])
theta=np.array([0,1])
fs=qdot1.compute_f(r,theta,0,ftype='t')
"""

#fs=qdot1.compute_f(r_p,theta_p,phi_p,ftype='t')


### Pilotage -----------------------------------------------------

lmax=10
eps=1e-10
N=2
h=5
npts=360
npts2 = int(npts/2)
### Initialisation -----------------------------------------------
# X=np.arange(0,lmax+h,h)
X=np.arange(2,16,2)
Y=np.zeros(X.size)
indice=0

### Boucle -------------------------------------------------------
for l in X:

	print(l)
	# qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1,1]),kd_z=np.array([0,2]),kd_y=np.array([0,0]),kp=1.25,nmax=l)
	qdot1 = qsa_arb.QdotSphereArray(alpha=0,ka=np.array([1,1]),kd_z=np.array([0,6]),kd_y=np.array([0,6]),kp=1.25,nmax=l,opt2=0)
	# qdot1 = qsa_arb.QdotSphereArray(alpha=45,ka=1,kd_z=0,kd_y=0,kp=1.25,nmax=l)
	ka=qdot1.ka

	#### sphère 1
	r1_p=np.array([ka[0]+eps]*npts)
	r1_m=np.array([ka[0]-eps]*npts)
	theta0=np.linspace(0,np.pi,npts2)
	theta1 = np.hstack([theta0,np.flipud(theta0)])
	phi=np.array([np.pi/2]*npts2+[-np.pi/2]*npts2)

	x,y,z = spu.sphere2cart(r1_p,theta1,phi=phi)
	r_p1_p,theta_p1_p,phi_p1_p = spu.cart2sphere(x,y+qdot1.d_py[0],z+qdot1.d_pz[0])

	x,y,z = spu.sphere2cart(r1_m,theta1,phi=phi)
	r_p1_m,theta_p1_m,phi_p1_m = spu.cart2sphere(x,y+qdot1.d_py[0],z+qdot1.d_pz[0])


	#### sphère 2
	r2_p=np.array([ka[1]+eps]*npts)
	r2_m=np.array([ka[1]-eps]*npts)

	x,y,z = spu.sphere2cart(r2_p,theta1,phi=phi)
	r_p2_p,theta_p2_p,phi_p2_p = spu.cart2sphere(x,y+qdot1.d_py[1],z+qdot1.d_pz[1])

	x,y,z = spu.sphere2cart(r2_m,theta1,phi=phi)
	r_p2_m,theta_p2_m,phi_p2_m = spu.cart2sphere(x,y+qdot1.d_py[1],z+qdot1.d_pz[1])

	#### total
	r_tot_p=np.hstack((r_p1_p,r_p2_p))
	theta_tot_p=np.hstack((theta_p1_p,theta_p2_p))
	phi_tot_p=np.hstack((phi_p1_p,phi_p2_p))
	r_tot_m=np.hstack((r_p1_m,r_p2_m))
	theta_tot_m=np.hstack((theta_p1_m,theta_p2_m))
	phi_tot_m=np.hstack((phi_p1_m,phi_p2_m))

	### only sphere 2
	# r_tot_p		= r_p2_p
	# theta_tot_p	= theta_p2_p
	# phi_tot_p	= phi_p2_p
	# r_tot_m		= r_p2_m
	# theta_tot_m	= theta_p2_m
	# phi_tot_m	= phi_p2_m

	### only sphere 1
	# r_tot_p		= r_p1_p
	# theta_tot_p	= theta_p1_p
	# phi_tot_p	= phi_p1_p
	# r_tot_m		= r_p1_m
	# theta_tot_m	= theta_p1_m
	# phi_tot_m	= phi_p1_m

	# calcul différence
	f_plus=qdot1.compute_f(r_tot_p,theta_tot_p,phi_tot_p,ftype='t')
	f_moins=qdot1.compute_f(r_tot_m,theta_tot_m,phi_tot_m,ftype='t')
	f_tot=abs(f_plus-f_moins)
	s=np.sum(f_tot)

	Y[indice]=s
	indice+=1


args = {'npts':300,'fz':np.real,'def_args':0,'caxis':[-1,1],'cmap':'jet',
 'r':(-1.5*ka,4*ka)*2,'short_title':1,'pOpt':'tXG','fonts':{'title':15}}
args['r']=(-2,8,-2,8)
#args['r']=(0.5,3.5)*2
fig,ax = qdot1.show_f(opts='t ',opt='p',**args);

### plot ---------------------------------------------------------------------------
# tle = r'Convergence error of $\Psi_{|r=a_p^+} - \Psi_{|r=a_p^-}$'
dsp.stddisp([X, np.log(Y),'b-o'],lw=2,#title=tle,
	fonts={'lab':30,'tick':20},
	labs=[r'$\nu_{max}$','$log_{10}(|err|)$']);
# plt.figure()
# plt.plot(X, np.log(Y), linewidth = 3)
# plt.grid()
plt.show()
