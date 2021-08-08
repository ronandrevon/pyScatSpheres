import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
import random
from pyScatSpheres import spherical_utils as spu


plt.close('all')
"""
indic=1 : Validate by checking continuity of the wave function and radial derivatives at interfaces
indic=2 : optional convergence test using HardSphereArrayBase.test_convergence()
- Example case N=5,10,20; ka=1; kd=2.5ka+epsnp.random.rand(N); nref=1.1 where esp varies in 0-0.5.
indic=4 : Show far field with show_ff for each values of eps
indic=3 : Show sigma(eps) for N=5,10 et 20 using HardSphereArrayBase.get_s()
"""
### Paramètres de controle ----------------------------------------------------
indic=4


### Continuity check -----------------------------------------------------------
if (indic==1):
	qdot1 = qsa.QdotSphereArray(N=3,ka=np.array([np.pi]*3),kp=1.25,nmax=5,kd=np.array([0,7,10]))
	fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(ncols=2,nrows=3)
	args = {'npts':200,'fz':np.real,'def_args':0,'caxis':[-1,1],'short_title':1,'pOpt':'','fonts':{'title':15}}
	qdot1.show_f(opts='i ',ax=ax11,**args,r=(-4,5,-5,35)); 
	qdot1.show_f(opts='s ',ax=ax21,**args,r=(-4,5,-5,35)); 
	qdot1.show_f(opts='t ',ax=ax31,**args,r=(-4,5,-5,35)); 
	qdot1.show_f(opts='iG',ax=ax12,**args,r=(-4,5,-5,35)); 
	qdot1.show_f(opts='sG',ax=ax22,**args,r=(-4,5,-5,35));
	qdot1.show_f(opts='tG',ax=ax32,**args,r=(-4,5,-5,35));

	lmin=1
	lmax=30
	eps=0.000000001
	h=1

	X=np.arange(lmin,lmax+h,h)
	Y=np.zeros(X.size)
	indice=0

	for l in X:
		#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1,1]),kd_z=np.array([2,6]),kd_y=np.array([0,7]),kp=1.25,nmax=l)
		#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1,1,2,2]),kd_z=np.array([2,4,3,5]),kd_y=np.array([0,5,-2,1]),kp=1.25,nmax=l)

		qdot1=qsa.QdotSphereArray(N=2,ka=np.array([np.pi]*2),kp=1.25,nmax=l,kd=np.array([0,7]))
		ka=qdot1.ka
		# sphère 1
		r1_p=np.array([ka[0]+eps]*360)
		r1_m=np.array([ka[0]-eps]*360)
		theta1=np.arange(0,2*np.pi,2*np.pi/360)

		x,y,z = spu.sphere2cart(r1_p,theta1,phi=np.pi/2)
		r_p1_p,theta_p1_p,phi_p1_p = spu.cart2sphere(x,y,z+qdot1.d_p[0])

		x,y,z = spu.sphere2cart(r1_m,theta1,phi=np.pi/2)
		r_p1_m,theta_p1_m,phi_p1_m = spu.cart2sphere(x,y,z+qdot1.d_p[0])


		# sphère 2
		r2_p=np.array([ka[1]+eps]*360)
		r2_m=np.array([ka[1]-eps]*360)

		x,y,z = spu.sphere2cart(r2_p,theta1,phi=np.pi/2)
		r_p2_p,theta_p2_p,phi_p2_p = spu.cart2sphere(x,y,z+qdot1.d_p[1])

		x,y,z = spu.sphere2cart(r2_m,theta1,phi=np.pi/2)
		r_p2_m,theta_p2_m,phi_p2_m = spu.cart2sphere(x,y,z+qdot1.d_p[1])

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
	print(Y)
	plt.figure()
	plt.plot(X, np.log(Y), linewidth = 3,color='red')
	plt.grid()
	plt.title('log de la valeur absolue de lerreur en fonction de lmax')
	plt.show()





### Convergence test ----------------------------------------------------------
if (indic==2):
	qdot1 = qsa.QdotSphereArray(N=3,ka=np.array([np.pi]*3),kp=1.25,nmax=5,kd=np.array([0,7,10]))
	qdot1.test_convergence(nmaxs=np.array([1,2,3,4,5,6,7]))
"""
### Section efficace ---------------------------------------------------------------
Example case N=5,10,20; ka=1; kd=2.5ka+epsnp.random.rand(N); nref=1.1 where esp varies in 0-0.5.
- Show sigma(eps) for N=5,10 et 20 using HardSphereArrayBase.get_s()
"""
if (indic==3):
	#-------------------------------------------------------
	# Paramètres
	#-------------------------------------------------------
	ListN=[5,10,20] #5,10,20
	nrefp=1.1
	h=0.05
	EPS=np.arange(0,0.5+h,h)
	SEFF=[]
	for Np in ListN :
		kap=np.ones(Np)

		for eps in EPS:
			kdp=eps*np.random.rand(Np)+2.5*kap*np.arange(Np)
			#print(kdp)
			qdot1 = qsa.QdotSphereArray(N=Np,ka=kap,kd=kdp,kp=nrefp)
			#-------------------------------------------------------
			# Calcul de la section efficace
			#-------------------------------------------------------
			seff=qdot1.get_s()
			SEFF.append(seff/(kap**2))

		#-------------------------------------------------------
		# Section efficace
		#-------------------------------------------------------	
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

### Far field ---------------------------------------------------------------
Example case N=5,10,20; ka=1; kd=2.5ka+epsnp.random.rand(N); nref=1.1 where esp varies in 0-0.5.
- Show far field with show_ff for each values of eps
"""
if (indic==4):
    N=20
    nrefp=1.1
    h=0.05
    EPS=np.arange(0,0.5+h,h)
    kap=np.ones(N)
    fig, ax = plt.subplots()
    for eps in EPS:
        kdp=eps*np.random.rand(N)+2.5*kap*np.arange(N)
        #print(kdp)
        qdot1 = qsa.QdotSphereArray(N=N,ka=kap,kd=kdp,kp=nrefp)

        #theta_d = np.rad2deg(theta)
        theta  = np.linspace(0,np.pi,100)
        ff = qdot1.get_ff(theta)
        #plts,ls = [], ['--','-'][leg]
        #plts+= [[theta_d,np.abs(ff)     ,['b'          ,ls],['',r'$||$'  ][leg] ]]
        #title='Scattering amplitude for %s' %qdot1._params()
        ax.plot(theta,ff,label='eps=%.2f'%(eps))
        legend =ax.legend(loc='upper right', shadow=True,)
    plt.title('Far field N=20')
    plt.grid()
    plt.show()