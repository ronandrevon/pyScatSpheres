import importlib as imp
import time,pickle
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import glob_colors as colors
from utils import displayStandards as dsp;imp.reload(dsp)
# from pyScatSpheres import displayStandards as dsp
from pyScatSpheres import harmo_sphe as hs
from pyScatSpheres import spherical_utils as spu
plt.close('all')
N,nmax,ka,kp = 4,6, 2,1.25#1

d=3
# d=3/2*np.sqrt(2)
# kd_z,kd_y = [0,d*ka]
# qdot1 = qsa_arb.QdotSphereArray(ka=np.array(ka),kd_z=np.array(kd_z),kd_y=np.array(kd_y),kp=kp,nmax=nmax)

#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1.5,2]),kd_z=np.array([2,6]),kd_y=np.array([0,7]),kp=1.25,nmax=10)
qdot1 = qsa_arb.QdotSphereArray(ka=np.array([1,2,1,2,1]),kd_z=np.array([2,4,3.5,4,3]),kd_y=np.array([0,5,-2,1,5]),kp=1.25,nmax=7,alpha=0,opt2=0)
#qdot1 = qsa.QdotSphereArray(ka=np.array([1,2,1,2,1]),kd=np.array([2,4,3.5,4,3]),kp=1.25,nmax=10,N=5)

#qdot1 = qsa_arb.QdotSphereArray(ka=np.array([2.4,4,2.4]),kd_z=np.array([5,-4.5,4.5]),kd_y=np.array([0,4.5,4.5]),kp=1.25,nmax=7,opt2=0)


#qdot1=qsa.QdotSphereArray(ka=np.array([1,2,1,2]),kd=np.array([2,4,3.5,4]),kp=1.1,nmax=2,copt=1,N=4)
# print(qdot1.ap)
#print("\n")
#print(qdot2.ap)
# print("\n")
# print(qdot1.bp)
#print("\n")
#print(qdot2.bp)


#fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(ncols=2,nrows=3)
# qdot1.show_f(npts=200,opts='tTG',fz=np.real,r=(0,7,-5,5),caxis=[-1,1]);
args = {'npts':500,'fz':np.real,'def_args':0,'caxis':[-1,1],
    'short_title':1,'pOpt':'t','fonts':{'title':15}}

# r,theta,y,z = spu.polar_mesh2D(1,100,(-3,3,-2,5))
#qdot1.show_f(opts='i ',r=(-4,4,-2,7),name='/home/tarik/Downloads/test.png',opt='p',**args);
ti1=time.time()
#qdot1.show_f(opts='t ',r=(-4,15,-5,15),name='/home/tarik/Downloads/test.png',opt='p',**args);
tf1=time.time()-ti1

print(tf1)
# pNmax = qdot1.N*(qdot1.nmax**2)
# Pab = qdot1.P[pNmax+np.arange(pNmax),np.arange(pNmax)]
# Pba = qdot1.P[np.arange(pNmax),pNmax+np.arange(pNmax)]
# print(np.diag(qdot1.P));print(Pab);print(Pba)
# qdot1.show_f(opts='s ',ax=ax21,**args);
# qdot1.show_f(opts='t ',ax=ax31,**args);
qdot1.show_ff()