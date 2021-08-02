import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
from pyScatSpheres import harmo_sphe as hs
from pyScatSpheres import spherical_utils as spu

import time
import pickle
plt.close('all')

qdot1 = qsa_arb.QdotSphereArray(ka=1,kd_z=np.array([1,3]),kd_y=np.array([0,1]),kp=1.1,nmax=2)
#qdot2=qsa.QdotSphereArray(ka=1,kd=np.array([2,2,2,4]),kp=1.1,nmax=2,copt=1,N=4)
# print(qdot1.ap)
#print("\n")
#print(qdot2.ap)
print("\n")
# print(qdot1.bp)
#print("\n")
#print(qdot2.bp)


fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(ncols=2,nrows=3)
# qdot1.show_f(npts=200,opts='tTG',fz=np.real,r=(0,7,-5,5),caxis=[-1,1]);
args = {'npts':200,'fz':np.real,'def_args':0,'caxis':[-1,1],
    'short_title':1,'pOpt':'t','fonts':{'title':15}}

# r,theta,y,z = spu.polar_mesh2D(1,100,(-3,3,-2,5))
qdot1.show_f(opts='i ',r=(-4,4,-2,7),name='/home/tarik/Downloads/test.png',opt='ps',**args);
# qdot1.show_f(opts='s ',ax=ax21,**args);
# qdot1.show_f(opts='t ',ax=ax31,**args);
