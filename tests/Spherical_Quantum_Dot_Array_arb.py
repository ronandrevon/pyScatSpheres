import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
from pyScatSpheres import harmo_sphe as hs
import time
import pickle

qdot1 = qsa_arb.QdotSphereArray(ka=1,kd_z=np.array([1,3,5]),kd_y=np.array([0,0,0]),kp=1.1,nmax=1)
qdot2=qsa.QdotSphereArray(ka=1,kd=np.array([1,3,5]),kp=1.1,nmax=1,copt=1,N=3)
#print(qdot1.ap)
#print(qdot2.ap)
#print("\n")
print(qdot1.bp)
print(qdot2.bp)

#pour nmax ca marche, mais quand j'augmente le nb d'atomes ca marche plus... grrrr