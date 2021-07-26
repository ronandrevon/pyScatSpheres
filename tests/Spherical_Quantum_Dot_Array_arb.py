import importlib as imp
import numpy as np
import matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
import harmo_sphe as hs
import time
import pickle

qdot1 = qsa_arb.QdotSphereArray_arb(N=3,ka=1,kdz=np.array([1,2,3,4]),kdy=np.array([0,0,0,0]),kp=1.1)

print(qdot1.ap)
print(qdot1.bp)
