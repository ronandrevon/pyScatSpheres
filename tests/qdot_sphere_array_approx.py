import importlib as imp
from utils import *
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)

s1 = qsa.QdotSphereArray(N=1,ka=10,kd=0 ,kp=1.1,nmax=15,solve=True);s1.show_ff()
s2 = qsa.QdotSphereArray(N=2,ka=10,kd=22,kp=1.1,nmax=15,solve=True);s2.show_ff()
s2 = qsa.QdotSphereArray(N=2,ka=10,kd=25,kp=1.1,nmax=15,solve=True);s2.show_ff()
s2 = qsa.QdotSphereArray(N=2,ka=10,kd=30,kp=1.1,nmax=15,solve=True);s2.show_ff()
