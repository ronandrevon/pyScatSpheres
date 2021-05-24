#!/usr/bin/env python
# coding: utf-8

# # Spherical quantum dot linear array 
# - [Single sphere](#Single-sphere)
#     - [Simple example](#Simple-example)
# <!-- - [Two spheres](#Two-spheres) -->
# <!-- - [Linear array of spheres](#Linear-array-spheres) -->
# 
# The module for solving the linear array of spherical quantum dot is `qdot_sphere_array` of the `pyScatSpheres` package.

# In[1]:


import importlib as imp
import numpy as np,matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
# import IPython
# %matplotlib notebook


# ## Single sphere
# The limit case of a an array of only 1 sphere can be used as a sanity check for which the analytical solution of the coefficients is known.
# 
# ###  Simple example 
# Chosing a normalised radius $ka=2$, a wave number $k_p=1.5$ 
# 
# Since we use a small normalised radius, only the first few modes appreciably contribute to field expansion. We set `nmax=5`to indicate that only 5 modes will be used to solve the system.

# In[2]:


qdot1 = qsa.QdotSphereArray(N=1,ka=np.pi,kp=1.25,nmax=5)


# The far field pattern can be displayed with `show_ff()`. 
# 
# Since the field is complex, it is possible using characters in the `fopts` argument to display the the real part `r`, the imaginary part `i`, the phase `a`, the magnitude `m`, and the magnitude squared `2`.

# In[3]:


fig,ax = plt.subplots()
qdot1.show_ff(ax=ax,fopts='2');


# The field around the sphere can be obtained with `show_f`.
# 
# By default, the absolute value of the field is displayed 
# 
# It is possible to use the `opts` argument to display the the total field `t`, the scattered field `s` and the incident field `i`. 
# The radial derivative of the field can be displayed by adding the `G` option.

# In[6]:


fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(ncols=2,nrows=3)
# qdot1.show_f(npts=200,opts='tTG',fz=np.real,r=(0,7,-5,5),caxis=[-1,1]); 
args = {'npts':200,'fz':np.real,'def_args':0,'caxis':[-1,1],'short_title':1,'pOpt':'','fonts':{'title':15}}
qdot1.show_f(opts='i ',ax=ax11,**args); 
qdot1.show_f(opts='s ',ax=ax21,**args); 
qdot1.show_f(opts='t ',ax=ax31,**args); 
qdot1.show_f(opts='iG',ax=ax12,**args); 
qdot1.show_f(opts='sG',ax=ax22,**args);
qdot1.show_f(opts='tG',ax=ax32,**args);

