#!/usr/bin/env python
# coding: utf-8

# ## 2 quantum dot spheres

# In[1]:


import importlib as imp
import numpy as np,matplotlib.pyplot as plt
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import glob_colors as colors
from pyScatSpheres import displayStandards as dsp
# get_ipython().run_line_magic('matplotlib', 'notebook')


# In[2]:


qdot1 = qsa.QdotSphereArray(N=2,ka=1,kp=1.25,nmax=4,kd=3)


# In[4]:


# fig,((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(ncols=2,nrows=3)
# fig,((ax11,ax12),(ax21,ax22)) = plt.subplots(ncols=2,nrows=2)
# qdot1.show_f(npts=200,opts='tTG',fz=np.real,r=(0,7,-5,5),caxis=[-1,1]);
args = {'npts':400,'fz':np.real,'def_args':0,'caxis':[-1,1],'short_title':1,'pOpt':'','fonts':{'title':15}}
# qdot1.show_f(opts='i ',ax=ax11,**args);
# qdot1.show_f(opts='s ',ax=ax21,**args);
# qdot1.show_f(opts='tG ',**args);
# qdot1.show_f(opts='tG',ax=ax12,**args);
qdot1.show_f(opts='tPG',idp=0,**args);
qdot1.show_f(opts='tPG',idp=1,**args);

plt.show()
# In[ ]:
