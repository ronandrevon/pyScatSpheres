import importlib as imp
import time,pickle
import numpy as np,matplotlib.pyplot as plt
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

# kd_z = np.array([0.0001]+[d*ka]*(N-1))#np.array([0.01,d*ka])
kd_z = np.array([d*ka]*N)#np.array([0.01,d*ka])
kd_y = kd_z
# kd_y = np.array([1]+[0]*(N-1))
# qdot0 = qsa_arb.QdotSphereArray(alpha=0 ,ka=np.array([ka]*N),kd_z=np.array([0,3*ka]),kd_y=np.array([0,0]),kp=kp,nmax=nmax)
# t0=time.time()
qdot2 = qsa_arb.QdotSphereArray(opt2=1,alpha=45,ka=ka,kd_z=np.array([0,1,2,1])*2*ka,kd_y=np.array([0,2,1,1])*2*ka,kp=kp,nmax=nmax)
# qdot2 = qsa_arb.QdotSphereArray(opt2=1,copt=1,alpha=0,ka=np.array([ka]*N),kd_z=kd_z,kd_y=kd_y,kp=kp,nmax=nmax)
# print(time.time()-t0)
# t0=time.time()
# qdot2_0 = qsa_arb.QdotSphereArray(alpha=45,ka=np.array([ka]*N),kd_z=kd_z,kd_y=kd_y,kp=kp,nmax=nmax,opt2=0)
# print(time.time()-t0)

# qdot2 = qsa_arb.QdotSphereArray(alpha=45,ka=np.array([ka]*N),kd_z=kd_z,kd_y=kd_y,kp=kp,nmax=10)
# qdot1 = qsa_arb.QdotSphereArray(alpha=45,ka=ka,kd_z=kd_z[0],kd_y=kd_y[0],kp=kp,nmax=nmax)
#qdot1=qsa.QdotSphereArray(ka=np.array([1,2,1,2]),kd=np.array([2,4,3.5,4]),kp=1.1,nmax=2,copt=1,N=4)
# print(qdot1.ap)
#print("\n")
#print(qdot2.ap)
# print("\n")
# print(qdot1.bp)
#print("\n")
#print(qdot2.bp)


args = {'npts':500,'fz':np.real,'def_args':0,'caxis':[-1,1],'cmap':'jet',
    'r':(0,N*ka*(1+d))*2,'short_title':1,'pOpt':'tXG','fonts':{'title':15}}
args['r']=(0,25)*2
# fig,(ax1,ax2,ax3) = dsp.create_fig('f',rc=[1,3])
# fig,ax = qdot2.show_f(opts='i ',opt='p',**args);
# fig,ax = qdot2.show_f(opts='s ',opt='p',**args);
fig,ax = qdot2.show_f(opts='t ',opt='p',**args);

# ls,ms = np.tile(spu.get_ls_ms(nmax+1),N)
# idx0 = abs(qdot2.ap)>1
# print('l',ls[idx0])
# print('m',ms[idx0])
#
# idx1 = abs(qdot2.L )>1e-3
# ls,ms = np.tile(spu.get_ls_ms(nmax+1),2*N)
# print('l',ls[idx1])
# print('m',ms[idx1])

# for i in range(N):qdot2.show_f(opts='tP',idp=[i],opt='p',**args);
# fig,ax = qdot2_0.show_f(opts='t ',opt='p',**args);

# fig,(ax1,ax2,ax3) = dsp.create_fig('f',rc=[1,3])
# fig,ax3 = dsp.stddisp(opt='')
# qdot2.show_f(opts='iT',ax=ax1,opt='p',**args);
# qdot2.show_f(opts='sP',ax=ax2,idp=[0],xylims=args['r'],opt='p',**args);
# qdot2.show_f(opts='sP',ax=ax3,idp=[1],xylims=args['r'],**args);

# ax3.plot(range(4),range(4),'r-');fig.canvas.draw()
# dsp.stddisp(im=[np.log10(np.maximum(abs(qdot2.T),1e-10))],pOpt='im',cmap='seismic')

# A1 = np.linalg.inv(qdot2.P-qdot2.T)
# fig,(ax1,ax2) = dsp.create_fig('f',rc='12')
# dsp.stddisp(ax=ax1,title='re(T1)',im=[np.real(A1)],pOpt='t',caxis=[-1,1],cmap='seismic',opt='')
# dsp.stddisp(ax=ax2,title='im(T1)',im=[np.imag(A1)],pOpt='t',caxis=[-1,1],cmap='seismic')

# dsp.stddisp(title='lmax=%d' %qdot1.nmax,im=[np.abs(qdot1.T)],pOpt='t',caxis=[0,1],cmap='Blues',opt='p')
# dsp.stddisp(title='lmax=%d' %qdot2.nmax,im=[np.abs(qdot2.T)],pOpt='t',caxis=[0,1],cmap='Blues',opt='p')

# fig,(ax1,ax2) = dsp.create_fig('f',rc='12')
# dsp.stddisp(ax=ax1,title='re(T2)',im=[np.real(qdot2.T)],pOpt='t',caxis=[-10,10],cmap='seismic',opt='')
# dsp.stddisp(ax=ax2,title='im(T2)',im=[np.imag(qdot2.T)],pOpt='t',caxis=[-1,1],cmap='seismic')

# pNmax = qdot1.N*(qdot1.nmax**2)
# Pab = qdot1.P[pNmax+np.arange(pNmax),np.arange(pNmax)]
# Pba = qdot1.P[np.arange(pNmax),pNmax+np.arange(pNmax)]
# print(np.diag(qdot1.P));print(Pab);print(Pba)
# qdot1.show_f(opts='s ',ax=ax21,**args);
# qdot1.show_f(opts='t ',ax=ax31,**args);
