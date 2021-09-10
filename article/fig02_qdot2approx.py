from utils import *
import matplotlib.patches as patches
import scipy,copy,time
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import spherical_utils as spu   ;imp.reload(spu)

name='qdot2_approx_img.png'
opt='ps'

nmax,ka,k0=5,6,1
kd  = 4.5*ka
nka = 5
nka2 = 2.5*nka
npts = 200
fig,ax = dsp.stddisp(figsize='f',opt='')

### plane waves
y0,z0 = np.meshgrid(np.linspace(-1.5*ka,1.5*ka,npts),np.linspace(-1.25*ka,0,npts))
z0-=1.5*ka
Z0 = np.exp(1J*k0*z0).real#*np.exp(-0.01*y**2)
dsp.stddisp(ax=ax,im=[y0,z0,Z0],lw=2,pOpt='Xe',caxis=[-1,1],opt='')
dsp.stddisp(ax=ax,im=[y0+nka2*ka,z0,Z0  ],lw=2,pOpt='Xe',caxis=[-1,1],opt='')

### spherical waves
r0 = np.linspace(1.1*ka,2.5*ka,100)
theta_d = np.linspace(-60,60,100)
theta0  = np.deg2rad(theta_d)
r,theta = np.meshgrid(r0,theta0)
phi = np.pi/2
qdot1 = qsa.QdotSphereArray(N=1,ka=ka,kp=kp,nmax=nmax,solve=0);qdot1.solve(opt2=1)
f = qdot1.compute_f(r,theta,phi,ftype='s')
x,y,z = spu.sphere2cart(r,theta,phi)
Z=f.real/f.real.max()
dsp.stddisp(ax=ax,im=[y,z+kd,Z],lw=2,pOpt='Xe',caxis=[-1,1],opt='')
dsp.stddisp(ax=ax,im=[y+nka*ka,z,Z],lw=2,pOpt='Xe',caxis=[-1,1],opt='')

### solution to spherical wave
vl = qdot1.vl
bp2a = vl*np.sum(qdot1.bp*spu.a_ln(nmax,nmax,spu.hn1,kd,0,0).T,axis=1)
qdot1.bp = bp2a
f1a = qdot1.compute_f(r,theta,phi,ftype='s')
Z1a=f1a.real/f1a.real.max()
dsp.stddisp(ax=ax,im=[y+nka*ka,z+kd,Z1a],lw=2,pOpt='Xe',caxis=[-1,1],opt='')

#### last sphere
qdot2 = qsa.QdotSphereArray(N=2,ka=ka,kd=kd,kp=kp,nmax=nmax,solve=1)
# theta_d2 = np.linspace(-30,30,100)
# theta02 = np.deg2rad(theta_d2)
r2,theta2,phi2 = spu.cart2sphere(0,y,z+kd)
f2 = qdot2.compute_f(r2,theta2,phi2,ftype='s')
Z2 = f2.real/f.real.max()
dsp.stddisp(ax=ax,im=[y+nka2*ka,z+kd,Z2],lw=2,pOpt='Xe',caxis=[-1,1],opt='')


#### patches arrows, and texts
kad=-2.5*ka
pp=[]
pp = [patches.FancyArrowPatch((kad, 0), (kad, kd), arrowstyle='<|-|>',
    mutation_scale=20,color='k',alpha=0.8,lw=2)]
pp+=[patches.Circle((kap,kdp),ka,color='b',alpha=0.25,ec='k',lw=2)
    for kap,kdp in zip([0,nka*ka,nka2*ka,nka2*ka],[kd,kd,0,kd])]
txt = [[1.2*kad,kd/2,'$kd$','k']]
txt+= [[0.75*nka/2*ka,kd/2,'$+$','k']]
n2,nl = (nka2+nka)/2,2
txt+= [[1.05*n2*ka,1.4*kd/2,'?','k']]
plt.arrow((n2-0.4)*ka,kd/2,0.8*nl*ka,0,width=1)
dsp.stddisp(texts=txt,fonts={'text':30},patches=pp,lw=2,pOpt='Xe',ax=ax,
    opt=opt,name=name)
