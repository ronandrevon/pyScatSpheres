from utils import*
import scipy.special as spe
import pyScatSpheres.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp        ;imp.reload(dsp)
from time import time
path='../../docs/figures/'
opts=''
plt.close('all')
lam=1
k0 =2*np.pi/lam
# lam = 2*np.pi/k0

b = 1
d=3*lam
x,y,z,r,theta,phi = spu.spherical_mesh2D(plane='yz',npts=300,r=(-d,d,-d,d))
# x0 = np.linspace(0.1,2*lam,50)
# y,z = np.meshgrid(np.linspace(-d,d,200),np.linspace(-d,d,200))
x = 0*z
r = np.sqrt(x**2+y**2+z**2)
rb = np.sqrt(x**2+y**2+(z-b)**2)
Z0 = np.exp(1J*k0*rb)/(k0*rb)
Zw = np.zeros(Z0.shape,dtype=complex)
za = abs(z)
if opts=='W'  :
    kk = np.linspace(-k0,k0,200)
    dk = kk[1]-kk[0]
    kxs,kys = np.meshgrid(kk,kk)
    kzs = np.sqrt(k0**2-kxs**2-kys**2+0J)
    kxs,kys,kzs = kxs.flatten(),kys.flatten(),kzs.flatten()
    idx = (abs(kzs.imag)==0) & (abs(kzs.real)>0) # >0 | (kzs.real>0)
    kxs = kxs[idx]
    kys = kys[idx]
    kzs = kzs[idx]
    for i,(kx,ky,kz) in enumerate(zip(kxs,kys,kzs)):
        Zw += np.exp(1J*(kx*x+ky*y+kz*za))/kz
    Zw *= dk**2/(2J*np.pi)
elif opts=='S':
    rho = np.sqrt(x**2+y**2)
    dkp = 0.1
    kps = np.arange(0,10)*dkp
    kzs = np.sqrt(k0**2-np.array(kps,dtype=complex)**2)
    for kp,kz in zip(kps,kzs):
        Zw += kp/kz*spe.jn(0,kp*rho)*np.exp(-1J*kz*z)
    Zw *=1J**dkp
else:
    lmax = 20
    # Zw=1J*k0*spu.hn1(0,k0*r)
    rplus  = np.maximum(r,b)
    rminus = np.minimum(r,b)
    for l in range(lmax)  :
        yl = spe.sph_harm(0,l,phi,theta)*np.conj(spe.sph_harm(0,l,np.pi/2,0))
        # Zw+=spu.hn1(l,k0*r)*yl
        Zw+=spu.hn1(l,k0*rplus)*spu.jn(l,k0*rminus)*yl
        # Zw+=spu.jn(l,rminus)*yl
    Zw *=1J*k0


dsp.stddisp(im=[y,z,np.imag(Z0)],title='Z0',pOpt='im',caxis=[-1,1],cmap='seismic')
dsp.stddisp(im=[y,z,np.imag(Zw)],title='Zw',pOpt='im',caxis=[-1,1],cmap='seismic')
dsp.stddisp(im=[y,z,np.log10(np.abs((Zw-Z0)/Z0))],title='Zw',pOpt='im',caxis=[-5,0],cmap='jet')
