from utils import*
import scipy.special as spe
import pyScatSpheres.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp        ;imp.reload(dsp)
path='../../docs/figures/'
plt.close('all')

nmax=10
dp = np.array([0,2,2])
l,m = 5,-4

x,y,z,r,theta,phi = spu.spherical_mesh2D(plane='yz',npts=500,r=(0,5)*3)


xp,yp,zp = dp[:,None]
kd_p,Theta_p,Phi_p = spu.cart2sphere(xp,yp,zp)
a_lmnumu = spu.get_almnumu(l,m,kd_p,Theta_p,Phi_p, Nmax=nmax)
xp,yp,zp = dp
r_p,theta_p,phi_p = spu.cart2sphere(x-xp,y-yp,z-zp)

ls,ms   = spu.get_ls_ms(nmax)
psi_lm  = np.zeros(r.shape,dtype=complex)
print('expansion')
for i,(nu,mu) in enumerate(zip(ls,ms)):
    print(i)
    Ylm = spe.sph_harm(mu,nu,phi_p,theta_p)
    jl = spu.jn(nu,r_p)
    psi_lm += a_lmnumu[0,i]*jl*Ylm

psi_lm0 = spu.hn1(l,r)*spe.sph_harm(m,l,phi,theta)

ka=1
t = np.linspace(0,2*np.pi,100)
ct,st = np.cos(t),np.sin(t)
plts1 = [[y0+ka*ct,z0+ka*st,'k',''] for y0,z0 in zip([0,yp],[0,zp])]
plts2 = [[y0+ka*ct,z0+ka*st,'k',''] for y0,z0 in zip([0,-yp],[0,-zp])]
# dsp.stddisp(plts1,im=[y,z,abs(psi_lm0)      ],lw=2,pOpt='im',caxis=[0,0.5],title=r'$\psi_{%d%d}^{(0)}, \nu_{max}=%d$' %(l,m,nmax))
# dsp.stddisp(plts2,im=[y-yp,z-zp,abs(psi_lm )],lw=2,pOpt='im',caxis=[0,0.5],title=r'$\psi_{%d%d}      , \nu_{max}=%d$' %(l,m,nmax))

dsp.stddisp(plts2,im=[y-yp,z-zp,np.log10(abs(psi_lm-psi_lm0))],lw=2,pOpt='im',caxis=[-5,0],title=r'$\log_{10}|\psi_{%d%d}-\psi_{%d%d}^{(0)}|, \nu_{max}=%d$' %(l,m,l,m,nmax))
