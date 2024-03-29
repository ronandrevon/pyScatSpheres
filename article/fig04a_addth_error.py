from utils import*
import scipy.special as spe,time
# from pyScatSpheres.harmo_sphe import harmonique
import pyScatSpheres.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp        ;imp.reload(dsp)
path='../docs/figures/'
plt.close('all')

nmax=10
dp = np.array([0,3,5])
l,m = 4,2

x,y,z,r,theta,phi = spu.spherical_mesh2D(plane='yz',npts=250,r=(-3,10)*2)


xp,yp,zp = dp[:,None]
kd_p,Theta_p,Phi_p = spu.cart2sphere(xp,yp,zp)
a_lmnumu = spu.get_almnumu(l,m,kd_p,Theta_p,Phi_p, Nmax=nmax)
xp,yp,zp = dp
r_p,theta_p,phi_p = spu.cart2sphere(x-xp,y-yp,z-zp)

ls,ms   = spu.get_ls_ms(nmax)
psi_lm  = np.zeros(r.shape,dtype=complex)
print('expansion')
time0 = time.time()

# Ylms = harmonique(nmax,theta_p,'tab',phi_p)
jl  = [spu.jn(nu,r_p) for nu in range(nmax)]
for i,(nu,mu) in enumerate(zip(ls,ms)):
    # print(i,nu,mu)
    # Ylm = Ylms[nu,mu]
    Ylm = spe.sph_harm(mu,nu,phi_p,theta_p)
    psi_lm += a_lmnumu[0,i]*jl[nu]*Ylm
print(time.time()-time0)

# jl  = [spu.jn(nu,r_p) for nu in range(nmax)]
# alnu = spu.a_ln(l,nmax,spu.hn1,kd_p[0],Theta_p[0],Phi_p[0])[-1]
# for i,nu in enumerate(range(nmax)):
#     # print(i,nu,mu)
#     Ylm = spe.sph_harm(0,nu,phi_p,theta_p)
#     psi_lm += alnu[nu]*jl[nu]*Ylm

psi_lm0 = spu.hn1(l,r)*spe.sph_harm(m,l,phi,theta)

ka=1
t = np.linspace(0,2*np.pi,100)
ct,st = np.cos(t),np.sin(t)
plts1 = [[y0+ka*ct,z0+ka*st,'k',''] for y0,z0 in zip([0,yp],[0,zp])]
plts2 = [[y0+ka*ct,z0+ka*st,'k',''] for y0,z0 in zip([0,-yp],[0,-zp])]
# dsp.stddisp(plts1,im=[y,z,psi_lm0.real     ],lw=2,pOpt='im',caxis=[-0.5,0.5],title=r'$\psi_{%d%d}^{(0)}, \nu_{max}=%d$' %(l,m,nmax),cmap='seismic')
# dsp.stddisp(plts2,im=[y-yp,z-zp,psi_lm.real],lw=2,pOpt='im',caxis=[-0.5,0.5],title=r'$\psi_{%d%d}      , \nu_{max}=%d$' %(l,m,nmax),cmap='seismic')

fig,ax=dsp.stddisp(plts2,im=[y-yp,z-zp,np.log10(abs((psi_lm-psi_lm0)/psi_lm0))],
    lw=2,pOpt='im',caxis=[-5,0],labs=['$y_p$','$z_p$'],#title=r'$\log_{10}|\psi_{%d%d}-\psi_{%d%d}^{(0)}|, \nu_{max}=%d$' %(l,m,l,m,nmax))
    fonts={'lab':30,'tick':20},name=path+'addth_error.png',opt='ps')
