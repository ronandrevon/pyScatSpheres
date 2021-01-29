import importlib as imp
import time
from utils import*
import scipy.special as spe
import optics.scattering.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp            ;imp.reload(dsp)
path='../../docs/figures/'
plt.close('all')

import sage.all
import py3nj

def plot_psi_lm(l,m,zl='in',zq='',bbt={}, phi=0,cart_opt=True,npts=100,r=(1e-3,5),fz=np.abs,name='',caxis=None,**kwargs):
    ''' plot psi_lm
    l,m : indices of psi
    zl : str - 'in'(jn), 'out'(hn1)
    '''
    r,theta,y,z = spu.polar_mesh2D(cart_opt,npts,r)
    z_l = {'in':spu.jn,'out':spu.hn1}[zl]
    tle = r'$\psi_{%d,%d}^{(%s)}$' %(l,m,zl)

    zlzq = zl
    if zq :
        zlzq += '-'+zq
        z_q = {'in-in':spu.jn,'out-out':spu.jn ,'out-in':spu.hn1,'in-out':spu.hn1}[zlzq]
        z_l = {'in-in':spu.jn,'out-out':spu.hn1,'out-in':spu.jn ,'in-out':spu.hn1}[zlzq]
        r_d,theta_d,phi_d,Nmax = [bbt[k] for k in ['r_d','theta_d','phi_d','Nmax']]
        # xd,yd,zd        = spu.sphere2cart(r_d,theta_d,phi_d)
        # rp,thetap,phip  = spu.cart2sphere(-xd,y-yd,z-zd)
        psi_lm = spu.translation_addth_scalar(r,theta,phi,l,m,zl=z_l,zq=z_q,**bbt)
        tle += r', %s, $r_d=%.1f$, $\theta_d=%.1f^{\circ}$, $\phi_d=%.1f^{\circ}$,$N_{max}$=%d'%(zlzq,r_d,np.rad2deg(theta_d),np.rad2deg(phi_d),Nmax)
    else:
        psi_lm = z_l(l,r)*spe.sph_harm(m,l,phi,theta)

    nameP = '%d%d_%s.png' %(l,m,zlzq)
    if not caxis : caxis= [0,min(2,fz(psi_lm).max())]
    dsp.stddisp(im=[y,z,fz(psi_lm)],title=tle,labs=['$y$','$z$'],caxis=caxis,name=name+nameP,**kwargs)
    return y,z,psi_lm





###########################################################################################################
#                                               Tests
###########################################################################################################
opts= '' #i(in-in), o(out-out)g(gaunt)
name,opt=path+'psi','p'

Nmax=1
ls = np.hstack([[n]*(n+1) for n in np.arange(Nmax)])
ms = np.hstack([np.arange(n+1) for n in np.arange(Nmax)])

if 'i' in opts:
    phi,npts,R = 0, 100, 2
    xylims,xylims0 = 2*R*np.array([1e-3,1,0,1]),R*np.array([-1,1,-1,1])
    bbt = {'r_d':6,'theta_d':2*np.pi/3,'phi_d':np.pi, 'Nmax':10}
    for l,m in zip(ls,ms): plot_psi_lm(l=l,m=m,zl='in',zq=''  ,bbt={} ,phi=phi,cart_opt=True,npts=npts,r=xylims0,
        opt=opt,name=name,xylims=xylims0,imOpt='c',axPos='V')
    for l,m in zip(ls,ms): plot_psi_lm(l=l,m=m,zl='in',zq='in',bbt=bbt,phi=phi,cart_opt=True,npts=npts,r=xylims,
        opt=opt,name=name,xylims=xylims,imOpt='c',axPos='V')

if 'o' in opts:
    kd = 2
    bbt = {'r_d':kd,'theta_d':0*np.pi/3,'phi_d':0*np.pi, 'Nmax':10}
    xylims0 = np.hstack([[0,1], -kd+np.array([-1,1])]); r0=xylims0
    r,xylims = (3.,10),[0,4,-2,2]
    caxis = [0,0.25]
    for l,m in zip(ls,ms): plot_psi_lm(l=l,m=m,zl='out',zq=''   ,bbt={} ,phi=phi,cart_opt=True,npts=npts,r=r0,
        opt=opt,name=name,xylims=xylims0,imOpt='c',caxis=caxis,axPos='V')
    for l,m in zip(ls,ms): plot_psi_lm(l=l,m=m,zl='out',zq='out',bbt=bbt,phi=phi,cart_opt=False,npts=npts,r=r,
        opt=opt,name=name,xylims=xylims ,imOpt='c',caxis=caxis,axPos='V')

### wigner speed
if 'g' in opts:
    w3j = sage.functions.wigner.wigner_3j
    l,n,m,p = 4,10, 0,1
    qs = np.arange(abs(n-l),n+l)
    t = time.time()
    print([float(w3j(l,n,q,m,p,m-p)) for q in qs       ]);print('sage  = %.1E' %(time.time()-t));t=time.time()
    print(py3nj.wigner3j(2*l,2*n,2*qs,2*m,2*p,2*(m-p))  );print('py3nj = %.1E' %(time.time()-t));t=time.time()


kdp,n,Nmax = 2.0,1, 5
phi = np.pi/2
# Annu,Bnnu = spu.get_vec_trans_coeff_lin0(kdp,n,Nmax)
# print(np.abs(Annu).max(),np.abs(Bnnu).max())
# print('Annu:',Annu);print('Bnnu:',Bnnu)
# print('jn(n,ka)',[spu.jn(n,7) for n in range(1,Nmax+1)])

# Mstr = r"$\boldsymbol{M}_{%d1}^{(3)}(r,\theta,\phi=%d^{\circ})$" %(n,np.rad2deg(phi))
# fM = lambda kr,theta,phi:np.real(spu.get_MNn1(kr,theta,phi, Nmax=1, zn=spu.hn1,znp=spu.hn1p,pol='y')[0][0])
Mstr = r"$\boldsymbol{M}_{%d1}^{(3)}(r,\theta,\phi)$, $kd_p=$%.1f, $\nu_{max}=$%d" %(n,kdp,Nmax)
fM = lambda kr,theta,phi:spu.get_MNn1_t(kdp,n,kr,theta,phi,Nmax,zn=spu.jn,znp=spu.jnp,pol='y',fz=np.real)[0]
spu.plot_E3d_plane(fM,cmpts=[0,1,2],cart_opt='E',r=(0.5,5),npts=100,phi=phi,fieldname=Mstr,
    opt='ps',name=path+'Mt3_11',caxis=[-0.5,0.5])






# xylims=[-1,1]*3
# x,y,z = spu.mesh_3D(npts=5,r=(-1,1))
# kr,theta,phi = spu.cart2sphere(x,y,z)

# M = spu.Esphere2cart(Mr,theta,phi)
# C = M[0]*M[0]+M[1]*M[1]+M[2]*M[2]

# name=path+'vec_trans'
# cmap,mm = 'jet',C.max()
# fig,ax = dsp.stddisp(rc='3d',opt='')
# ax.quiver(x,y,z,M[0],M[1],M[2],color='r',length=0.25,linewidth=3)#,normalize=True)
# ax.scatter(x,y,z,c=C,vmin=0,vmax=mm,cmap=cmap)
# tle = r'$|| M_{%d1}^{(3)} ||^2$, $\nu_{max}=%d$' %(n,Nmax)
# nameE = '_M%d1_%d.png' %(n,Nmax)
# dsp.stddisp(fig=fig,ax=ax,labs=['$x$','$y$','$z$'],title=tle,imOpt='c',caxis=[0,mm],cmap=cmap,
#     name=name+nameE,xylims=xylims,axPos='V')#,**kwargs)
