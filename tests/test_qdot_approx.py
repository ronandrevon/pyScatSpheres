from utils import *
import matplotlib.patches as patches
import scipy,copy,time
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import spherical_utils as spu   ;imp.reload(spu)
from scipy.integrate import trapz,quad
import scipy.special as spe
plt.close('all')
ka,kp=6,1.01
N=10
nmax = 10
kd = 10
kds = np.arange(25,30,5)
opts = ''

qdot1 = qsa.QdotSphereArray(N=1,ka=ka,kp=kp,nmax=nmax,solve=0)
qdot1.solve(opt2=1)
# qdot1.test_convergence([10,11,12,13,15])
t0 = time.time()
qdot1 = qsa.QdotSphereArray(N=10,ka=ka,kp=kp,nmax=nmax,solve=1)
print(time.time()-t0)

if 'f' in opts:
    qdot2   = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=1)
    qdot2_0 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=0)
    qdot2.test_convergence([5,7,9])#,11])


    theta_d = np.linspace(0,60,2000)
    theta = np.deg2rad(theta_d)
    ff =  qdot2.get_ff(theta)
    ffu = qdot2_0.get_ff(theta)
    ffmax = abs(ff).max()
    plts = [[theta_d,abs(ff)/ffmax,'b-','exact'],[theta_d,abs(ffu)/ffmax,'r--','approx']]
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],lw=2)
    idx = abs(ff)/ffmax>1e-3
    err = abs((ff[idx]-ffu[idx])/ff[idx]).mean()
    print(err)
    # if err<0.05:qdot2.show_f()

if 'c' in opts:
    qdot2u = [qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=1) for kd in kds]

    cs = dsp.getCs('Spectral',nmax+1)

    ll = np.arange(nmax+1)
    jll = 1#1J**ll
    bp0 = qdot1.bp/jll
    bp1 = np.array([qdot.bp[:nmax+1] for qdot in qdot2u])/jll
    bp2 = np.array([qdot.bp[1*(nmax+1):2*(nmax+1)] for qdot in qdot2u])*np.exp(-1J*kds)[:,None]/jll
    bp3 = np.array([qdot.bp[2*(nmax+1):3*(nmax+1)] for qdot in qdot2u])*np.exp(-2J*kds)[:,None]/jll

    # for qdot in qdot2u:qdot.solve(optT=1)
    # bpT = np.array([qdot.bp[nmax+1:] for qdot in qdot2u])*np.exp(-1J*kds)[:,None]/jll

    # alnu = spu.a_ln(nmax,nmax,spu.hn1,kd,0,0)
    # bps = np.sum(qdot1.bp)
    vl = qdot1.vl
    bp2a = np.array([
        qdot1.bp*np.exp(1j*kd)+vl*np.sum(qdot1.bp*spu.a_ln(nmax,nmax,spu.hn1,kd,0,0).T,axis=1) for kd in kds])*np.exp(-1J*kds)[:,None]/jll
    bp3a = np.array([
        qdot1.bp*np.exp(2j*kd)+vl*np.sum(qdot1.bp*(spu.a_ln(nmax,nmax,spu.hn1,2*kd,0,0).T+spu.a_ln(nmax,nmax,spu.hn1,kd,0,0).T),axis=1) for kd in kds])*np.exp(-2J*kds)[:,None]/jll

    plts,legElt=[],{}
    plts+= [[bp0.real[i],bp0.imag[i],[cs[i],'o'],''] for i in range(nmax+1)]     ;legElt['$b_p^{(0)}$']=['k','o' ]
    # plts+= [[bp1.real[:,i],bp1.imag[:,i],[cs[i],'v-'],''] for i in range(nmax+1)];legElt['$b_1}$'     ]=['k','v-']
    # plts+= [[bp2.real[:,i],bp2.imag[:,i],[cs[i],'s-' ],''] for i in range(nmax+1)];legElt['$b_2}$'      ]=['k','s-']
    # plts+= [[bp2a.real[:,i],bpa.imag[:,i],[cs[i],'x--'],''] for i in range(nmax+1)];legElt['$b_2^{(a)}$' ]=['k','x--']
    plts+= [[bp3.real[:,i],bp3.imag[:,i],[cs[i],'s-' ],''] for i in range(nmax+1)];legElt['$b_3}$'      ]=['k','s-']
    plts+= [[bp3a.real[:,i],bp3a.imag[:,i],[cs[i],'x--'],''] for i in range(nmax+1)];legElt['$b_3^{(a)}$' ]=['k','x--']
    # plts+= [[bpT.real[:,i],bpT.imag[:,i],[cs[i],'d--'],''] for i in range(nmax+1)];legElt['$b_2^{(T)}$' ]=['k','d--']
    # legElt = {'$b_p^{(0)}$' :['k','o'],'$b_1$' :['k','d-'],
    #     '$b_2$' :['k','s-'],'$b_2^{(a)}$' :['k','x--']}
    dsp.stddisp(plts,labs=[r"$\Re(b_p)$",r"$\Im(b_p)$"],lw=2,legElt=legElt)

if 'C' in opts:
    nu = 1
    theta0  = np.linspace(0,np.pi,10)
    theta_p = np.hstack([theta0,np.flipud(theta0)])
    phi_p   = np.hstack([np.pi/2*np.ones(theta0.shape),-np.pi/2*np.ones(theta0.shape)])
    x,y,z   = spu.sphere2cart(ka*np.ones(theta_p.shape),theta_p,phi_p)
    kr,theta,phi =spu.cart2sphere(x,y,z+kd)

    psi_nu0 = spu.hn1(nu,kr)*spe.sph_harm(0,nu,phi,theta)

    ul,vl = qdot1.ul,qdot1.vl
    n_p  = qdot1.kp
    alnu = spu.a_ln(nu,nmax,spu.hn1,kd,0,0)[-1]
    # Yl0 = np.array([ spe.sph_harm(0,l,phi_p,theta_p) for l in range(nmax)])
    psi_out,psi_in=np.zeros(psi_nu0.shape,dtype=complex),np.zeros(psi_nu0.shape,dtype=complex)
    psi_nu0p=np.zeros(psi_nu0.shape,dtype=complex)
    for l in range(nmax):
        Yl0 = spe.sph_harm(0,l,phi_p,theta_p)
        psi_in   += alnu[l]*ul[l]*spu.jn(l,n_p*ka)*Yl0
        psi_out  += alnu[l]*vl[l]*spu.hn1(l,ka)*Yl0
        psi_nu0p += alnu[l]      *spu.jn(l,ka) *Yl0
    # scat = ([y,z+kd,psi_nu0.real,'o'],
    #         [y,z,psi_out.real,'s'])
            # (y,z,psi_in.real,psi_in.real,'d')]
    # dsp.stddisp(scat=scat)
    print(abs(psi_nu0+psi_out-psi_in).sum())
    # print(abs(psi_nu0-psi_nu0p).sum())


def draw(**kwargs):
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


    #patches arrows, and texts
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
    dsp.stddisp(texts=txt,fonts={'text':30},patches=pp,lw=2,pOpt='Xe',ax=ax,**kwargs)
if 'D' in opts:draw(opt='ps',name='../docs/figures/qdot2_approx_img.png')
