import importlib as imp
from optics.scattering import hard_sphere_arrayApprox as hsaA   ;imp.reload(hsaA)
from optics.scattering import hard_sphere_array as hsa          ;imp.reload(hsa)
from optics.scattering import gui as hsa_gui                    ;imp.reload(hsa_gui)
from optics.scattering import spherical_utils as spu            ;imp.reload(spu)
from optics.scattering import spherical_utils as sphu           ;imp.reload(sphu)
import scipy.special as spe
from utils import *
plt.close('all')
path='../../docs/figures/'
impath='../../docs/resources/hamid1990/'
name=path+'hardSphereArray_'


def plot_hkl1(kd=np.linspace(1,30,300),ls=np.arange(10),fz=np.abs,**kwargs):
    hlinf = lambda l,kd:1j**(l+1)*np.exp(1J*kd)/kd
    cs = dsp.getCs('Spectral',ls.size)
    plts  = [[kd,fz(spu.hn1(l,kd)),cs[i],'$l=%d$' %l] for i,l in enumerate(ls)]
    plts += [[kd,fz(hlinf(l,kd))  ,[cs[i],'--'],''   ] for i,l in enumerate(ls)]

    dsp.stddisp(plts,labs=['$kd$','$h_l(kd)$'],lw=2,**kwargs)

def test_addth(opts='ht'):
    lmax, nmax = 4,10
    kd,theta_d = 1.5,1*np.pi
    ap = 0.5
    npts = 200

    t = np.linspace(-np.pi/2,np.pi/2,100)
    ct,st = np.cos(t), np.sin(t)
    ll = np.arange(lmax)
    mm = np.ones((lmax))
    z_q,z_r,fz_q = spu.hn1,spu.jn,spu.hn1
    # hankel in pth sphere
    if 'h' in opts:
        kr_p,theta_p,yp,zp = spu.polar_mesh2D(npts=npts,r=(ap,3),cart_opt=0)
        psi = np.array([z_q(l,kr_p)*spe.sph_harm(0,l,0,theta_p) for l in range(lmax)])
        mm = [np.abs(psi[l]).max() for l in ll]

        plts = [[ap*ct, ap*st,'k-','']]
        for l in ll:dsp.stddisp(plts,im=[yp,zp,np.abs(psi[l])], lw=2,labs=['$y_p$','$z_p$'],
            xylims=1.5*np.array([0,1,-1,0]),imOpt='c',axPos='E',caxis=[0,mm[l]],
            title=r"$\psi_{%d}(r_p,\theta_p)$" %l,fonts={'title':30},
            name=name+'_psi%d.png' %l, opt='p')

    # addition theorem
    if 't' in opts:
        aln = np.zeros((lmax+1,nmax),dtype=complex)
        for l in range(lmax):aln[l,:] = [spu.a_lmnp(l,0,n,0, fz_q,kd,theta_d,0) for n in range(nmax)]
        aln = a_ln(lmax,nmax, fz_q,kd,theta_d)

        kr,theta,y,z = spu.polar_mesh2D(npts=npts,r=(1e-3,1.5),cart_opt=1)
        zn = np.array([z_r(n,kr) for n in range(nmax)])
        Yn  = np.array([spe.sph_harm(0,n,0,theta) for n in range(nmax)])
        psi_l = np.zeros((lmax,)+kr.shape,dtype=complex)
        for l in range(lmax):
            for n in range(nmax):
                psi_l[l] += aln[l,n]*zn[n]*Yn[n]

        plts = [ [ap*ct, kd+ap*st,'k-',''] ]
        for l in ll:dsp.stddisp(plts,im=[y,z,np.abs(psi_l[  l])], lw=2,labs=['$y$','$z$'],
            xylims=1.5*np.array([0,1,0,1]),imOpt='c',axPos='E',caxis=[0,mm[l]],
            title=r"$\psi_{%d}(r,\theta)$ $N_{max}=%d$" %(l,nmax),fonts={'title':30},
            name=name+'_psi%d_t.png' %l, opt='p')

################################################################################
#### single sphere tests
################################################################################
def singleSphere(kas,nmax=4,idl=[],ida=1,opts='csf',fopts='m',name='',**kwargs):
    '''Plot single sphere
    - opts : c(cp), f(scattering amplitudes), s(scattering cross section)
    - fopts : r(real), i(imag), m(mag), a(angle), 2(mag2)
    '''
    s1 = [hsa.HardSphereArray(N=1,ka=ka,kd=0,nmax=nmax) for ka in kas]
    cp0 = np.array([s.solve(copt=0,v=0) for s in s1])

    if 'f' in opts:
        ida=np.arange(0,kas.size,ida)
        nkas = ida.size

        npts  = 361
        theta = np.linspace(0,np.pi,npts)
        f = np.zeros((nkas,npts),dtype=complex)
        for l in range(nmax):
            Yl = spe.sph_harm(0,l,0,theta)
            for i,ia in enumerate(ida):
                f[i] += s1[ia].Cp[l]*(-1J)**(l+1)*Yl

        # for i,ia in enumerate(ida):f[i]/=kas[ia]
        cs = dsp.getCs('jet',nkas)
        theta = np.rad2deg(theta)
        csb,csr,csg = dsp.getCs('Blues',nkas),dsp.getCs('Reds',nkas),dsp.getCs('Greens',nkas)
        plts=[]
        if 'r' in fopts:plts += [[theta,np.real(f[i])  ,[csb[i],'--'],"$ka=%.1f$" %kas[ia]] for i,ia in enumerate(ida)]
        if 'i' in fopts:plts += [[theta,np.imag(f[i])  ,[csr[i],'--'],"$ka=%.1f$" %kas[ia]] for i,ia in enumerate(ida)]
        if 'a' in fopts:plts += [[theta,np.angle(f[i]) ,[csg[i],'- '],"$ka=%.1f$" %kas[ia]] for i,ia in enumerate(ida)]
        if 'm' in fopts:plts += [[theta,np.abs(f[i])   ,[csb[i],'- '],"$ka=%.1f$" %kas[ia]] for i,ia in enumerate(ida)]
        if '2' in fopts:plts += [[theta,np.abs(f[i])**2,[csr[i],'- '],"$ka=%.1f$" %kas[ia]] for i,ia in enumerate(ida)]
        dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,logOpt='y',
            title='Scattering amplitudes for different $ka$',xylims=[0,180,0.09,1e3],
            name=name+'fka.svg',**kwargs)

    if 's' in opts:
        sig_th = np.sum(np.abs(cp0)**2,axis=1)
        plts = [[kas,sig_th/kas**2,'b','']]
        # numerical integral checled was identical
        # sig_nm = np.array([s.get_s(npts=361) for s in s1])
        # plts+= [[kas,sig_nm/kas**2,'b--','']]
        dsp.stddisp(plts,labs=[r"$ka$",r"$\sigma/(ka)^2$"],lw=2,logOpt='x',xylims=[kas[0],kas[-1],0,13],
            name=name+'ska.svg',**kwargs)


    if 'c' in opts:
        idl=np.arange(0,nmax,idl)
        # if not idl.size:idl=np.arange(nmax)
        cs = dsp.getCs('Spectral',idl.size)
        plts  = [[ cp0[:,l].real ,cp0[:,l].imag,[cs[i],'-o' ],'$l=%d$' %l] for i,l in enumerate(idl)]
        dsp.stddisp(plts,labs=["$Re$","$Im$"],lw=2,
            title='$c_{p;l}(ka=%.1f-%.1f)$' %(kas[0],kas[-1]),
            name=name+'cpl.svg',**kwargs)

    return s1


################################################################################
##### multiple sphere tests
################################################################################
def sphere_array2(N=2,nmax=7,ka=2.0,kd=8,opts='T',npts=300,nZ=(1,1),nY=1, **kwargs):
    ''' plot the near field amplitudes
    - opts : T(total field), P(field individual spheres)
    - nR   : intor 2-tuple : extent of plot before and after the spheres in ka units (normalized radius)
    '''
    if isinstance(nZ,int) or isinstance(nZ,float):nZ = (nZ,nZ)
    s2 = hsa.HardSphereArray(N,ka,kd,nmax)
    r = [0,(1+nY)*ka,-(nZ[0]+1)*ka,kd+(1+nZ[1])*ka]
    if 'P' in opts:
        Cp = s2.solve(copt=0)
        s2.show_f(cart_opt=True,npts=npts,r=r,opts=opts.replace('T',''),fz=np.real,
            xylims=r,caxis=[-1,1],
            **kwargs)
    if 'T' in opts:
        Cp = s2.solve(copt=1)
        s2.show_f(cart_opt=True,npts=npts,r=r,opts=opts.replace('P',''),fz=np.real,
            xylims=r,caxis=[-1,1],
            **kwargs)
    return s2

def sphere_array_convergence(nmaxs=[2,3,4,5,7,10],N=2,ka=2.0,kd=8,npts=361,name=name, **kwargs):
    theta = np.linspace(0,np.pi,npts)
    theta_d = np.rad2deg(theta)
    ns = np.array(nmaxs).size
    f = np.zeros((ns,npts),dtype=complex)
    s2 = hsa.HardSphereArray(N,ka,kd)
    for i,nmax in enumerate(nmaxs):
        print('solving nmax=%d' %nmax)
        s2.solve(nmax=nmax,copt=1,v=0)
        f[i] = s2.get_ff(theta)

    cs = dsp.getCs('Spectral',ns)
    plts = [[theta_d,np.abs(f[i]), cs[i] ,'$n_{max}=%d$' %nmax] for i,nmax in enumerate(nmaxs)]
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,
        title='convergence test $ka=%.1f$, $kd=%.1f$' %(ka,kd),xylims=['x',0,180],
        name=name+'fconv.svg',**kwargs)

def sweep_ka(N=2,nmax=7,kas=[0.5,1,2,4],nkds=50,kdr=(2,10),npts=361,Nmax=5,fz=np.abs,opts='',name='', **kwargs):
    '''
    - opts : 'c'(Cp(kd;ka)), 's'(sigma) 'f'(scattering amplitudes) A(approx)
    '''
    kas = np.array(kas)
    nkas = kas.size
    theta = np.linspace(0,np.pi,npts)
    theta_d = np.rad2deg(theta)
    ct,st,dt = np.cos(theta),np.sin(theta), theta[1]-theta[0]

    sig,sig0 = np.zeros((nkas,nkds)),np.zeros((nkas,nkds))
    # f  = np.zeros((nkas,nkds,npts),dtype=complex)
    # f0 = np.zeros((nkas,nkds,npts),dtype=complex)
    for i,ka in enumerate(kas):
        kds = np.linspace(max(kdr[0],2*ka),kdr[1],nkds)
        # kds = np.linspace(kdr[0]*ka,kdr[1]*ka,nkds)
        ss = [hsa.HardSphereArray(N,ka,kd,nmax) for kd in kds]
        print('...%d : solving ka=%.1f...' %(i,ka))
        if 'A' in opts:
            ssA = [hsaA.HardSphereArray(N,ka,kd,nmax,solve=True) for kd in kds]
            cp0 = np.array([s._cpl() for s in ssA])
        else:
            cp0 = np.array([s.solve(copt=0,v=0) for s in ss])
        cp  = np.array([s.solve(copt=1,v=0) for s in ss])


        if 'c' in opts:
            print('...plotting coefficients...')
            plotCp(ss,cp,cp0,Nmax=Nmax,name=name+'ka%d' %i, **kwargs)

        if 'f' in opts or 's' in opts:
            print("...Computing scattering amplitudes...")
            f,f0 = np.zeros((nkds,npts),dtype=complex),np.zeros((nkds,npts),dtype=complex)
            for l in range(nmax):
                Yl = spe.sph_harm(0,l,0,theta)
                for id,kd in enumerate(kds):
                    for p in range(N):
                        idp = p*nmax
                        ep = np.exp(-1J*kd*p*ct)
                        f[id]  +=  cp[id,idp+l]*ep*(-1J)**(l+1)*Yl
                        f0[id] += cp0[id,idp+l]*ep*(-1J)**(l+1)*Yl
            if 'f' in opts:
                print("...plotting scattering amplitude...")
                csD = dsp.getCs('Spectral',nkds)
                for id,kd in enumerate(kds):
                    plts = [[theta_d,fz( f[id]), csD[id]      ,'kd=%.1f' %kd] for id,kd in enumerate(kds)]
                    plts+= [[theta_d,fz(f0[id]),[csD[id],'--'],''           ] for id,kd in enumerate(kds)]
                    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,
                        title='Scattering amplitudes for $ka=%.1f$' %ka,xylims=['x',0,180],
                        name=name+'fka%d.svg' %i,**kwargs)

            if 's' in opts:
                # for l in range(nmax):
                #     for p in range(N):
                #         idq = np.arange(p+1,N)*nmax
                        # sig[i,:] += np.abs(cp[:,l+nmax*p])**2 ## sig_pI = cp[:,l+nmax*p]*np.sum(np.conj(cp[:,l+idq]),axis=1);print(sig_pI.imag.max(()))
                        # sig[i,:] += 2*np.real(cp[:,l+nmax*p]*np.sum(np.conj(cp[:,l+idq]),axis=1))
                print("...integrating cross section...")
                for id,kd in enumerate(kds):
                    # print(f[id].shape,sig[i,id])
                    # print(np.sum(np.abs(f[id,:])**2))
                    sig[i,id]  = 2*np.pi*np.sum(np.abs( f[id,:])**2*st)*dt
                    sig0[i,id] = 2*np.pi*np.sum(np.abs(f0[id,:])**2*st)*dt
    if 's' in opts:
        print("...display cross sections")
        # print(sig);print(sig0)
        cs,plts = dsp.getCs('Spectral',nkas),[]
        for i,ka in enumerate(kas):
            plts += [[kds, sig[i]/ka**2,[cs[i],'-o' ],'$ka=%.1f$' %ka]]
            plts += [[kds,sig0[i]/ka**2,[cs[i],'--o'],'']]
        dsp.stddisp(plts,labs=[r"$kd$",r"$\sigma/(ka)^2$"],lw=2,
            #logOpt='x',#xylims=[kas[0],kas[-1],0,13],
            name=name+'ska.svg',**kwargs)

    return sig

def plotCp(ss,cp,cp0,Nmax=5,name=name,**kwargs):
    ''' display trajectory of Cp in complex plane for each sphere'''
    nds,N,nmax,ka = len(ss),ss[0].N,ss[0].nmax, ss[0].ka
    kds = np.array([s.kd for s in ss])

    # print('...plotting coefficients...')
    cs = dsp.getCs('Spectral',nmax)
    for p in range(N):
        idp = p*nmax
        plts  = [[ cp[:,idp+n].real ,cp[:,idp+n].imag,[cs[n],'-o' ],'$l=%d$' %n] for n in range(Nmax)]
        plts += [[cp0[:,idp+n].real,cp0[:,idp+n].imag,[cs[n],'--s'],''] for n in range(Nmax)]
        # print('Cp=%d' %p)
        dsp.stddisp(plts,labs=["$Re$","$Im$"],lw=2,
            title='$c_{p;l}(ka=%.1f,kd=%.1f-%.1f)$ sphere p=%d' %(ka,kds[0],kds[-1],p),
            name=name+'cpl%d.svg' %p,**kwargs)


################################################################################
##### multiple sphere tests
################################################################################
def sphere_array2Approx(N=2,nmax=7,ka=2.0,kd=8,fopts='ue', **kwargs):
    '''
    fopts : 'u' (uncoupled), e'(exact)
    '''
    theta  = np.linspace(0,np.pi,361)

    s2a = hsaA.HardSphereArray(N,ka,kd,nmax)
    s2  = hsa.HardSphereArray(N,ka,kd,nmax)

    cp  =  s2.solve(copt=1)#;print(cp)
    cpa = s2a.solve(copt=1)#;print(cpa)
    fe  = s2.get_ff(theta)
    fa = s2a.get_ff(theta)


    if 'u' in fopts:
        cp  =  s2.solve(copt=0)#;print(cp)
        cpa = s2a.solve(copt=0)#;print(cpa)
        fe0 = s2.get_ff(theta)
        fa0 = s2a.get_ff(theta)


    theta_d = np.rad2deg(theta)
    plts= [[theta_d,np.abs(fa) ,'r-' ,'$f^{approx}_{coupled}$']]
    if 'u' in fopts:plts+= [[theta_d,np.abs(fa0),'r--','$f^{approx}_{uncoupled}$']]
    if 'e' in fopts:
        plts+= [[theta_d,np.abs(fe0),'b--','$f^{exact}_{uncoupled}$']]
        if 'u' in fopts:plts+= [[theta_d,np.abs(fe) ,'b-' ,'$f^{exact}_{coupled}$']]


    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,
        title='Scattering amplitudes for $ka=%.1f$, $kd=%.1f$' %(ka,kd),xylims=['x',0,180],
        **kwargs)

    return s2a



########################################################################################################################
# test_addth(opts='ht')
# plot_hkl1(kd=np.linspace(1,50,300),ls=np.arange(0,51,5),fz=np.abs,xylims=[1,50,1e-2,1e3],logOpt='y',
#     name=path+'hardSphereArray1_hl.svg',opt='ps')
df_name = 'data/hard_sphere_array2.pkl'
# df = hsa.sweep_ka(df_name,N=2,nmax=15,kas=[0.5,1,2,4,10],kds=np.linspace(1,30,150))
# df = hsa.sweep_ka(df_name,N=2,nmax=10,kas=[0.5,1,2,4],kds=None,nkds=2,kdr=(2,10),v=1)


name=path+'hardSphereArray1_'
# s1 = singleSphere(kas=np.logspace(-1,1.75,101),nmax=60,opts='s',ida=10,fopts='m',
#     name=name,opt='p')
# s1 = singleSphere(kas=np.linspace(30,40,50) ,nmax=50,idl=5,opts='c' ,
#     name=name, opt='ps')
# s1 = singleSphere(kas=np.linspace(0.5,5,50) ,nmax=10,idl=1,opts='c' ,
#     name=name+'1', opt='ps')
# s1 = singleSphere(kas=np.linspace(5,20,50) ,nmax=30,idl=4,opts='c' ,
#     name=name+'2', opt='ps')

name=path+'hardSphereArray2_'
# s2 = sphere_array2(N=2,nmax=7,ka=2.5,kd=8.6,opts='stT',npts=600,nZ=(1,2),nY=4, name=name,opt='p')
# sphere_array_convergence(nmaxs=np.arange(1,6),N=2,ka=0.5,kd=10,npts=361,name=name, opt='p')
# sphere_array_convergence(nmaxs=np.arange(3,22,3),N=2,ka=10.0,kd=20,npts=361,name=name, opt='p')
# ss = sweep_ka(N=2,nmax=7 ,kas=[0.5,1,2,4]                , opts='c',name=name ,opt='ps')
# ss = sweep_ka(N=2,nmax=10,kas=[1,2,5],kdr=(10,30),nkds=50, opts='s',name=name ,opt='ps')
# ss = sweep_ka(N=2,nmax=15,kas=[2,5,10],kdr=(2,10),nkds=3 , opts='f',name=name ,opt='ps')

### approx
ka,kd,nmax=2,10,5
name=path+'hardSphereArray2approx_'
## s1=hsaA.HardSphereArray(N=1,ka=ka,nmax=nmax);s1.solve(copt=0);s1.show_ff(npts=361,fopts='rim',opt='p')
## s2a = hsaA.HardSphereArray(2,ka,kd,nmax,solve=True);print(s2a.L)
# name=path+'hardSphereArray2approx_'
# s2a=sphere_array2Approx(N=2,ka=1.0,kd=20,nmax=10,opts='ue' ,name=name+'fka0.svg',opt='ps')
# s2a=sphere_array2Approx(N=2,ka=5.0,kd=40,nmax=40,fopts='ue',name=name+'fka1.svg',opt='ps')
# ss = sweep_ka(N=2,nmax=20 ,kas=[2.0],nkds=10,kdr=(10,11),Nmax=5,opts='cA',name=name,opt='p')

# plt.close('all')

s1=hsa.HardSphereArray(N=10,ka=2,kd=10,nmax=5,solve=True);
s1.show_ff(fopts='m',opt='p')
s1.show_f(opts='tT',npts=500,fz=np.real,caxis=[-1,1],name=name+'10_',opt='p')
s1.show_f(opts='sT',npts=500,fz=np.real,caxis=[-1,1],name=name+'10_',opt='p')
# s1.test_convergence(nmaxs=[3,4,5])
