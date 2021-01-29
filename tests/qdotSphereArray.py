import importlib as imp
import scipy,copy
from optics.scattering import qdot_sphere_array as qsa ;imp.reload(qsa)
from optics.scattering import spherical_utils as spu   ;imp.reload(spu)
from optics.scattering import spherical_utils as sphu  ;imp.reload(sphu)
# from optics.scattering import gui_base as gui          ;imp.reload(gui)
from optics.scattering import gui_base as gui          ;imp.reload(gui)
from utils import *
from scipy.integrate import trapz,quad
import scipy.special as spe


plt.close('all')
path='../../docs/figures/'
nameB=path+'qdotSphereArray'


def solve_batch(df_name,kas,kps,kds=[0],N=1,nmax=4):
    kas,kps,kds = np.meshgrid(kas,kps,kds)
    kas,kps,kds = kas.flatten(),kps.flatten(),kds.flatten()
    cols = ['N','ka','kp','kd','nmax','sigma','ap','bp']
    df = pd.DataFrame(columns=cols)
    for ka,kp,kd in zip(kas,kps,kds):
        nmax = max(int(np.ceil(1.3*ka)),int(ka)+4)
        print(ka,kp,kd,nmax)
        s = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax,solve=1,copt=1)
        sig = s.get_s(npts=1000)
        df=df.append(dict(zip(cols,[s.N,s.ka,s.kp,s.kd,s.nmax,sig,s.ap,s.bp])),ignore_index=True)
    df.to_pickle(df_name)
    print(colors.green+'DataFrame saved:\n'+ colors.yellow+df_name+colors.black)


def singleSphere_f(kas,kp=1.5,nmax=4,fopts='m',npts=361,name='',**kwargs):
    '''Plot single sphere amplitude
    - fopts : r(real), i(imag), m(mag), a(angle), 2(mag2)
    '''
    s1 = [qsa.QdotSphereArray(N=1,ka=ka,kp=kp,kd=0,nmax=nmax,solve=1,copt=0) for ka in kas]
    nkas  = kas.size
    theta = np.linspace(1e-3,np.pi,npts)
    f = np.zeros((nkas,npts),dtype=complex)
    for l in range(nmax):
        Yl = spe.sph_harm(0,l,0,theta)
        for i,s in enumerate(s1):
            f[i] += s.bp[l]*(-1J)**(l+1)*Yl

    # for i,ka in enumerate(kas):f[i]/=ka
    cs = dsp.getCs('jet',nkas)
    theta = np.rad2deg(theta)
    csb,csr,csg = dsp.getCs('Blues',nkas),dsp.getCs('Reds',nkas),dsp.getCs('Greens',nkas)
    plts=[]
    if 'r' in fopts:plts += [[theta,np.real(f[i])  ,[csb[i],'--'],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
    if 'i' in fopts:plts += [[theta,np.imag(f[i])  ,[csr[i],'--'],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
    if 'a' in fopts:plts += [[theta,np.angle(f[i]) ,[csg[i],'- '],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
    if 'm' in fopts:plts += [[theta,np.abs(f[i])   ,[csb[i],'- '],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
    if '2' in fopts:plts += [[theta,np.abs(f[i])**2,[csr[i],'- '],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,#logOpt='y',
        title='Scattering amplitudes $n_{ref}=%.4f$ for different $ka$' %kp,xylims=['x',0,180],#,1e-2,1e2],
        name=name+'fka.svg',**kwargs)
    return f

def singleSphere_s(df_name,name='',**kwargs):
    df = pd.read_pickle(df_name)
    kas,kps = df.ka.unique(),df.kp.unique()

    plts,cs = [],dsp.getCs('Spectral',kps.size)
    for ik,kp in enumerate(kps):
        dfn = df[df.kp==kp]
        sig_ka = np.array([np.sum(np.abs(r.bp)**2)/r.ka**2 for i,r in dfn.iterrows()]) #;print(sig_th)
        plts += [[kas,sig_ka,cs[ik],'$n_{ref}=%.4f$' %kp]]
    dsp.stddisp(plts,labs=[r"$ka$",r"$\sigma/(ka)^2$"],lw=2,#xylims=['x',kas[0],kas[-1]],
        name=name+'ska.svg',**kwargs)

    # sig = np.array([np.sum(np.abs(s.bp)**2) for i,s in df.iterrows()]) #;print(sig_th)
    # df['sigma'] = sig
    # df.to_pickle(df_name)
    # print(colors.green+'DataFrame saved:\n'+ colors.yellow+df_name+colors.black)

def plot_sigma(df_name,ika=slice(None,None,2),name='',**kwargs):
    df = pd.read_pickle(df_name)
    kas,kds = df.ka.unique(),df.kd.unique()

    plts,cs = [],dsp.getCs('Spectral',kas.size)
    for ik,ka in enumerate(kas[ika]):
        dfn = df[df.ka==ka]
        plts += [[dfn.kd/ka,dfn.sigma/ka**2,cs[ik],'$ka=%.1f$' %ka]]
    dsp.stddisp(plts,labs=[r"$kd/ka$",r"$\sigma/(ka)^2$"],lw=2,#xylims=['x',kas[0],kas[-1]],
        name=name+'ska.svg',**kwargs)




##############################################################################
#### N=1 tests
##############################################################################
name=nameB+'1_'
opts = '' # 'n'(near field) 's'(sigma) 'f'(amplitude) 'g'(gui)
####near field
if 'n' in opts:
    ka = 3.0
    r,npts  = ka*np.array([1e-3,2,-1.5,1.5]),400
    s1 = qsa.QdotSphereArray(N=1,ka=ka,kp=1.5,kd=0,nmax=10,solve=1)
    s1.show_f(cart_opt=True,npts=npts,r=r,xylims=r,opts='istTG',fz=np.real,caxis=[-1,1],name=name,opt='sc')
    s1.show_f(cart_opt=True,npts=npts,r=r,xylims=r,opts='istT' ,fz=np.real,caxis=[-1,1],name=name,opt='sc')
####far field plots
if 'f' in opts:
    # kps = np.array([np.sqrt(1+sig)])
    kps = np.array([1.0001,1.01,1.1,1.2])
    # kas = np.array([10])
    kas = np.array([0.5,2,5,10,30])
    for i,kp in enumerate(kps):
        f1 = singleSphere_f(kas,kp,nmax=30,fopts='m',name=name+'n%d_'%i,opt='c')
#### scattering cross section plot
if 's' in opts:
    df_name = 'data/qdotSphereArray1.pkl'
    kps = np.array([1.0001,1.01,1.1,1.2])
    kas = np.linspace(0.1,40,100)
    # solve_batch(df_name,kas,kps,nmax=40)
    s1 = singleSphere_s(df_name,name=name+'log_',logOpt='y',xylims=[0,40,1e-8,20],opt='p')

    df_name = 'data/qdotSphereArray1_res.pkl'
    kps = np.array([1.1,1.2,1.5,3])
    kas = np.linspace(0.1,3,150)
    # ssolve_batch(df_name,kas,kps,nmax=10)
    # s1 = singleSphere_s(df_name,name=name+'lin_',opt='p')
#### gui
if 'g' in opts:gui.GUI_handler(df_name = 'data/qdotSphereArray1.pkl')



##############################################################################
#### N=2 tests
##############################################################################
name=nameB+'2_'
opts = '' #1(N=1 test) u(unique)M(Matrix)F(Field) c(convergence) s(solve batch) g(gui)

if '1' in opts:
    s1 = qsa.QdotSphereArray(N=1,ka=1.,kp=1.5,kd=0.5,nmax=4,solve=True,copt=0)
    s0 = qsa.QdotSphereArray(N=1,ka=1.,kp=1.5,kd=0.5,nmax=4,solve=True,copt=1)
    print(s1.ap);print(s0.ap)
    print(s1.bp);print(s0.bp)


if 'u' in opts:
    ka,kp,kd,nmax = 4,1.1,12, 10
    s1 = qsa.QdotSphereArray(N=2,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=True,copt=1)
    if 'M' in opts:
        h,k = np.meshgrid(range(20),range(20))
        M = np.log10(np.abs(s1.P)+1)
        M = np.log10(np.abs(s1.T)+1)
        M = np.log10(np.abs(s1.P-s1.T)+1)
        dsp.stddisp(im=[h,np.flipud(k),M],imOpt='c',axPos='V',xyTicks=[1,1],xylims=[0,20,0,20],caxis=[0,M.max()],cmap='Blues')
    if 'F' in opts:
        r = (1e-3,2*ka,-2*ka,kd+2*ka)
        # s1.show_f(r=r,xylims=r)
        # s1.show_f(opts='sP',idp=[0],r=r,xylims=r,fz=np.real,caxis=[-1,1])
        # s1.show_f(opts='sP',idp=[1],r=r,xylims=r,fz=np.real,caxis=[-1,1])
        # s1.show_f(opts='tT',r=r,xylims=r,fz=np.real,caxis=[-1,1],name=name,opt='p')
        # s1.show_f(opts='tTG',r=r,xylims=r,fz=np.real,caxis=[-1,1])
        # s1.show_ff(xylims=['x',0,180])
    if 'P' in opts:
        s1 = qsa.QdotSphereArray(N=2,ka=1,kp=1.5,kd=5,nmax=5,solve=True,copt=1)
        r = (1e-3,2*s1.ka,-2*s1.ka,s1.kd+2*s1.ka)
        s1.show_f(opts='tT',r=r,xylims=r,fz=np.real,caxis=[-1,1],name=name+'0_',opt='ps')
        s1.show_ff(name=name+'0_',opt='ps')
        s1 = qsa.QdotSphereArray(N=2,ka=5,kp=1.1,kd=12.5,nmax=8,solve=True,copt=1)
        r = (1e-3,2*s1.ka,-2*s1.ka,s1.kd+2*s1.ka)
        s1.show_f(opts='tT',r=r,xylims=r,fz=np.real,caxis=[-1,1],name=name+'1_',opt='ps')
        s1.show_ff(name=name+'1_',opt='ps')
        s1 = qsa.QdotSphereArray(N=2,ka=15,kp=1.01,kd=35,nmax=20,solve=True,copt=1)
        r = (1e-3,2*s1.ka,-2*s1.ka,s1.kd+2*s1.ka)
        s1.show_f(opts='tT',r=r,xylims=r,fz=np.real,caxis=[-1,1],name=name+'2_',opt='ps')
        s1.show_ff(name=name+'2_',opt='ps')
## convergence test
if 'c' in opts:
    nmaxs = np.unique(np.array(np.ceil(np.array([1.25,1.4,1.5])*ka),dtype=int))
    s1 = qsa.QdotSphereArray(N=2,ka=ka,kp=kp,kd=kd*ka,nmax=15,solve=False,copt=1)
    s1.test_convergence(nmaxs=nmaxs,name=name,opt='p')

if 's' in opts:
    kas = np.linspace(0.5,20,20)
    kds = np.linspace(2,8,20)
    kps = np.array([1.0001,1.01,1.1])
    for i,kp in enumerate(kps[:1]):
        df_name = 'data/qdotSphereArray2_kp%d.pkl' %i
        # solve_batch(df_name,kas,[kp],kds,N=2)
        plot_sigma(df_name,ika=slice(None,None,2),name=name+'kp%d' %i,opt='p')
if 'g' in opts:
    hg = gui.GUI_handler(df_name='data/qdotSphereArray2_kp2.pkl')#,hardSphere=False)


##############################################################################
#### N=2 Approx comparison
##############################################################################
name=nameB+'2approx_'
opts = 'N'
if 'd' in opts:
    N,ka,kp,nmax = 5,2,1.2,5
    kds = [2.5,4,10]
    for i,kd in enumerate(kds):
        s2 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax,solve=True)
        s0 = copy.copy(s2); s0.solve(copt=0)
        fig,ax = s0.show_ff(fopts='m',opt='',leg=0)
        s2.show_ff(ax=ax,fopts='m',legElt={'uncoupled':'k--'},name=name+'kd%d_' %i,opt='sc')

if 'N' in opts:
    ka,kd,kp,nmax = 3,2.5,1.005,5
    Ns = [1,10,20,30,50]
    s2s,s0s = [],[]
    for i,N in enumerate(Ns):
        s2 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax);s2.solve(v=1)
        s0 = copy.copy(s2); s0.solve(copt=0)
        fig,ax = s0.show_ff(npts=2000,fopts='m',opt='',leg=0)
        s2.show_ff(npts=2000,ax=ax,fopts='m',legElt={'uncoupled':'k--'},name=name+'N%d_' %i,opt='p')
        s2s+=[s2]
        s0s+=[s0]
