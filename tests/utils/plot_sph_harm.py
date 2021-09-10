import importlib as imp
from utils import*
from optics.scattering.spherical_utils import*
import scipy.special as spe
import optics.scattering.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp            ;imp.reload(dsp)
path='../../docs/figures/'


def plot_spherical_Bessel(opt='j'):
    r = np.linspace(0.01,20,300)
    N = range(10)#[0,1]#,2,5,10]#,20,30]
    plts,cs =[], dsp.getCs('Spectral',len(N))
    if 'j' in opt:
        psi_jn  = [jn(n,r) for n in N]
        psi_jnp  = [spu.jnp(n,r) for n in N]
        # psi_jnp1 = [np.gradient(jn(n,r),r) for n in N]

        pltsJn = [[r,np.real(psi_jn[i] ),[cs[i],'-'],'$n=%d$' %n] for i,n in enumerate(N)]
        # pltsJn = [[r,np.real(psi_jnp[i]),[cs[i],'--'],''] for i,n in enumerate(N)]
        # pltsJn+= [[r,np.real(psi_jnp1[i]),[cs[i],'o'],''] for i,n in enumerate(N)]
        dsp.stddisp(pltsJn,lw=2,labs=['$r$','$j_n(r)$'],xylims=[-1,20,-0.5,1.1])

    if 'h' in opt:
        psi_hn1a = [spu.hn1(n,r) for n in N]
        psi_hn1b = [spu.jn(n,r)+1J*spe.spherical_yn(n,r) for n in N]


        plts += [[r,Hn[i].real,[cs[i],'--'],'$n=%d$' %n] for i,n in enumerate(N)]
        plts += [[r,Hn[i].imag,[cs[i],'-.'],''] for i,n in enumerate(N)]
        pltsHn = [[r,np.abs(psi_hn1a[i]),[cs[i],'-'],'$n=%d$' %n] for i,n in enumerate(N)]
        pltsHn+= [[r,np.abs(psi_hn1b[i]),[cs[i],'--'],''] for i,n in enumerate(N)]
        dsp.stddisp(pltsHn,lw=2,labs=['$r$','$h_n^{(1)}(r)$'],xylims=[-1,20,0,10])

    if 'y' in opt:
        psi_yn  = [spe.spherical_yn(n,r) for n in N]
        pltsYn = [[r,np.real(psi_yn[i]),[cs[i],'-'],'$n=%d$' %n] for i,n in enumerate(N)]
        dsp.stddisp(pltsYn,lw=2,labs=['$r$','$y_n(r)$'],xylims=[-1,20,-1.65,0.4])



def plot_real_sh_harm(m,n,npts=100,**kwargs):
    '''plots REAL Spherical harmonics from Ynm'''
    theta, phi = np.meshgrid(np.linspace(0,np.pi,npts),np.linspace(0,2*np.pi,npts))
    if m==0 :
        ymn = np.real(spe.sph_harm(m,n,phi,theta))
    elif m>0:
        ymn = np.sqrt(2)*(-1)**m*np.real(spe.sph_harm(m,n,phi,theta))
    elif m<0:
        ymn = np.sqrt(2)*(-1)**m*np.imag(spe.sph_harm(-m,n,phi,theta))
    r = np.abs(ymn)
    c = ( np.sign(ymn) +1)/2
    x,y,z = spu.sphere2cart(r,theta,phi)
    y_str=r'$Y_{%d}^{%d}$' %(n,m)
    dsp.stddisp(scat=[x,y,z,c],labs=['$x$','$y$','$z$'],title='real spherical harmonics %s' %y_str,rc='3d',**kwargs)

def plot_vec_sh_harm(Nmax,npts=121,ntp=15,MNopt=False,Ymn_opt=False,**kwargs):
    '''plots vectorial Spherical harmonics from Ulm,Vlm'''
    theta, phi = np.linspace(0,np.pi,npts),np.linspace(0,2*np.pi,npts)
    idx = slice(0,None,ntp),slice(0,None,ntp)
    theta3d,phi3d = np.meshgrid(theta,phi)
    theta_3d,phi_3d = theta3d[idx],phi3d[idx]
    ct,st,cp,sp = np.cos(theta_3d),np.sin(theta_3d),np.cos(phi_3d),np.sin(phi_3d)
    jmPlm,dtPlm = spu.get_jmdtPlm(Nmax,theta3d,phi3d,Ymn_opt,v=1)


    print('plot...')
    for l in range(1,Nmax+1):
        idl = int(l*(l+1)/2)
        for m in range(l+1):
            print('idlm:',idl+m)
            Ut = np.real(jmPlm[idl+m])
            Up = np.real(-dtPlm[idl+m])

            r3d = np.sqrt(Ut**2+Up**2)
            x,y,z = sphere2cart(r3d,theta3d,phi3d)
            fig,ax = dsp.stddisp(rc='3d',opt='',
                surfs=[[x,y,z,{'color':'b','alpha':0.65}]])

            if ntp :
                r_3d,Ut_3d,Up_3d = r3d[idx],Ut[idx],Up[idx]
                x,y,z = sphere2cart(r_3d,theta_3d,phi_3d)
                Ex =  Ut_3d*ct*cp-sp*Up_3d
                Ey =  Ut_3d*ct*sp+cp*Up_3d
                Ez = -Ut_3d*st
                ax.quiver(x,y,z,Ex,Ey,Ez,color='r',linewidth=3,normalize=True,length=0.25),#scale=0.5)

            y_str=r'$U_{%d}^{%d}$' %(l,m)
            dsp.standardDisplay(ax=ax,labs=['$x$','$y$','$z$'],title='vector spherical harmonics %s' %y_str,
            name=path+'UV%d%d.png' %(l,m),**kwargs)


def plot_MN(Nmax,npts=100,phi=0,Ymn_opt=False,**kwargs):
    r,theta = np.meshgrid(np.linspace(0.5,2,npts),np.linspace(0,np.pi,npts))
    z,y = r*np.cos(theta), r*np.sin(theta)
    Mnm,Nmn = get_MNnm(r,theta,phi, Nmax, zn=spu.hn1,znp=spu.hn1p,pol='xy')
    for l in range(1,Nmax+1):
        idl = int(l*(l+1)/2)
        for m in range(l+1):
            for i in range(3):
                print('l=%d,m=%d,i=%d' %(l,m,i))
                N = np.real(Nmn[idl+m,i])
                cmpt = [r'$e_r$',r'$e_{\theta}$', r'$e_{\phi}$'][i]
                tle = r'$N_{%d}^{%d}$'%(l,m) +' component '+cmpt
                dsp.stddisp(im=[y,z,N],labs=['$y$','$z$'],title=tle,
                    name=path+'M%s%d%d.png' %(['r','t','p'][i],l,m),**kwargs)



plt.close('all')
plot_spherical_Bessel(opt='j')
# plot_real_sph_harm()
# plot_vec_sh_harm(Nmax=1,npts=121,ntp=60,
#     xylims=0.5*np.array([-1,1,-1,1,-1,1]),pOpt='X',
#     opt='p')
# plot_MN(Nmax=1,npts=100,phi=np.pi/2)#,
    # xylims=0.5*np.array([-1,1,-1,1,-1,1]),pOpt='X',
    # opt='p')
