import importlib as imp
from utils import*
from optics.scattering.spherical_utils import*
import scipy.special as spe
import optics.scattering.spherical_utils as spu ;imp.reload(spu)
import utils.displayStandards as dsp            ;imp.reload(dsp)
path='../../docs/figures/'


####################################################################################################################
#                       Plane wave
####################################################################################################################
def plot_plane_wave_scalar(lam=1,d_p=2.5,alpha=0,R=5,npts=50,Nmax=10,fz=np.real,**kwargs):
    r'''plots Plane Wave Ex = exp(jkz) expansion upon Spherical waves.
    - lam : wavelength
    - d_p : distance to origin
    - alpha : incident angle degrees
    - R     : radius extent of the plot
    - nmax  : max order expansion
    - npts  : nb points along z
    '''
    yz0 = 1e-3
    y,z = np.meshgrid( R*np.linspace(yz0,1,npts), R*np.linspace(yz0,1,npts))
    Ex = spu.get_plane_wave_scalar(y,z,lam,d_p,alpha,Nmax,fz)

    a_str=['',r'\cos\alpha'][np.float(alpha)>0]
    tle  = r'$E_x=e^{jkz%s}$, $d_p=$%.1f, $\alpha$=%d$^{\circ}$, R=%1.f, $N_{max}$=%d ' %(a_str,d_p,np.rad2deg(alpha),R,Nmax)
    labs = ['$y$','$z$']
    dsp.stddisp(im=[y,z,Ex],labs=labs,title=tle,
        imOpt='cv',axPos='V',pOpt='pt',#xyTicks=[[],[]]
        **kwargs)

def plot_plane_wave_vector3d(lam=1,d_p=0,alpha=0,R=5,npts=50,Nmax=10,fz=np.real,cmpt=0,pol='y',name='',**kwargs):
    r'''plots Plane Wave expansion upon Vector Spherical waves in 3D
    A quiver plot is done for electric field map
    - cmpt : A scatter is plot for E[compt]
    '''
    cmap,ee = 'jet',1e-3
    fig,ax = dsp.stddisp(rc='3d',opt='')

    x0,y0,z0 = np.meshgrid( R*np.linspace(ee,1,3), R*np.linspace(ee,1,3),R*np.linspace(ee,1,5))
    r,theta,phi = spu.cart2sphere(x0,y0,z0)
    Esph = spu.get_plane_wave_vector(r,theta,phi,lam,alpha,d_p,Nmax,fz,pol=pol,v=True)
    E = spu.Esphere2cart(Esph,theta,phi)
    u = np.zeros(r.shape);u=E[0,:,:,:]
    v = np.zeros(r.shape);v=E[1,:,:,:]
    w = np.zeros(r.shape);w=E[2,:,:,:]
    ax.quiver(x0,y0,z0,u,v,w,color='r',length=0.25,linewidth=3)#,normalize=True)

    mm = np.abs([E.min(),E.max()]).max()
    x,y,z = np.meshgrid([0],R*np.linspace(ee,1,npts), R*np.linspace(ee,1,npts))
    r,theta,phi = spu.cart2sphere(x,y,z)
    Esph = spu.get_plane_wave_vector(r,theta,phi,lam,alpha,d_p,Nmax,fz,pol=pol,v=True)
    E0 = spu.Esphere2cart(Esph,theta,phi)
    ax.scatter(x,y,z,c=E0[cmpt],vmin=-mm,vmax=mm,cmap=cmap)
    #
    x,y,z = np.meshgrid( R*np.linspace(ee,1,npts), [R] ,R*np.linspace(ee,1,npts))
    r,theta,phi = spu.cart2sphere(x,y,z)
    Esph = spu.get_plane_wave_vector(r,theta,phi,lam,alpha,d_p,Nmax,fz,pol=pol,v=True)
    E1 = spu.Esphere2cart(Esph,theta,phi)
    ax.scatter(x,y,z,c=E1[cmpt],vmin=-mm,vmax=mm,cmap=cmap)
    #
    x,y,z = np.meshgrid( R*np.linspace(ee,1,npts), R*np.linspace(ee,1,npts),[0,R/2])#[R/4,3*R/4])
    r,theta,phi = spu.cart2sphere(x,y,z)
    Esph = spu.get_plane_wave_vector(r,theta,phi,lam,alpha,d_p,Nmax,fz,pol=pol,v=True)
    E2 = spu.Esphere2cart(Esph,theta,phi)
    ax.scatter(x,y,z,c=E2[cmpt],vmin=-mm,vmax=mm,cmap=cmap)

    tle = r'scatter plot for $E_{%s}$' %(['x','y','z'][cmpt])
    nameE = 'E%s.png' %(['x','y','z'][cmpt])
    dsp.stddisp(fig=fig,ax=ax,labs=['$x$','$y$','$z$'],title=tle,imOpt='c',caxis=[-mm,mm],cmap=cmap,
        name=name+nameE,**kwargs)


def plot_plane_wave_vector(lam=1,d_p=2.5,alpha=0,R=5,npts=50,Nmax=10,cmpts=[0,1,2],fz=np.real,cart_opt=True,phi=0,pol='y',name='',**kwargs):
    r'''plots Plane Wave Ex = exp(jkz) expansion upon Vector Spherical waves.
    - lam : wavelength
    - d_p : distance to origin
    - alpha : incident angle degrees
    - R     : radius extent of the plot
    - nmax  : max order expansion
    - npts  : nb points along z
    '''
    yz0=1e-3
    y,z = np.meshgrid( R*np.linspace(yz0,1,npts), R*np.linspace(yz0,1,npts))
    r,theta = np.sqrt(y**2+z**2),np.arctan(y/z)

    #get plane wave
    Esph = spu.get_plane_wave_vector(r,theta,phi,lam,alpha,d_p,Nmax,fz,pol=pol,v=True)
    if cart_opt:
        E = spu.Esphere2cart(Esph,theta,phi)
    else:
        E = Esph

    #display
    params = r'$\lambda$=%.1f, $d_p=$%.2f, $\alpha=$%d$^{\circ}$, $N_{max}=%d$ ' %(lam,d_p,np.rad2deg(alpha),Nmax)
    mm = np.abs([E.min(),E.max()]).max()
    for cmpt in cmpts:
        if cart_opt:
            tle  = r'Field $E_{%s}$, ' %['x','y','z'][cmpt]
            nameE = 'E%s.png' %(['x','y','z'][cmpt])

        else:
            tle  = r'Field $E_{%s}$, ' %['r',r'\theta',r'\phi'][cmpt]
            nameE = 'E%s.png' %(['r','t','p'][cmpt])
        tle += params
        dsp.stddisp(im=[y,z,E[cmpt]],labs=['$y$','$z$'],title=tle,
            imOpt='cv',caxis=[-mm,mm],axPos='V',
            name=name+nameE,**kwargs)



#################################################################################################
# Plane wave tests
#################################################################################################
plt.close('all')
pol='y'
opt='p'

plot_plane_wave_scalar(npts=100,Nmax=50 ,lam=1,R=5,alpha=10, opt=opt,name=path+'Exi_sphere1.png')
# plot_plane_wave_scalar(npts=100,nmax=10 ,lam=1,R=5, opt='p',name=path+'Exi_sphere1.png')
# plot_plane_wave_scalar(npts=100,nmax=20 ,lam=1,R=5, opt='p',name=path+'Exi_sphere2.png')
# plot_plane_wave_scalar(npts=100,nmax=50 ,lam=1,R=5, opt='p',name=path+'Exi_sphere3.png')
# plot_plane_wave_scalar(lam=1,alpha=20,d_p=1,R=10,nmax=20,npts=100, opt='p',name=path+'Exi_alpha_sphere1.png')
# plot_plane_wave_scalar(lam=1,alpha=20,d_p=1,R=10,nmax=30,npts=100, opt='sc',name=path+'Exi_alpha_sphere2.png')
# plot_plane_wave_scalar(lam=1,alpha=20,d_p=1,R=10,nmax=40,npts=100, opt='sc',name=path+'Exi_alpha_sphere3.png')
# plot_plane_wave_scalar(lam=1,alpha=20,d_p=1,R=10,nmax=50,npts=100, opt='sc',name=path+'Exi_alpha_sphere4.png')
### vector plane wave cartesian and spherical components
for cart_opt in [True,False]:plot_plane_wave_vector(lam=1,cart_opt=cart_opt,
    R=3,npts=50,Nmax=30,d_p=0,alpha=0, phi=np.pi/2,pol=pol,
    xylims=3*np.array([0,1,0,1]),
    opt=opt,name=path+'pw_vec')
### vector plane wave Ey varying d_p
for i,d_p in enumerate([0,0.25,0.5,0.75]): plot_plane_wave_vector(lam=1,d_p=d_p,cmpts=[1],
    alpha=0,R=2,cart_opt=True,phi=np.pi/2,pol=pol,npts=50,Nmax=20,
    xylims=2*np.array([0,1,0,1]),
    opt=opt,name=path+'pw_vec_d%d' %i)

### vector plane wave 3d plot Ex,Ey,Ez as colorbar
for cmpt in [0,1,2]:plot_plane_wave_vector3d(lam=1,cmpt=cmpt,pol=pol,d_p=0,alpha=0,npts=50,R=1,Nmax=20,
    pOpt='tX',xyTicks=[1,1,1],xylims=1*np.array([0,1,0,1,0,1]),axPos='L',
    opt=opt,name=path+'pw_vec3d')
