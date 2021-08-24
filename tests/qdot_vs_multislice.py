from utils import * ;imp.reload(dsp);imp.reload(cst)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from multislice import multislice as ms;
from multislice import mupy_utils as mut;
from multislice import postprocess as pp;
import scipy.fft as fft

plt.close('all')
path = 'data/qdot/'
opts = 'p' #t(transmission function), p(pattern), r(run)
nz = 10

keV = 50  #keV
Za = 12
r0,V0 = {
    100 :[0.1,4],
    1 :[0.1,0.5],
    2 :[0.2,0.5],
    3 :[0.3,0.5],
    4 :[0.4,0.5],
    5 :[0.5,0.5],
    6 :[0.1,1.0],
    7 :[0.2,1.0],
    8 :[0.3,1.0],
    9 :[0.4,1.0],
    10:[0.1,1.0],
    12:[0.1,2.0],
    13:[0.1,3.0],
    14:[0.1,4.0],
}[Za]

kdka = 5   #
sig,lam = cst.keV2sigma(keV), cst.keV2lam(keV)
E = cst.lam2keV0(lam)/cst.meff(keV) #cst.lam2E(lam)/1e3
E = keV
kp = np.sqrt(1+V0/E)#*np.sqrt(cst.meff(keV))
ka = 2*np.pi/lam*r0#*0.9725
# kp=1.033
# kp = 1.033 #V=0.5
# kp = 1.061  #V=1.0,rr=0.9725,Za=6
# kp = 1.068#V=1.0,rr=0.9725,Za=7
# kp = 1.068#V=1.0,rr=0.9725,Za=8,
kd  = kdka*ka
xyz = path+'qdot.xyz'
Nxy = 2048
nmax = int(np.ceil(ka*1.1))

Npad = 200
#### make xyz
if 'r' in opts:
    a0  = kdka*r0
    lat_params = (Npad*r0,Npad*r0,kdka*r0)
    lat_vec = np.diag(lat_params)
    pattern = np.array([[Za,Npad/2*r0,Npad/2*r0,0, 1,0.05]])
    mut.make_xyz(xyz,pattern,lat_vec,lat_params)
    multi = ms.Multislice(path,keV=keV,data=xyz,NxNy=Nxy,repeat=[1,1,nz],
        dz=a0,i_slice=1,fopt='f');
    multi.p.wait()

multi= pp.load(path,'')

#### check Vz and T are as predicted
if 't' in opts:
    xylims = [Nxy/2-30,Nxy/2+30]*2
    r0s = 50*r0-np.arange(Nxy)*100*r0/Nxy
    x,y = np.meshgrid(r0s,r0s)
    r2 = x**2+y**2
    Vz = np.zeros(r2.shape)
    Vz[r2<=r0**2] = 4*np.pi*V0*np.sqrt(r0**2-r2[r2<=r0**2])

    # Vzms = np.loadtxt(path+'vz')
    # Vzms = Vzms[:,::2]+1J*Vzms[:,1::2]
    # dsp.stddisp(im=[Vzms.real/1e3]            ,caxis=[0,7],xylims=xylims,pOpt='im',title='TEMSIM Re T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[Vz] ,caxis=[0,7],xylims=xylims,pOpt='im',title='py Re T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[Vzms.real/1e3-Vz],caxis=[-1e-5,1e-5],pOpt='im',title='err Re T iz=%d' %iz)#,**kwargs)

    # Tms = np.loadtxt(path+'vz')
    iz=0
    Tms = np.loadtxt(path+'translayer.%s' %str(iz).zfill(3) )
    Tms = Tms[:,::2]+1J*Tms[:,1::2]

    T = np.exp(1J*sig*Vz)

    dsp.stddisp(im=[Tms.real]       ,caxis=[-1,1],xylims=xylims,pOpt='im',title='TEMSIM Re T iz=%d' %iz)#,**kwargs)
    dsp.stddisp(im=[T.real]         ,caxis=[-1,1],xylims=xylims,pOpt='im',title='py Re T iz=%d' %iz)#,**kwargs)
    dsp.stddisp(im=[Tms.real-T.real],caxis=[-1,1],xylims=xylims,pOpt='im',title='err Re T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[Tms.imag],pOpt='im',title='TEMSIM Im T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[T.imag],pOpt='im',title='py Im T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[Tms.imag-T.imag],pOpt='im',title='err Im T iz=%d' %iz)#,**kwargs)
    # dsp.stddisp(im=[abs(Tms-T)],pOpt='im',caxis=[0,1],title='abs(err T) iz=%d' %iz)#,**kwargs)


if 'p' in opts:
    optsP = '' #a(analytical),f(fft)
    nz=10
    Nmax = 500
    s = np.s_[0,:Nmax]
    # qx,qy,Ip = multi.pattern(iz=0,Iopt='',out=1)
    # qx0 = qx[s]
    ax,by = multi.cell_params[:2]
    # Nh,Nk = multi.repeat[:2];
    dq = 1/(ax)
    qx0 = np.arange(Nmax)*dq
    thetas = 2*np.arcsin(qx0*lam/2)
    thetas_d=np.rad2deg(thetas)
    plts,cs,ms = [],dsp.getCs('jet',nz),'x'
    for iz in range(nz):
        # im = np.loadtxt(path+'qdot_autoslic_pattern.txt.%s' %str(iz).zfill(3))
        # real,imag = im[:,0:-1:2],im[:,1::2];#print(real.max())
        # I = real**2+imag**2
        I = np.load(path+'qdot_autoslic_pattern%s.npy' %str(iz).zfill(3))
        # qx,qy,Ip = multi.pattern(iz=iz,Iopt='',out=1)
        I0 = I[s]
        plts += [[thetas_d,I0/I0[1],[cs[iz],'-'+ms],'']]

    qdot = [qsa.QdotSphereArray(N=iz,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,
        copt=iz>1) for iz in np.arange(nz)+1 ]
    Iff = lambda ff:(abs(ff)/abs(ff).max())**2
    for iz in range(nz):
        plts += [[thetas_d,Iff(qdot[iz].get_ff(thetas)),[cs[iz],'--'+ms],'']]

    if 'a' in optsP:
        from scipy.integrate import trapz,quad
        sig0 = sig
        sig0 = np.pi/(lam*E)
        vz = lambda r2:2*V0*np.sqrt(r0**2-r2)
        fvr = lambda rho,q :np.real((np.exp(1J*sig0*vz(rho**2))-1)*spe.jv(0,2*np.pi*q*rho)*rho)
        fvi = lambda rho,q :np.imag((np.exp(1J*sig0*vz(rho**2))-1)*spe.jv(0,2*np.pi*q*rho)*rho)
        fqr = 2*np.pi*1J/lam*np.array([quad(fvr,0,r0,(q0,))[0] for q0 in qx0])
        fqi = 2*np.pi*1J/lam*np.array([quad(fvi,0,r0,(q0,))[0] for q0 in qx0])
        fq = fqr+1J*fqi
        plts += [[thetas_d,Iff(fq),'k:','analytical']]

    if 'f' in optsP:
        Tms = np.loadtxt(path+'translayer.%s' %str(0).zfill(3) )
        Tms = Tms[:,::2]+1J*Tms[:,1::2]
        fq2 = abs(fft.fft2(Tms)[s])**2;I=fq2/fq2[1]
        plts += [[thetas_d,I,'k-.','fft2']]

    legElt = {'MS':'k-','pyscat':'k--'}
    legElt.update({'iz=%d' %iz: [cs[iz],ms] for iz in range(nz)})
    dsp.stddisp(plts,labs=[r'$\theta$(degrees)','I'],legElt=legElt,lw=2,
        xylims=[0,thetas_d.max(),0,0.02])

    # dsp.stddisp(im=[I],caxis=[0,10000],pOpt='im')
    # multi.pattern(iz=0,Iopt='s',caxis=[0,10000],cmap='jet')
