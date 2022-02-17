from utils import *                                ;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from scipy.integrate import trapz,quad
import scipy.special as spe
plt.close('all')
name='figures/TmatrixSingleF.eps'

#### parameter sets to show scattering amplitudes
ka1,kp1 = 3,1.1
ka2,kp2 = 10,1.02
# lam = cst.keV2lam(200)
k0,K = 1,1/(2*np.pi)#2*np.pi/lam,1/lam
npts=180*8

#### Born and MS approximations
thetas = np.linspace(1e-3,np.pi,npts)
theta  = np.rad2deg(thetas)
q      = 2*K*np.sin(thetas/2)
plts = []
born = lambda qa:(np.sin(qa)/(qa)-np.cos(qa))/qa**2
def get_MS(q,ka,kp):
    eps = kp**2-1
    sigV0 = k0*eps/2
    fvr = lambda rho,q :np.real((1-np.exp(2J*sigV0*np.sqrt(ka**2-rho**2)))*spe.jv(0,2*np.pi*q*rho)*rho)
    fvi = lambda rho,q :np.imag((1-np.exp(2J*sigV0*np.sqrt(ka**2-rho**2)))*spe.jv(0,2*np.pi*q*rho)*rho)
    fqr = np.array([quad(fvr,0,ka,(q0,))[0] for q0 in q])
    fqi = np.array([quad(fvi,0,ka,(q0,))[0] for q0 in q])
    fq = fqr+1J*fqi
    return fq
ff = lambda f:np.abs(f)/np.abs(f).max()


#### first set
qdot1 = qsa.QdotSphereArray(N=1,ka=ka1,kp=kp1,kd=0,copt=0,solve=1,nmax=int(ka1)+3)
f1    = ff(qdot1.get_ff(thetas))
f1b   = ff(born(ka1*2*np.pi*q))
f1m   = ff(get_MS(q,ka1,kp1))
plts += [[theta,f1 ,'g-','$ka=%.1f$' %ka1]]
plts += [[theta,f1b,'g--','']]
plts += [[theta,f1m,'g-.','']]

#### second set
qdot2 = qsa.QdotSphereArray(N=1,ka=ka2,kp=kp2,kd=0,copt=0,solve=1,nmax=int(ka2)+3)
f2    = ff(qdot2.get_ff(thetas))
f2b   = ff(born(ka2*2*np.pi*q))
f2m   = ff(get_MS(q,ka2,kp2))
plts += [[theta,f2 ,'b-','$ka=%d$' %ka2]]
plts += [[theta,f2b,'b--','']]
plts += [[theta,f2m,'b-.','']]


legElt = {'T-matrix':'k-','Born':'k--','MS':'k-.'}
dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],
    legElt=legElt,
    lw=2,fonts = '2',xylims=[0,180,-1e-2,1.01],
    name=name,opt='p')
